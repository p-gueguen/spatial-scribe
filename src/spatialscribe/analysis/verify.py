"""Self-verification of cell-type labels against canonical markers (advisory).

What it does
------------
Answers "do the labels actually agree with the markers?" - the self-check half of the project's
"self-verify, then re-run what fails" vision. For each assigned cell type it asks whether that
type's cells score highest on THAT type's canonical markers (marker-argmax agreement), corroborates
with a one-vs-rest AUC, classifies WHY a type fails (panel gap / mislabel / hidden substructure),
and proposes concrete, engine-grounded fixes (subcluster / abstain / merge a confusable pair).

It reuses the existing engine and adds NO new marker scoring: the panel-restricted per-lineage
scores come from :func:`annotate.score_marker_sets`, the AUC from
:func:`eval_metrics.marker_program_fidelity`, the panel-adequacy context from ``panel_check``, and
the optional spatial corroborator from :func:`spatial.neighborhood_sanity`.

**It never mutates labels** - it writes ``adata.uns['annotation_verification']`` and returns advice.

How to use it
-------------
>>> res = verify_annotation(adata, marker_sets=ctx.markers())
>>> res["failed"]                       # worst-first list of types that disagree with their markers
>>> suggest_reruns(res)                 # ordered, advisory remediation (dry-run by default)

Depends on
----------
numpy; :mod:`annotate`, :mod:`markers`, :mod:`eval_metrics`, :mod:`panel_check`, :mod:`config`,
and (only when ``neighborhood=True``) :mod:`spatial`. Heavy imports stay inside functions.
"""

from __future__ import annotations

# Defaults; the live values are read from docs/research/annotation_qc_thresholds.yaml via config.get.
ARGMAX_AGREE_MIN = 0.5         # < -> a type's cells mostly do NOT top their own program
AUC_MIN = 0.6                  # one-vs-rest AUC < -> poor marker separation
MIN_CELLS = 20                 # skip types too small to audit reliably
CONFUSION_DOMINANT_FRAC = 0.5  # a single confuser >= this share of the miss -> mislabel/merge
NEAR_CHANCE_AUC = 0.55         # AUC at/below this is ~chance -> low signal / depth, not resolvable substructure

# Major cell compartments, used ONLY to refuse a biologically-unsound cross-compartment merge (never to
# assign identity). Merging within a compartment (e.g. two epithelial sublineages) is defensible
# coarsening; merging an immune granulocyte INTO tumour epithelium or fibroblasts is not - a panel that
# cannot separate them is a reason to abstain, not to dissolve one lineage into another.
_COMPARTMENT: dict[str, str] = {
    # immune / haematopoietic
    "T cell": "immune", "NK cell": "immune", "B/Plasma": "immune", "Myeloid": "immune",
    "Mast": "immune", "Microglia": "immune",
    # epithelial / tumour
    "Epithelial/Tumor": "epithelial", "Basal/Myoepithelial": "epithelial", "Keratinocyte": "epithelial",
    "Malignant/Melanocyte": "epithelial",
    # stromal / mesenchymal
    "Stromal/CAF": "stromal", "Adipocyte": "stromal", "Mural/Pericyte": "stromal",
    # vascular
    "Endothelial": "endothelial",
    # neural / glial
    "Excitatory neuron": "neural", "Inhibitory neuron": "neural", "Astrocyte": "neural",
    "Oligodendrocyte": "neural", "OPC": "neural", "Ependymal": "neural", "Glial/Neural": "neural",
}


def _same_compartment(a: str, b: str) -> bool:
    """True unless a and b are in DIFFERENT major lineages (so a merge across them is refused).

    Prefers the **Cell Ontology** (CellGuide `cl-basic.obo`): two types may merge only if they share a
    major-lineage anchor (leukocyte / epithelial / stromal / endothelial / neural / muscle) - so CD4
    and CD8 T cells coarsen to "T cell", but a mast cell (immune) is never merged into tumour
    epithelium or a fibroblast. This generalises to arbitrary fine subtypes for free (a CD8 effector
    memory T cell still resolves to the leukocyte anchor). Falls back to the hardcoded
    :data:`_COMPARTMENT` map only when the ontology cannot decide (offline, or a label with no CL id)."""
    try:
        from . import cellguide as _cg
        v = _cg.labels_mergeable(str(a), str(b)).get("mergeable")
        if v is not None:
            return bool(v)
    except Exception:
        pass
    ca, cb = _COMPARTMENT.get(str(a)), _COMPARTMENT.get(str(b))
    return not (ca is not None and cb is not None and ca != cb)


def _cfg(key: str, default):
    from . import config
    return config.get("verify", key, default=default)


def verify_annotation(adata, marker_sets: dict[str, list[str]] | None = None,
                      cluster_key: str = "cell_type", neighborhood: bool = False,
                      min_cells: int | None = None) -> dict:
    """Per-cell-type marker-agreement audit of ``adata.obs[cluster_key]`` (writes uns, no mutation).

    For each type: ``argmax_agreement`` (fraction of its cells whose top marker-program IS that
    type), whether its own program is the section-mean top, a one-vs-rest AUC, and - for failures -
    the dominant confuser and a cause (``panel_gap`` / ``mislabel`` / ``heterogeneous``). A type
    with no on-panel markers is ``unscoreable`` (honest, not a failure). Returns the full result.
    """
    import numpy as np

    from . import annotate as _an
    from . import eval_metrics as _em
    from . import markers as _m
    from . import panel_check as _pc

    marker_sets = marker_sets or _m.LINEAGE_MARKERS
    argmax_min = float(_cfg("argmax_agreement_min", ARGMAX_AGREE_MIN))
    auc_min = float(_cfg("auc_min", AUC_MIN))
    min_cells = int(min_cells if min_cells is not None else _cfg("min_cells", MIN_CELLS))
    conf_dom = float(_cfg("confusion_dominant_frac", CONFUSION_DOMINANT_FRAC))
    thresholds = {"argmax_agreement_min": argmax_min, "auc_min": auc_min,
                  "min_cells": min_cells, "confusion_dominant_frac": conf_dom}

    result: dict = {"cluster_key": cluster_key, "per_type": {}, "failed": [],
                    "section_agreement": None, "mean_auc": None, "n_types_scored": 0,
                    "thresholds": thresholds}

    colmap = _an.score_marker_sets(adata, marker_sets)   # {lineage: score_col} - REUSE, no new scoring
    names = list(colmap)
    labels = adata.obs[cluster_key].astype(str).to_numpy()
    if not names:
        adata.uns["annotation_verification"] = result
        return result

    S = adata.obs[[colmap[nm] for nm in names]].to_numpy(dtype=float)
    names_arr = np.array(names)
    argmax_type = names_arr[S.argmax(1)]

    # AUC corroborator (reuse) - restrict marker sets to on-panel genes.
    present_sets = _m.present(set(adata.var_names), marker_sets)
    mpf = _em.marker_program_fidelity(adata, cluster_key, present_sets)
    mpf_per = mpf.get("per_type", {})

    # panel-adequacy context (guarded against an h5ad round-trip that corrupts panel_check).
    pc = adata.uns.get("panel_check")
    pc = pc if _pc.is_valid(pc) else None
    coverage = pc.get("coverage", {}) if pc else {}
    # Map each type to the types it is ACTUALLY confusable with on this panel (a genuine confusable
    # pair, i.e. neither has a private on-panel marker). Only a real pair justifies a merge; a low
    # marker-coverage COUNT alone does not (that is a depth/sparse-marker problem, not a confusion).
    confusable_partners: dict[str, set[str]] = {}
    if pc:
        for p in pc.get("confusable_pairs", []):
            pair = p.get("pair", []) if isinstance(p, dict) else []
            if len(pair) == 2:
                a0, b0 = str(pair[0]), str(pair[1])
                confusable_partners.setdefault(a0, set()).add(b0)
                confusable_partners.setdefault(b0, set()).add(a0)

    per: dict[str, dict] = {}
    scored_agreements: list[float] = []
    for t in list(dict.fromkeys(labels.tolist())):
        mask = labels == t
        ncell = int(mask.sum())
        base = {"n_cells": ncell, "argmax_agreement": None, "self_program_is_top": None,
                "self_rank": None, "auc": None, "cohens_d": None, "confused_with": None,
                "confusion_frac": 0.0, "cause": None}
        if t not in colmap:
            per[t] = {**base, "status": "unscoreable"}     # panel has no markers for this label
            continue
        if ncell < min_cells:
            per[t] = {**base, "status": "skipped_small"}
            continue
        am = argmax_type[mask]
        agr = float(np.mean(am == t))
        scored_agreements.append(agr)
        miss = am[am != t]
        confused_with, conf_frac = None, 0.0
        if len(miss):
            vals, counts = np.unique(miss, return_counts=True)
            k = int(counts.argmax())
            # Share OF THE MISS (not of all the type's cells): confusion_dominant_frac gates on how
            # much of the *disagreement* one sibling owns, so a mostly-correct type with a single
            # dominant confuser among its few misses is still classified mislabel, not heterogeneous.
            confused_with, conf_frac = str(vals[k]), float(counts[k] / len(miss))
        row_means = S[mask].mean(0)
        order = np.argsort(-row_means)
        self_is_top = bool(names_arr[order[0]] == t)
        self_rank = int(np.where(names_arr[order] == t)[0][0]) + 1
        aucd = mpf_per.get(t, {})
        auc = aucd.get("auc")
        cohens_d = aucd.get("cohens_d")
        near_chance = float(_cfg("near_chance_auc", NEAR_CHANCE_AUC))
        failed = (agr < argmax_min) and ((not self_is_top) or (auc is not None and auc < auc_min))
        cause, merge_target = None, None
        if failed:
            partners = confusable_partners.get(t, set())
            cov_gap = isinstance(coverage.get(t), dict) and coverage[t].get("status") in ("red", "amber")
            low_signal = (auc is None) or (auc <= near_chance)
            if partners:
                # a REAL panel confusable pair -> coarsen, but NEVER across major compartments (do not
                # dissolve an immune type into tumour/stroma); a cross-compartment-only confusion abstains.
                same = sorted(p for p in partners if _same_compartment(t, p))
                if same:
                    cause, merge_target = "panel_gap", same[0]
                else:
                    cause = "low_signal"
            elif (not self_is_top) and conf_frac >= conf_dom and confused_with is not None:
                # cells express ANOTHER program (own program is NOT the top) with one dominant confuser.
                cause = "mislabel"
            elif cov_gap or low_signal:
                # too few on-panel markers, OR AUC ~chance: depth/panel-limited, not a mislabel or a
                # resolvable substructure -> abstain the low-signal cells, never a confident label/merge.
                cause = "low_signal"
            else:
                # own program tops, adequate markers AND signal, cells still split -> real substructure.
                cause = "heterogeneous"
        per[t] = {"n_cells": ncell, "argmax_agreement": round(agr, 3),
                  "self_program_is_top": self_is_top, "self_rank": self_rank,
                  "auc": (float(auc) if auc is not None else None),
                  "cohens_d": (float(cohens_d) if cohens_d is not None else None),
                  "confused_with": confused_with, "confusion_frac": round(conf_frac, 3),
                  "merge_target": merge_target,
                  "status": "fail" if failed else "pass", "cause": cause}

    # Optional spatial corroborator (heavy squidpy import stays inside the branch).
    if neighborhood and "spatial" in getattr(adata, "obsm", {}):
        try:
            from . import spatial as _sp
            ns = _sp.neighborhood_sanity(adata, cluster_key=cluster_key)
            for s in ns.get("suspicious", []):
                t = str(s.get("type"))
                if t in per:
                    per[t]["spatially_incoherent"] = True
        except Exception:
            pass

    failed_list = sorted(
        [{"cell_type": t, "argmax_agreement": per[t]["argmax_agreement"], "auc": per[t]["auc"],
          "confused_with": per[t]["confused_with"], "merge_target": per[t].get("merge_target"),
          "cause": per[t]["cause"], "n_cells": per[t]["n_cells"]}
         for t in per if per[t]["status"] == "fail"],
        key=lambda r: (r["argmax_agreement"] if r["argmax_agreement"] is not None else 0.0))

    result["per_type"] = per
    result["failed"] = failed_list
    result["section_agreement"] = float(np.mean(scored_agreements)) if scored_agreements else None
    result["mean_auc"] = mpf.get("mean_auc")
    result["n_types_scored"] = len(scored_agreements)
    adata.uns["annotation_verification"] = result
    return result


def suggest_reruns(verify_result: dict) -> list[dict]:
    """Map each failed type to an ordered, engine-grounded remediation PLAN (advisory - it never
    executes anything and never mutates labels).

    Causes -> actions: ``panel_gap`` -> merge the confusable pair (advisory, no capability);
    ``heterogeneous`` -> ``subcluster`` (safe rerun, writes obs['subtype'], never overwrites
    cell_type); ``mislabel`` -> re-run ``annotate`` so ``apply_confidence`` abstains those cells,
    then an advisory relabel-to-``confused_with``. Each item carries the ``capability`` + ``params``
    to run; ACTUALLY running a chosen rerun is the UI/adapter layer's job (it calls
    ``capabilities.run`` on a button click), keeping this pure function plan-only.
    """
    order_rank = {"panel_gap": 0, "heterogeneous": 1, "low_signal": 2, "mislabel": 3}
    items: list[dict] = []
    for f in verify_result.get("failed", []):
        t, cause, cw = f["cell_type"], f.get("cause"), f.get("confused_with")
        mt = f.get("merge_target")
        if cause == "panel_gap":
            # merge ONLY toward a genuine, same-compartment confusable partner (from panel_check's
            # confusable_pairs) - never toward the marker-argmax residual, and never across compartments.
            onto = None
            try:
                from . import cellguide as _cg
                onto = _cg.labels_mergeable(str(t), str(mt))
            except Exception:
                onto = None
            ground = ""
            if onto and onto.get("shared_compartments"):
                ground = (f" Both resolve to the same {'/'.join(onto['shared_compartments'])} lineage in "
                          f"the Cell Ontology ({onto.get('cl_a')} / {onto.get('cl_b')}), so coarsening them "
                          f"is biologically defensible.")
            items.append({"cell_type": t, "cause": cause, "action": "merge", "capability": None,
                          "params": None, "with": mt, "advisory": True, "ontology": onto,
                          "reason": (f"The panel has no private on-panel marker separating '{t}' from "
                                     f"'{mt}' (a confusable pair). Consider merging them into one coarser "
                                     f"label.{ground}"),
                          "code": f"adata.obs['cell_type'] = adata.obs['cell_type'].replace({{'{t}': '{t}+{mt}'}})"})
        elif cause == "low_signal":
            # too few on-panel markers, or AUC ~chance (depth-limited): abstain the disagreeing cells -
            # honest, never a confident label, and never a cross-lineage merge into the dominant compartment.
            items.append({"cell_type": t, "cause": cause, "action": "abstain", "capability": None,
                          "params": None, "advisory": False,
                          "reason": (f"'{t}' cannot be confidently separated on this panel at this depth "
                                     f"(few markers or near-chance AUC). Abstain the low-signal cells rather "
                                     f"than forcing a confident label or merging into another lineage."),
                          "code": f"verify._abstain_mislabelled(adata, cluster_key, '{t}', markers)"})
        elif cause == "heterogeneous":
            items.append({"cell_type": t, "cause": cause, "action": "subcluster",
                          "capability": "subcluster", "params": {"cell_type": t}, "advisory": False,
                          "reason": (f"'{t}' cells do not agree on one marker program (no single "
                                     f"confuser); likely hidden substructure. Subcluster it to resolve."),
                          "code": f"cap.run(adata, 'subcluster', {{'cell_type': '{t}'}}, ctx)"})
        elif cause == "mislabel":
            items.append({"cell_type": t, "cause": cause, "action": "abstain",
                          "capability": None, "params": None, "advisory": False,
                          "reason": (f"Most '{t}' cells score highest on '{cw}' markers. Abstain the "
                                     f"disagreeing cells (honest - never a confident wrong label)."),
                          "code": f"verify._abstain_mislabelled(adata, cluster_key, '{t}', markers)"})
            items.append({"cell_type": t, "cause": cause, "action": "relabel", "capability": None,
                          "params": None, "with": cw, "advisory": True,
                          "reason": f"If '{t}' is a mislabel, relabel these cells to '{cw}' after review.",
                          "code": f"adata.obs['cell_type'] = adata.obs['cell_type'].replace({{'{t}': '{cw}'}})"})
        else:
            items.append({"cell_type": t, "cause": cause, "action": "review", "capability": None,
                          "params": None, "advisory": True,
                          "reason": f"'{t}' failed marker agreement for an unclassified reason; review.",
                          "code": ""})
    items.sort(key=lambda it: order_rank.get(it.get("cause") or "", 3))
    return items


# --------------------------------------------------------------------------- #
# The executor: self-verify, then RE-RUN what fails (the loop half of the vision).
# --------------------------------------------------------------------------- #
# SAFE fixes run automatically: subclustering a heterogeneous type (writes obs['subtype'], never
# touches cell_type), and ABSTAINING the disagreeing cells of a mislabelled type (honest - it sets
# them to an abstention label so the type can LEAVE the failing set, never a new confident wrong
# lineage). Still ADVISORY (never auto-applied): merging a confusable pair, or relabelling to the
# confuser - a dominant sibling program is not proof of identity.
_AUTO_ACTIONS = {"subcluster", "abstain"}


def _abstain_mislabelled(adata, cluster_key: str, cell_type: str,
                         marker_sets: dict, label: str | None = None) -> int:
    """Abstain the cells of ``cell_type`` whose top marker-program is NOT ``cell_type`` (the
    mislabelled subset) by setting them to an abstention label. Honest by construction: it never
    relabels to the confuser. Reuses :func:`annotate.score_marker_sets` (no new scoring). Returns
    the number abstained.
    """
    import numpy as np
    import pandas as pd

    from . import annotate as _an

    label = label or _an.NOT_ASSIGNED
    colmap = _an.score_marker_sets(adata, marker_sets)
    names = list(colmap)
    if cell_type not in colmap or not names:
        return 0
    scores = adata.obs[[colmap[nm] for nm in names]].to_numpy(dtype=float)
    argmax = np.asarray(names)[scores.argmax(1)]
    labels = adata.obs[cluster_key].astype(str).to_numpy()
    mask = (labels == cell_type) & (argmax != cell_type)
    if not mask.any():
        return 0
    s = adata.obs[cluster_key].astype(str)
    s[mask] = label
    adata.obs[cluster_key] = pd.Categorical(s)
    return int(mask.sum())


def _union_markers_for_merge(marker_sets: dict, applied: list[dict]) -> dict:
    """After-merge marker sets: drop each group's per-member programs and add ONE union program
    keyed by the joint label, so the after verify SCORES the merged type (a cell that was A or B
    now argmaxes to 'A / B') instead of silently dropping it as unscoreable. Measurement-only."""
    ms = dict(marker_sets)
    for a in applied:
        union: list[str] = []
        seen: set[str] = set()
        for m in a["group"]:
            for g in marker_sets.get(m, []):
                if g not in seen:
                    seen.add(g)
                    union.append(g)
        for m in a["group"]:
            ms.pop(m, None)
        ms[a["label"]] = union
    return ms


def _lineage_partition(members: list[str]) -> list[list[str]]:
    """Split a panel-confusable GROUP into lineage-coherent subgroups: two members join only when
    :func:`_same_compartment` allows it (Cell-Ontology-grounded, hardcoded-map fallback offline).

    This is the guard against the live "optimize merging" bug - `panel_check.merge_groups` connects any
    two types with no private on-panel marker, so on a targeted panel it blobbed biologically distant
    types together (Oligodendrocyte + B cell, Fibroblast + pigmented epithelial cell). Members that no
    True edge connects come back as singletons (not merged). Unknown labels default to allowed (so the
    historic behaviour holds where the ontology cannot decide)."""
    parent = {m: m for m in members}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, a in enumerate(members):
        for b in members[i + 1:]:
            if _same_compartment(a, b):
                parent[find(a)] = find(b)
    comps: dict[str, list[str]] = {}
    for m in members:
        comps.setdefault(find(m), []).append(m)
    return list(comps.values())


def merge_confusable_types(adata, groups, cluster_key: str = "cell_type",
                           marker_sets: dict | None = None, dry_run: bool = False) -> dict:
    """Collapse each group of panel-confusable types into ONE joint label in ``obs[cluster_key]``,
    and report the annotation-quality delta (section marker-agreement + number of failing types)
    before vs after.

    Honest coarsening, not a guess: ``panel_check`` groups these types precisely because the panel
    has NO private on-panel marker to separate them, so one coarser label is more defensible than a
    coin-flip. Only members actually present in the column are merged; a group with <2 present
    members is skipped. The after verify scores each group's union program (see
    :func:`_union_markers_for_merge`) so the merged type is credited, not dropped. Never raises.

    ``dry_run`` returns the ontology-partitioned merges that WOULD be applied (``would_merge``)
    without touching ``obs`` and without the expensive before/after verify - the "preview before
    apply" path for the Curate-tab suggest button.
    """
    import pandas as pd

    labels = adata.obs[cluster_key].astype(str)
    present = set(labels.unique())
    applied: list[dict] = []
    remap: dict[str, str] = {}
    for g in groups:
        members = [m for m in dict.fromkeys(str(x) for x in g) if m in present]   # stable, present-only
        if len(members) < 2:
            continue
        # Lineage guard: never merge across major cell lineages. A panel_check group can span
        # compartments (no private on-panel marker != same lineage); split it into lineage-coherent
        # subgroups and merge each, so an immune type is never folded into neural/epithelial/stromal.
        for sub in _lineage_partition(members):
            if len(sub) < 2:
                continue
            joint = " / ".join(sub)
            for m in sub:
                remap[m] = joint
            applied.append({"group": sub, "label": joint,
                            "n_cells": int(labels.isin(sub).sum())})
    if dry_run:
        return {"status": "preview", "cluster_key": cluster_key,
                "would_merge": applied, "n_groups": len(applied)}
    before = verify_annotation(adata, marker_sets=marker_sets, cluster_key=cluster_key)
    if remap:
        adata.obs[cluster_key] = pd.Categorical(labels.replace(remap))
    ms_after = _union_markers_for_merge(marker_sets, applied) if (marker_sets and applied) else marker_sets
    after = verify_annotation(adata, marker_sets=ms_after, cluster_key=cluster_key)

    def _q(v):
        return {"section_agreement": v.get("section_agreement"), "mean_auc": v.get("mean_auc"),
                "n_failed": len(v.get("failed", []))}

    qb, qa = _q(before), _q(after)
    delta = (None if qb["section_agreement"] is None or qa["section_agreement"] is None
             else round(qa["section_agreement"] - qb["section_agreement"], 3))
    return {"status": "ok", "cluster_key": cluster_key, "merged": applied,
            "n_groups_merged": len(applied), "quality_before": qb, "quality_after": qa,
            "agreement_delta": delta, "n_failed_delta": qa["n_failed"] - qb["n_failed"]}


def _panel_gap_groups(adata, marker_sets: dict, cluster_key: str) -> list[list[str]]:
    """Groups to merge for a panel-gap remediation: prefer the cached ``panel_check`` merge_groups
    (full connected components - handles triplets like Stromal 1/2/3), else derive pairs from the
    current verify's ``panel_gap`` failures."""
    from . import panel_check as _pc

    pc = adata.uns.get("panel_check")
    if not _pc.is_valid(pc):
        try:
            pc = _pc.check_panel(list(adata.var_names), marker_sets=marker_sets)
            adata.uns["panel_check"] = pc
        except Exception:
            pc = None
    if _pc.is_valid(pc) and pc.get("merge_groups"):
        return [list(g) for g in pc["merge_groups"]]
    v0 = verify_annotation(adata, marker_sets=marker_sets, cluster_key=cluster_key)
    return [[f["cell_type"], f["merge_target"]] for f in v0.get("failed", [])
            if f.get("cause") == "panel_gap" and f.get("merge_target")]


def autorerun(adata, ctx, max_rounds: int = 2, marker_sets: dict | None = None,
              cluster_key: str = "cell_type", merge_confusable: bool = False) -> dict:
    """Self-verify -> auto-run the safe fix for each failing type -> re-verify, up to ``max_rounds``.

    Runs :func:`verify_annotation`, then for every failing type whose safe remediation is a
    ``subcluster`` (a ``heterogeneous`` cause), calls ``capabilities.run('subcluster', ...)`` to
    resolve its hidden substructure, and re-verifies. Each (type, action) is attempted at most once
    (so it converges - subclustering does not change ``cell_type``, so a type would otherwise re-flag
    forever). Advisory fixes (merge / relabel / re-annotate) are collected, never executed. Returns
    the per-round action log, the initial vs final failing types, and the advisory remainder. Never
    raises (a failed sub-run is recorded with its error and the loop continues).
    """
    from . import capabilities as _cap

    ms = marker_sets or ctx.markers(adata)
    # Optional panel-gap remediation FIRST: collapse the types the panel cannot separate into one
    # coarser label (the honest fix for a panel gap), so the subcluster/abstain loop below then runs
    # on the merged labels. Applied only when asked (default keeps the historic advisory behaviour).
    merge_result = None
    if merge_confusable:
        groups = _panel_gap_groups(adata, ms, cluster_key)
        if groups:
            merge_result = merge_confusable_types(adata, groups, cluster_key=cluster_key, marker_sets=ms)
    attempted: set[tuple] = set()
    v = verify_annotation(adata, marker_sets=ms, cluster_key=cluster_key)
    initial_failed = [f["cell_type"] for f in v.get("failed", [])]
    rounds: list[dict] = []
    for r in range(max(1, int(max_rounds))):
        plan = suggest_reruns(v)
        todo = [it for it in plan if it.get("action") in _AUTO_ACTIONS
                and (it["cell_type"], it["action"]) not in attempted]
        if not todo:
            break
        actions = []
        for it in todo:
            attempted.add((it["cell_type"], it["action"]))
            if it["action"] == "abstain":                 # direct, honest mutation (no capability)
                n = _abstain_mislabelled(adata, cluster_key, it["cell_type"], ms)
                actions.append({"cell_type": it["cell_type"], "action": "abstain", "capability": None,
                                "ok": n > 0, "n_abstained": n, "error": None})
            else:                                          # a safe capability rerun (subcluster)
                res = _cap.run(adata, it["capability"], it.get("params") or {}, ctx)
                actions.append({"cell_type": it["cell_type"], "action": it["action"],
                                "capability": it["capability"], "ok": bool(res.ok),
                                "error": (res.error.get("message") if res.error else None)})
        v = verify_annotation(adata, marker_sets=ms, cluster_key=cluster_key)
        rounds.append({"round": r + 1, "actions": actions,
                       "failed_after": [f["cell_type"] for f in v.get("failed", [])]})
    final_failed = [f["cell_type"] for f in v.get("failed", [])]
    advisory = [{"cell_type": it["cell_type"], "action": it["action"], "reason": it["reason"]}
                for it in suggest_reruns(v) if it.get("advisory")]
    return {"status": "ok", "max_rounds": int(max_rounds),
            "initial_failed": initial_failed, "final_failed": final_failed,
            "n_auto_actions": sum(len(rr["actions"]) for rr in rounds),
            "n_fixed": len(set(initial_failed) - set(final_failed)),
            "rounds": rounds, "advisory_remaining": advisory, "merge": merge_result}
