"""Self-verification of annotation (analysis/verify.py) - synthetic, deterministic, no network.

A deliberately mislabelled type must be flagged with the right confuser and cause; clean labels
must pass; suggest_reruns must be ordered/advisory and mutate nothing; the self_verify capability
must produce uns['annotation_verification'] and be JSON-able.
"""
from __future__ import annotations

import json

import anndata
import numpy as np
import pandas as pd


def _labelled(mislabel: bool = False, mystery: bool = False):
    """3 disjoint lineages (T cell / Myeloid / Endothelial). If mislabel, most 'T cell'-labelled
    cells actually express Myeloid markers (so 'T cell' should fail, confused_with 'Myeloid')."""
    import scanpy as sc

    from spatialscribe.analysis import markers as m

    markers = list(dict.fromkeys(
        m.LINEAGE_MARKERS["T cell"] + m.LINEAGE_MARKERS["Myeloid"] + m.LINEAGE_MARKERS["Endothelial"]))
    # Background genes so scanpy score_genes has a control pool (real sections have hundreds of
    # genes). They must SPAN the expression range - each gets its own rate 0.5..20 - so every
    # expression bin has non-marker controls (a uniformly-low background leaves the high bins
    # marker-only and score_genes raises "No control genes found").
    n_bg = 80
    genes = markers + [f"BG{i}" for i in range(n_bg)]
    gi = {g: i for i, g in enumerate(genes)}
    rng = np.random.default_rng(0)
    bg_rates = rng.uniform(0.5, 20.0, size=n_bg)

    def block(n, lineage):
        M = rng.poisson(0.3, size=(n, len(genes))).astype("float32")
        M[:, len(markers):] += rng.poisson(np.broadcast_to(bg_rates, (n, n_bg))).astype("float32")
        for g in m.LINEAGE_MARKERS[lineage]:
            if g in gi:
                M[:, gi[g]] += rng.poisson(30, size=n)
        return M

    if mislabel:
        specs = [(40, "T cell", "T cell"), (80, "Myeloid", "T cell"), (60, "Endothelial", "Endothelial")]
    else:
        specs = [(60, "T cell", "T cell"), (60, "Myeloid", "Myeloid"), (60, "Endothelial", "Endothelial")]
    blocks, labels = [], []
    for n, express, lab in specs:
        blocks.append(block(n, express))
        labels += [lab] * n
    if mystery:
        blocks.append(block(30, "T cell"))            # some cells with a label that has no markers
        labels += ["MysteryType"] * 30
    a = anndata.AnnData(X=np.vstack(blocks))
    a.var_names = genes
    a.obs["cell_type"] = pd.Categorical(labels)
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    return a


def test_verify_flags_mislabelled_type():
    from spatialscribe.analysis import verify

    a = _labelled(mislabel=True)
    res = verify.verify_annotation(a)
    failed = [f["cell_type"] for f in res["failed"]]
    assert "T cell" in failed
    assert res["per_type"]["T cell"]["argmax_agreement"] < 0.5
    assert res["per_type"]["T cell"]["confused_with"] == "Myeloid"
    assert res["per_type"]["Endothelial"]["status"] == "pass"
    assert "Endothelial" not in failed
    json.dumps(res)                                    # fully JSON-able


def test_verify_clean_passes():
    from spatialscribe.analysis import verify

    a = _labelled(mislabel=False)
    res = verify.verify_annotation(a)
    assert res["failed"] == []
    assert res["section_agreement"] > 0.9
    assert all(d["status"] == "pass" for d in res["per_type"].values() if d["status"] != "unscoreable")


def test_unscoreable_type_not_failed():
    from spatialscribe.analysis import verify

    a = _labelled(mislabel=False, mystery=True)
    res = verify.verify_annotation(a)
    assert res["per_type"]["MysteryType"]["status"] == "unscoreable"
    assert "MysteryType" not in [f["cell_type"] for f in res["failed"]]


def test_suggest_reruns_ordered_advisory_no_mutation():
    from spatialscribe.analysis import verify

    a = _labelled(mislabel=True)
    before = a.obs["cell_type"].astype(str).tolist()
    res = verify.verify_annotation(a)
    sug = verify.suggest_reruns(res)                   # apply=False by default
    assert sug, "expected at least one suggestion for a failed type"
    from spatialscribe.analysis import capabilities as cap
    for it in sug:
        assert it["capability"] in cap.REGISTRY or it["capability"] is None
    assert any(it["action"] in {"abstain", "subcluster", "merge", "relabel", "review"} for it in sug)
    assert a.obs["cell_type"].astype(str).tolist() == before   # dry-run mutates nothing


def test_panel_gap_cause_yields_merge():
    from spatialscribe.analysis import verify

    a = _labelled(mislabel=True)
    # A valid panel_check that marks T cell <-> Myeloid unresolvable on the panel.
    a.uns["panel_check"] = {
        "coverage": {"T cell": {"status": "amber", "n_present": 1, "n_markers": 7},
                     "Myeloid": {"status": "green", "n_present": 3, "n_markers": 7},
                     "Endothelial": {"status": "green", "n_present": 3, "n_markers": 4}},
        "confusable_pairs": [{"pair": ["T cell", "Myeloid"]}], "merge_groups": [["T cell", "Myeloid"]]}
    res = verify.verify_annotation(a)
    assert res["per_type"]["T cell"]["cause"] == "panel_gap"
    sug = verify.suggest_reruns(res)
    merge = [it for it in sug if it["cell_type"] == "T cell" and it["action"] == "merge"]
    assert merge and merge[0]["capability"] is None and merge[0]["advisory"] is True


def test_low_coverage_without_confusable_pair_abstains_not_merges():
    # A type with sparse on-panel markers (amber) but NO genuine confusable pair must route to
    # abstain (low_signal), NOT a merge into the marker-argmax residual. This is the Mast case: the
    # old code merged Mast into Epithelial/Tumor with a false "no private marker" reason.
    from spatialscribe.analysis import verify

    a = _labelled(mislabel=True)
    # panel_check with amber coverage for the failing type but an EMPTY confusable_pairs list.
    a.uns["panel_check"] = {
        "coverage": {"T cell": {"status": "amber", "n_present": 1, "n_markers": 7},
                     "Myeloid": {"status": "green", "n_present": 3, "n_markers": 7},
                     "Endothelial": {"status": "green", "n_present": 3, "n_markers": 4}},
        "confusable_pairs": [], "merge_groups": []}
    res = verify.verify_annotation(a)
    # Without a genuine confusable pair, amber coverage must NOT be classified panel_gap (the old code
    # merged on amber coverage alone), and no merge action may be emitted for the failing type.
    assert res["per_type"]["T cell"]["cause"] != "panel_gap"
    sug = verify.suggest_reruns(res)
    tcell = [it for it in sug if it["cell_type"] == "T cell"]
    assert tcell
    assert not any(it["action"] == "merge" for it in tcell)   # no merge without a real confusable pair


def test_cross_compartment_confusable_pair_does_not_merge():
    # Even a genuine confusable pair must NOT merge across major compartments (immune into epithelial):
    # Mast <-> Epithelial/Tumor is unresolvable on a panel, but dissolving an immune granulocyte into
    # tumour epithelium is unsound - it must abstain instead.
    from spatialscribe.analysis import verify

    a = _labelled(mislabel=True)
    a.uns["panel_check"] = {
        "coverage": {"T cell": {"status": "amber", "n_present": 1, "n_markers": 7},
                     "Myeloid": {"status": "green", "n_present": 3, "n_markers": 7},
                     "Endothelial": {"status": "green", "n_present": 3, "n_markers": 4}},
        # T cell (immune) declared confusable ONLY with Endothelial (a different compartment).
        "confusable_pairs": [{"pair": ["T cell", "Endothelial"]}],
        "merge_groups": [["T cell", "Endothelial"]]}
    res = verify.verify_annotation(a)
    assert res["per_type"]["T cell"]["cause"] == "low_signal"           # cross-compartment -> not panel_gap
    sug = verify.suggest_reruns(res)
    assert not any(it["cell_type"] == "T cell" and it["action"] == "merge" for it in sug)


def test_self_verify_capability_produces_and_jsonable(processed_adata, ctx):
    from spatialscribe.analysis import capabilities as cap

    res = cap.run(processed_adata, "self_verify", {}, ctx)
    assert res.ok, res.error
    assert "annotation_verification" in processed_adata.uns
    json.dumps(res.value, default=str)
    assert "suggestions" in res.value


def test_self_verify_prereq_enforced(raw_adata, ctx):
    from spatialscribe.analysis import capabilities as cap

    res = cap.run(raw_adata, "self_verify", {}, ctx)      # raw has no cell_type
    assert not res.ok
    assert res.error["error_type"] == "prerequisite_missing"


def test_verify_key_and_capability_registered():
    from spatialscribe.analysis import capabilities as cap
    from spatialscribe.analysis import keys

    assert keys.Uns.ANNOTATION_VERIFICATION in keys.all_keys()
    assert "self_verify" in cap.copilot_names()


def test_autorerun_abstains_a_mislabel_so_it_leaves_the_failing_set():
    # The enhanced re-run loop: a failed marker-check (mislabel) must ACTUALLY act, not re-flag the
    # type forever. It abstains the disagreeing cells -> the type leaves the failing set, and the
    # correction is an honest abstention, never a new confident (wrong) lineage.
    from spatialscribe.analysis import annotate as an
    from spatialscribe.analysis import capabilities as cap
    from spatialscribe.analysis import verify

    a = _labelled(mislabel=True)
    ctx = cap.RunContext(tissue="melanoma", use_llm=False)
    before = a.obs["cell_type"].astype(str).to_numpy().copy()

    res = verify.autorerun(a, ctx, max_rounds=2)
    assert "T cell" in res["initial_failed"]
    assert "T cell" not in res["final_failed"], res          # it LEFT the failing set
    assert res["n_fixed"] >= 1

    after = a.obs["cell_type"].astype(str).to_numpy()
    changed = after != before
    assert changed.any(), "expected the mislabelled cells to be re-labelled"
    assert all(an.is_abstention(x) for x in after[changed])  # honest: abstained, never relabelled
    assert int((after == "Endothelial").sum()) == 60         # the clean lineage is untouched


def test_autorerun_clean_labels_take_no_action():
    from spatialscribe.analysis import capabilities as cap
    from spatialscribe.analysis import verify

    a = _labelled(mislabel=False)
    ctx = cap.RunContext(tissue="melanoma", use_llm=False)
    res = verify.autorerun(a, ctx, max_rounds=2)
    assert res["final_failed"] == []
    assert res["n_auto_actions"] == 0
