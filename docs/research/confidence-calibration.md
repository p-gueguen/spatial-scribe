# Annotation confidence: what it means, and what it does not

`annotate.apply_confidence` writes a per-cell `annotation_confidence`. This document records what that
number was measured to be worth, on three real sections with independent expert ground truth, and what you
are therefore allowed to claim about it.

The short version: **the confidence score is an ordinal heuristic with almost no per-cell discrimination.
Post-hoc calibration makes the number honest, but it mostly just learns the section's average accuracy.**

## Two defects, only one of them fixable

A confidence score can fail in two independent ways, and the benchmark conflated them at first.

* **Calibration** asks *does the number mean `P(correct)`?* This is fixable after the fact: fit a monotone
  isotonic map from `(confidence, correct)` on labeled cells. `analysis/calibration.py` does this.
* **Discrimination** asks *does the number RANK correct cells above incorrect ones?* This is **not** fixable
  post hoc, because a monotone map cannot reorder cells. Abstention depends entirely on discrimination.

`calibration.auc` measures the second. Note isotonic can still *raise* AUC slightly, because it is
non-strictly increasing and its flat regions merge locally anti-correlated cells into ties. It never
reorders them.

## Measured, on three sections with independent ground truth

All labels coarsened onto one 6-lineage axis (T/NK, B/Plasma, Myeloid, Endothelial, Stromal, Epithelial);
`Other` excluded. Calibrator fit on one random half, every number below scored on the held-out half.

| section              | genes  | base acc | ECE raw → calibrated | pooled AUC | **within-lineage AUC** | **Brier skill** |
| -------------------- | -----: | -------: | -------------------- | ---------: | ---------------------: | --------------: |
| CosMx breast         | ~1,000 |    0.664 | 0.3315 → 0.0116      |     0.4985 |                 0.4910 |         +0.0200 |
| Atera (WTA) breast   | 18,028 |    0.899 | 0.4573 → 0.0042      | **0.4191** |                 0.4790 |         +0.0011 |
| Atera (WTA) cervical | 18,028 |    0.742 | 0.2751 → 0.0041      | **0.5780** |                 0.5712 |         +0.0167 |

These are internal validation sections; only the aggregate metrics below are shown, and no raw
data is redistributed in this repo. Ground truth: CosMx ships a 49-type annotation; both Atera demos ship 10x expert `cell_groups.csv`
(real cell types, not cluster ids), 99.8% / 99.0% of which map cleanly onto the lineage axis.

### The ECE headline is largely vacuous

A **constant** predictor that says "every cell is correct with probability = the section's accuracy" is
already near-perfectly calibrated: it scores ECE 0.0002–0.0034 on these sections, i.e. **better than
isotonic**. So a low `ece_after` proves almost nothing on its own.

Brier punishes a constant, so `brier_skill = 1 - brier_after / brier_baserate` is the honest headline. It is
**+0.020 / +0.001 / +0.017**. On both breast sections the fitted calibrator is literally degenerate: 100% of
held-out predictions land within 0.02 of the base rate (prediction sd 0.004–0.010).

`calibration.report` therefore always emits `ece_baserate`, `brier_baserate` and `brier_skill`, and the
`calibrate_confidence` capability says in plain words when skill < 0.01 that calibration has reduced to
reporting the section's accuracy for every cell.

### Pooled AUC is Simpson-confounded; use `auc_within`

Atera breast inverts (pooled AUC 0.4191, bootstrap 95% CI [0.4115, 0.4272]) for a specific reason: the
pipeline is **most confident on the lineage it labels worst and least confident on the one it nails**.

| predicted lineage | n      | accuracy | mean confidence |
| ----------------- | -----: | -------: | --------------: |
| Epithelial        | 29,212 |    0.986 |           0.439 |
| Endothelial       |  4,782 |    0.527 |           0.503 |
| T/NK              |  4,041 |    0.583 |           0.488 |

`corr(lineage mean-conf, lineage accuracy) = -0.79`. Decomposing: the across-lineage term carries AUC 0.134
while the within-lineage signal is a flat 0.479. So the pooled inversion is an **across-group artifact**,
not a per-cell property.

A per-cell gate can only exploit the *within-group* term, so that is the number to report.
`calibration.auc_within` rank-normalises confidence inside each predicted cell type first. Across the three
sections it is **0.491 / 0.479 / 0.571** — barely above chance anywhere.

> **Do not** use `corr(lineage-conf, lineage-acc)` as a diagnostic for whether the gate is usable. It looks
> compelling on two sections and is falsified by the third: CosMx has corr **+0.94** yet pooled AUC 0.4985.

### A calibrator does not transfer between sections

Fit on section A, apply to section B, for all six ordered pairs: transfer ECE regresses on
`|Δ base-rate|` with **slope 1.07, R² 0.982, intercept ≈ 0**, and matches "just ship A's base rate as a
constant" to within 0.008 on average. A transferred calibrator exports only the source section's accuracy —
the one quantity that requires B's labels.

**Consequence:** `calibrate_confidence` is deliberately label-bound and opt-in. There is no reference-transfer
path, and adding one would be dishonest without a covariate-shift caveat that swallows the result.

## The abstention gate

* The **FAIL cut** (`confidence < 0.25`) is effectively inert: it fired on 0, 2 and 7 cells across the three
  sections. The low-quality gate (`counts < floor` or `genes < RICH_PANEL_MIN_GENES` - 15 on these >=1000-gene
  panels, 5 on targeted) that shares the FAIL verdict is separate and legitimate.
* The **WARN cut** (`confidence < 0.5`) greys 63–78% of cells, and **on two of three sections the greyed
  cells are more accurate than the confident ones** (CosMx 0.721 vs 0.460; Atera breast 0.918 vs 0.873;
  Atera cervical is the only well-behaved one, 0.762 vs 0.799).

The sign is section-dependent, reproduces across an independent subsample seed *and* clustering resolution
(breast 0.4216 → 0.4264; cervical 0.5783 → 0.5738), and survives matching the base accuracy across sections
(breast 0.4197 at a matched 0.70). It is therefore **not** a clustering or class-imbalance artifact, and on an
unlabeled query section you cannot know which regime you are in.

**Open decision:** the WARN greying is display-only today. Given `auc_within` ≈ 0.5 everywhere, there is
almost no per-cell signal for it to display.

## The AQI: one section-level quality index (validated)

The per-cell confidence above barely ranks cells (`auc_within` ≈ 0.5). The **AQI**
(`spatial_anno_metrics.annotation_quality_index`, v0.4.1) answers the *other* question: is this whole
section's annotation any good, on a single 0–1 scale that means the same thing across tissues, panels and
depths? It is an **index** (section-comparable, monotone in quality), **not** `P(correct)` — the
no-transfer finding above forbids a universal reference-free accuracy curve.

**Formula:** `AQI = w_coh · min(A, soft-min_{p=-4}(C, M))`.
* **C** — contamination/purity (macro-median PMP × retention; panel-size invariant).
* **M** — marker-program fidelity (macro one-vs-rest AUC over `n ≥ 50` types, rescaled `2(auc-0.5)`).
* **A** — panel × depth adequacy (base-rate-normalised, depth-matched macro F1 from `panel_resolvability`
  on a reference). A **CEILING**, not a driver: you cannot type better than the panel resolves, but a rich
  panel does not by itself make the labels right.
* **H** — spatial/internal coherence, a bounded ≤15% multiplier `w_coh` (coherence ≠ correctness).

**Validated** (`quality_jobs_2026-07-10/validate_aqi.py`, three sections with independent expert GT, all
coarsened to a common 6-lineage axis). Section-level Spearman of each component vs true balanced accuracy:

| component | ρ(component, true accuracy) | role in the index |
|---|---:|---|
| M marker fidelity | **+1.0** | drives (soft-min) |
| C contamination | +0.5 | drives (soft-min) |
| A panel/depth adequacy | **−1.0** | **ceiling only** (`min`) |
| G cross-method agreement | +0.5 | **not a term** (→ abstention) |

A **anti-orders** accuracy because it is a resolvable ceiling — high on any deep panel regardless of the
labels actually assigned (CosMx has the highest A = 0.98 but the lowest accuracy, 0.33). Folding it into the
soft-min bottleneck, as the first design did (`soft-min(A,C,G,M)`), gave ρ = +0.5 and the **wrong** section
order. Making A a ceiling and leaning on C+M gives **ρ(AQI, true accuracy) = 1.0**, isotonic residual 0.000
vs 0.079 for the mean-confidence baseline. The index is then **voter-independent** (it does not change when
you add/remove reference methods).

**Cross-method agreement G is a within-section signal, not a cross-section term.** Its absolute level tracks
*reference quality*, not accuracy, so it does not transfer across sections — the same reason a calibrator
does not transfer (above). It is reported under `abstention {signal, n_voters, available}`, where it does
its real job: within a section it ranks correct-vs-wrong cells at **mean per-cell AUC 0.77**, and — crucially
— strongest exactly where accuracy is lowest (CosMx 0.33 accuracy → abstention AUC 0.91), which is where you
actually want to grey cells. It needs ≥3 diverse voters (`consensus.MIN_TRUST_METHODS`); with fewer, the
`regime` is `index_only` and the greying is advisory (matches the WARN-gate finding above).

**In the app:** `qc.run_funnel` populates `annotation_quality.aqi` on every labelled section — C+M always,
A only when a reference is already loaded (never triggers a multi-GB atlas load just for the ceiling), and
the abstention signal when ≥3 of RCTD/SingleR/scANVI/panhumanpy have voted. The Annotate → confidence tab
shows the AQI headline, what limits it (`argmin`), and whether the abstention signal is trustworthy.

## Using it

```python
from spatialscribe.analysis import capabilities as cap

# Needs a ground-truth obs column on the SAME label axis as pred_key. Without one it skips, by design.
res = cap.run(adata, "calibrate_confidence",
              {"truth_key": "_gt_lineage", "pred_key": "cell_type", "cal_frac": 0.5}, ctx)

res.value["ece_before"], res.value["ece_after"]     # calibration
res.value["brier_skill"]                            # did we learn anything beyond the base rate?
res.value["auc_within"]                             # can a per-cell gate work at all?
res.value["gate"]                                   # threshold picked on the fit half, scored held-out
```

Writes `obs['annotation_confidence_calibrated']` and `uns['calibration_report']`. It **never** overwrites the
raw `annotation_confidence`. `apply_confidence(use_calibrated=True)` opts the PASS/WARN/FAIL gate onto the
calibrated score; the default is unchanged.

## What you may claim

* ✅ "Confidence is an ordinal heuristic." Always true.
* ✅ "Calibrated against labels, the score means `P(correct)`" — only after `calibrate_confidence` ran with a
  `truth_key`, and only for that section.
* ❌ "Calibration improved the annotation." It does not change a single label.
* ❌ "Low ECE means the confidence is informative." Quote `brier_skill` and `auc_within` instead.
* ❌ "High-confidence cells are more likely correct." Section-dependent, and false on 2 of 3 sections tested.
* ✅ "The AQI is a section-level quality **index** that orders sections by true accuracy." ρ = 1.0 on 3 GT
  sections. Say "index", never `P(correct)` or an accuracy rate.
* ✅ "Cross-method agreement ranks correct cells **within** a section (abstention), not across sections."
* ❌ "A higher AQI means a higher fraction of cells are correct." It is monotone-in-quality, not a rate; the
  no-transfer finding forbids a universal reference-free accuracy curve.

## Gotchas

1. **Never report an ECE improvement without the base-rate null.** A constant beats isotonic on ECE here.
   Report `brier_skill` alongside, or the result is unfalsifiable.
2. **Never report a pooled AUC for a grouped prediction.** Report `auc_within` too; the pooled figure inverted
   on a real section purely through Simpson's paradox.
3. **Pick abstention thresholds on the fit half, score them on the held-out half.** Doing both on the same
   cells reports a cherry-picked gain. `calibrate_confidence` enforces this and labels the split.
4. **Deprecated codewords are not negative controls.** On Xenium Prime 5K they carry ~15% of every cell's
   counts; counting them as controls inflates `pct_counts_control`, zeroes the confidence and abstains ~95%
   of a section. Atera (18k) and CosMx ship none. Use `io.build_neg_control_mask` (strict), not
   `build_control_mask` (broad, for panel-gene selection).
5. **Never put panel adequacy (A) or cross-method agreement (G) in the AQI soft-min bottleneck.** A is a
   resolvable *ceiling* that ANTI-orders accuracy (ρ = −1.0 on the 3 GT sections — high on any deep panel
   regardless of the labels), so it may only `min`-cap. G's absolute level tracks reference quality and does
   not transfer across sections, so it is a *within-section* abstention signal, not a cross-section term. The
   index is `min(A, soft-min(C, M))`; C and M are the terms that track accuracy. Mixing A/G into the
   bottleneck (`soft-min(A,C,G,M)`) drops ρ(AQI, accuracy) from 1.0 to 0.5 and mis-orders the sections.

## Reproducing

`calibration.py` is pure and dependency-light (numpy, sklearn's `IsotonicRegression`, scipy `rankdata`).
Tests in `tests/test_calibration.py` pin every claim above, including a synthetic Simpson fixture where
pooled AUC < 0.25 while `auc_within` correctly reports ≈ 0.5.

See also: [cell-annotation-qc.md](cell-annotation-qc.md) (the 6-layer QC funnel),
[POPV_UNCERTAINTY.md](../POPV_UNCERTAINTY.md) (consensus agreement as the complementary signal).
