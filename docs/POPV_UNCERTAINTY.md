# popV consensus uncertainty in SpatialScribe

**What / why.** Ports the load-bearing idea of **popV** (Kimmel, Ergen, Yosef et al., *Nature Genetics*
2024, [10.1038/s41588-024-01993-3](https://doi.org/10.1038/s41588-024-01993-3)): across an ensemble of
annotation methods, the number that **agree** on a cell's label - the *consensus score* - tracks annotation
accuracy far better than any single method's self-reported certainty. popV's own finding: a consensus of
>=6/8 methods is >90% accurate; <=3/8 is <50%. It explicitly refuses to weight by individual methods'
certainties ("calibrated differently... futile"), and a **low consensus is the interpretable flag** for the
three cases that need review - an ambiguous continuum state, an **out-of-reference / novel** type, or a
wrong reference label.

**Why it fits SpatialScribe.** SpatialScribe already runs a natural popV ensemble (the marker+Claude cluster
label, plus per-cell RCTD / SingleR / scANVI / panhumanpy / TACCO), and already computed
`consensus_agreement` - but `apply_confidence` **ignored it** and instead down-weighted by single-method
self-scores (`scanvi_confidence`, `rctd_weight`), the exact thing popV proves is unreliable. So this port
both adds the well-calibrated signal and **corrects a backwards weighting**.

## What was built

| Piece | Where |
|-------|-------|
| Consensus scoring (`consensus_score`, `consensus_n_methods`, `consensus_agreement`, `consensus_reliability` bins) | `analysis/consensus.py` (`consensus_metrics`, `reliability_bin`) |
| Consensus folded into confidence as the **trusted primary factor**, but only with a diverse ensemble (`>= MIN_TRUST_METHODS` voted) | `analysis/annotate.py` (`apply_confidence`) |
| Plain-language **"annotation methods disagree"** abstention reason (low consensus + diverse ensemble) | `analysis/rejection.py` (`method_disagreement`) |
| Copilot readout: reliability-bin distribution + `pct_trusted_ensemble` | `analysis/capabilities.py` (`annotation_methods` tool) |
| Tunable thresholds (bins, min-methods, disagreement cutoff) | `docs/research/annotation_qc_thresholds.yaml` (`consensus_popv`) |

## Deliberate simplifications (be honest about these)

1. **Exact-label majority voting, not ontology-aware consensus.** popV weights agreement by cell-ontology
   distance. We match labels exactly. This is the same choice **LatchBio** made productionizing popV on PBMC
   data (simple majority voting scored ~92% vs the ontology consensus's ~93% - on par), and it avoids
   needing a cell-ontology graph.
2. **Calibration needs a diverse ensemble.** popV's guarantee depends on multiple *different* predictors.
   The reference-free demo runs only marker+Claude (2 methods), so the consensus is coarse and is surfaced
   but **not weighted** (the `MIN_TRUST_METHODS = 3` gate). It becomes load-bearing once the reference path
   (RCTD/SingleR/scANVI/TACCO) is on - exactly the regime popV was designed for.
3. **We port the idea, not the package.** Running popV itself needs a reference + GPU (scANVI/scArches
   training). We reuse the labels the existing subprocess annotators already emit, keeping it reference-free
   and CPU. Running popV-proper as an optional subprocess (like RCTD/scANVI) is a clean future extension.

## Reliability-weighted consensus: when naive majority is the wrong call (measured 2026-07-10)

popV is right that weighting by each method's **own self-reported certainty** is futile: those scores are
calibrated differently and are not comparable across methods. **Reliability** is a different quantity - a
method's accuracy measured against labels, on one shared scale - and it does matter.

Live 3-method run on a 3k-cell CosMx breast subsample, scored against the section's independent 49-type
ground truth on a shared lineage axis:

| annotator                | lineage accuracy |
| ------------------------ | ---------------: |
| RCTD-doublet             |        **0.923** |
| SingleR                  |            0.667 |
| panhumanpy               |            0.272 |
| *naive majority consensus* |      *0.858* |

The consensus is **worse than the best single method**, because two weak voters dilute the one that is
right. `consensus.weighted_vote` fixes this, opt-in, and the weight must be the **log-odds** of each
method's reliability, not its raw accuracy:

* linear weights still get it wrong: RCTD 0.923 < SingleR 0.667 + panhumanpy 0.272 = 0.939 → picks the
  wrong label. Measured on the real labels, linear-weighted scores **0.8522**, identical to naive majority.
* log-odds weights get it right: `logit(0.92) = 2.44 > logit(0.67) = 0.71 + max(0, logit(0.27)) = 0`.
  This is the optimal weighting for independent noisy voters (Nitzan–Paroush). A method no better than a
  coin flip clamps to weight 0 and is ignored rather than allowed to vote against the field.

Measured with weights fit on a train half and scored on the held-out half: **reliability-weighted 0.9200 =
RCTD alone 0.9200**, versus naive majority 0.8522.

Read that honestly: when one method dominates this hard, reliability weighting **collapses toward it**. It
does not manufacture a super-consensus that beats the best method; it stops the weak methods dragging the
strong one down. **Consensus only adds value when the methods are comparable.**

Use `consensus.reliability_from_labels(adata, cols, gt_col)` whenever a labeled calibration set exists; it
is the only path allowed to call the weighting calibrated. Otherwise the documented
`consensus.DEFAULT_RELIABILITY` prior applies, in which `scanvi_label` (0.35) and the cluster label
`cell_type` (0.60) are **unmeasured guesses**, flagged as such in the source. The default remains naive
majority, so nothing changes unless you ask for it.

## Interpretation surfaced to the user

`consensus_reliability` bins (agreement fraction, **recalibrated on real data - see below**): **very_high**
>=0.9, **high** >=0.7, **moderate** >=0.5, **low** below - with `single` when only one method voted. The
copilot narrates, e.g., "5 of 6 methods agree (high) - trust this call" vs "2 of 6 agree (low) - methods
disagree, which popV associates with an ambiguous state, a novel type, or a bad reference label; worth review."

## Real-data validation (2026-07-08)

Validated on a REAL labeled dataset, no simulation: a CELLxGENE hematopoietic-progenitor scRNA-seq set
(27,998 cells, 12 Cell-Ontology types along a differentiation continuum), split reference/calibration/test,
5-classifier ensemble (LogReg/RF/KNN/HistGB/MLP). Findings:

- **Consensus predicts accuracy (holds):** monotonic - 1 method agreeing -> 25% accurate, all 5 -> 90%.
- **Honest limit:** consensus is NOT a better error-*ranker* than single-method confidence here (AUPRC 0.89
  vs the best classifier's 0.95). With only 5 well-calibrated classifiers the coarse 1-5 consensus can't
  out-rank a continuous probability; popV's advantage needs a larger, more heterogeneous ensemble.
- **Calibration:** 0.6 agreement ("high" under the old guess) was only ~44% accurate, so the bins were
  tightened 0.875/0.6/0.4 -> 0.9/0.7/0.5 (provisional; one hard tissue - re-tune per dataset).

Run paired with scConform (whose conformal coverage held almost exactly on the same split). Scripts + report:
`/data/spatial-scribe/popv_scconform_validation/`.

## Sources
popV [Nat Genet 2024]; LatchBio "Benchmarking popV Ensemble Cell Type Annotations on CS Genetics PBMC Data"
(blog.latch.bio) - see `docs/research/sources/popV_annotation.md`.
