"""The malignant-caller gate: an explicit answer beats the tissue-keyword guess.

Nothing can know whether a section contains tumour before the malignant callers run, so the Data
step asks. Before this gate existed, the ONLY signal was substring-matching the free-text tissue
against ``cnv._TUMOUR_KEYWORDS``, which is wrong in both directions - these tests pin exactly that.
"""
from spatialscribe.analysis.cnv import is_tumour_context


def test_explicit_answer_overrides_the_keyword_guess():
    # Ticked on a tissue the keywords would never match, and unticked on one they always match.
    assert is_tumour_context("glioblastoma", is_tumour=True) == (True, "explicit")
    assert is_tumour_context("normal breast", is_tumour=False) == (False, "explicit")


def test_unanswered_falls_back_to_the_keyword_guess():
    run, gate = is_tumour_context("uveal melanoma", is_tumour=None)
    assert (run, gate) == (True, "tissue_keyword")
    run, gate = is_tumour_context("mouse brain", is_tumour=None)
    assert (run, gate) == (False, "tissue_keyword")


def test_the_keyword_guess_is_wrong_in_both_directions():
    """Documents WHY the checkbox exists. If someone "fixes" the keyword list, this test tells them
    what it used to get wrong; it is not asserting the mistakes are desirable."""
    # False positives: normal tissue trips malignant calling.
    for normal in ("normal breast", "healthy skin biopsy", "lung (healthy donor)"):
        assert is_tumour_context(normal, is_tumour=None)[0] is True
        assert is_tumour_context(normal, is_tumour=False)[0] is False   # the checkbox rescues it

    # False negatives: real tumours are silently skipped.
    for tumour in ("glioblastoma", "sarcoma", "lymphoma", "leukemia", "neuroblastoma", "myeloma"):
        assert is_tumour_context(tumour, is_tumour=None)[0] is False
        assert is_tumour_context(tumour, is_tumour=True)[0] is True     # the checkbox rescues it


def test_concordance_reports_which_gate_decided(processed_adata):
    """A skipped result must say WHY, so a report can never read a gate that did not fire as
    'this section has no malignant cells'."""
    from spatialscribe.analysis import cnv

    out = cnv.call_malignant_concordance(processed_adata, tissue="melanoma", is_tumour=False)
    assert out["status"].startswith("skipped")
    assert out["tumour_context"] is False
    assert out["gate"] == "explicit"
    assert "unticked" in " ".join(out["notes"])

    out = cnv.call_malignant_concordance(processed_adata, tissue="mouse brain", is_tumour=None)
    assert out["gate"] == "tissue_keyword"
