def test_build_query_pins_census_and_labels():
    import importlib.util, pathlib
    spec = importlib.util.spec_from_file_location(
        "fetch_reference",
        pathlib.Path(__file__).resolve().parents[1] / "scripts" / "fetch_reference.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    q = mod.build_query(tissue="skin of body", census_version="2025-01-30")
    assert q["species"] == "homo_sapiens"
    assert q["census_version"] == "2025-01-30"          # pinned = reproducible
    assert "cell_type_ontology_term_id" in q["column_names"]  # keep the CL id, not just free text
