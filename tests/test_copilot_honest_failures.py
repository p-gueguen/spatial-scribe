"""issue #1: the copilot's answer must report a failed / empty tool call whatever the model says."""
import anndata
import numpy as np
import pandas as pd

from spatialscribe.agent import tools as agent_tools
from spatialscribe.analysis import llm


def _adata():
    a = anndata.AnnData(X=np.ones((10, 2), dtype="float32"), var=pd.DataFrame(index=["CD3E", "EPCAM"]))
    a.obs["cell_type"] = pd.Categorical(["T cell"] * 5 + ["Epithelial/Tumor"] * 5)
    a.obs_names = [f"c{i}" for i in range(10)]
    return a


def test_copilot_appends_failure_note_over_model_success_claim(monkeypatch):
    calls = iter([
        {"text": "", "assistant": {"role": "assistant", "content": "x"},
         "tool_calls": [{"id": "t1", "name": "marker_dotplot", "input": {"genes": ["NOTAGENE"]}},
                        {"id": "t2", "name": "no_such_capability", "input": {}}]},
        {"text": "The dot-plot has been drawn and is on the panel.",     # the confabulation
         "assistant": {"role": "assistant", "content": "done"}, "tool_calls": []},
    ])
    monkeypatch.setattr(llm, "tool_chat", lambda *a, **k: next(calls))
    reply = agent_tools.run_copilot(_adata(), "draw the marker dotplot", use_llm=True)
    assert "marker_dotplot drew nothing" in reply
    assert "no_such_capability FAILED" in reply


def test_openai_tool_result_error_is_stated_in_text():
    msgs = [{"role": "user", "content": [{"type": "tool_result", "tool_use_id": "t1",
                                          "content": '{"message": "Connection refused"}', "is_error": True}]}]
    out = llm._openai_msgs("SYS", msgs)
    tool = [m for m in out if m["role"] == "tool"][0]
    assert tool["content"].startswith("ERROR") and "Connection refused" in tool["content"]
