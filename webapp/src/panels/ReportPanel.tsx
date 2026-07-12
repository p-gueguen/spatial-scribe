// ReportPanel: re-runnable analysis log rendered as copy-pasteable Python.
import React from "react";
import * as api from "../api";

export default function ReportPanel({ summary, sid }: { summary: any; sid?: string }) {
  const d = summary?.report;
  if (!d) return null;

  const log = Array.isArray(d.log) ? d.log : [];
  const lines = log.map((entry) =>
    entry && entry.code
      ? entry.code
      : 'cap.run(adata, "' + ((entry && entry.name) || "?") + '", ' +
        JSON.stringify((entry && entry.params) || {}) + ")"
  );

  return (
    <div className="panel">
      <h4>The whole run, packaged</h4>
      {sid && (
        <div className="pills" style={{ marginBottom: 10 }}>
          <a className="btn primary" style={{ textDecoration: "none", minHeight: "auto", padding: ".45rem .8rem" }}
             href={api.exportReportUrl(sid)}>Download HTML report</a>
          <a className="btn" style={{ textDecoration: "none", minHeight: "auto", padding: ".45rem .8rem" }}
             href={api.exportH5adUrl(sid)}>Annotated .h5ad</a>
          <a className="btn" style={{ textDecoration: "none", minHeight: "auto", padding: ".45rem .8rem" }}
             href={api.exportScriptUrl(sid)}>analysis.py</a>
        </div>
      )}
      <div className="pmuted">Every step you ran, as re-runnable code.</div>
      {lines.length === 0 ? (
        <div className="pmuted">Run some steps - they'll be captured here.</div>
      ) : (
        <div className="codeblock">{lines.join("\n")}</div>
      )}
    </div>
  );
}
