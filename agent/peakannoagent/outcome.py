"""Shared run completion policy for manifests, reports, CLI and the web app."""

from .config import RunConfig
from .evidence import EvidenceStore


def assess_run(store: EvidenceStore, config: RunConfig, agent_error=None) -> dict:
    required = ["peak_qc", "annotate_peaks", "enrich_pathways"]
    if not config.biology.is_broad_mark:
        required.append("find_motifs")
    if config.biology.rnaseq:
        required.append("integrate_rnaseq")
    latest = {ev.tool: ev for ev in store.items.values()}
    missing = [name for name in required if name not in latest]
    unsuccessful = [name for name in required if name in latest and latest[name].status != "ok"]
    incomplete_methods = [ev.tool for ev in latest.values() if ev.summary.get("failed_methods")]
    errors = [ev.tool for ev in latest.values() if ev.status == "error"]
    fundamental_failure = any(name in missing + unsuccessful for name in ("peak_qc", "annotate_peaks"))
    status = "failed" if fundamental_failure else (
        "partial" if missing or unsuccessful or errors or incomplete_methods or agent_error else "success")
    reasons = ([f"Required step not run: {name}" for name in missing] +
               [f"Required step {name}: {latest[name].error or latest[name].status}" for name in unsuccessful] +
               [f"Step {name}: {latest[name].error or 'error'}" for name in errors if name not in unsuccessful])
    if agent_error:
        reasons.append(f"Agent did not complete normally: {agent_error}")
    for name in incomplete_methods:
        reasons.append(f"Some {name} methods failed: {latest[name].summary['failed_methods']}")
    return {"status": status, "required_steps": required, "missing_steps": missing,
            "unsuccessful_steps": unsuccessful, "reasons": reasons}
