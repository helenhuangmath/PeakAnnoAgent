"""Evidence store: every tool result gets an ID (E1, E2, ...) that the narrative
must cite. Literature (PMID:n) and public datasets (GSE/SRX accessions) are
registered too, so the report can check that every citation is real."""

from __future__ import annotations

import json
import re
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

CITATION_RE = re.compile(r"\[((?:E\d+|PMID:\s*\d+|G[SD][EM]\d+|SRX\d+|DRX\d+|ERX\d+)"
                         r"(?:\s*[,;]\s*(?:E\d+|PMID:\s*\d+|G[SD][EM]\d+|SRX\d+|DRX\d+|ERX\d+))*)\]")

# Benign Bioconductor warnings that would only distract in the report.
NOISY_R_WARNINGS = ("out-of-bound range", "sequence levels in common", "select()' returned",
                    "built under R version")


def _as_list(value) -> list:
    """jsonlite unboxes length-1 vectors, so a single warning arrives as a bare string."""
    if value is None:
        return []
    return [value] if isinstance(value, str) else list(value)


@dataclass
class Evidence:
    id: str
    tool: str
    title: str
    status: str                      # ok | error | skipped
    params: dict
    summary: dict = field(default_factory=dict)
    tables: Dict[str, str] = field(default_factory=dict)
    figures: Dict[str, str] = field(default_factory=dict)
    error: Optional[str] = None
    warnings: List[str] = field(default_factory=list)
    rationale: str = ""
    elapsed_sec: float = 0.0
    versions: dict = field(default_factory=dict)


class EvidenceStore:
    def __init__(self):
        self.items: Dict[str, Evidence] = {}
        self.references: Dict[str, dict] = {}   # "PMID:123" -> metadata
        self.datasets: Dict[str, dict] = {}     # "GSE123" / "SRX123" -> metadata
        self.trace: List[dict] = []             # every tool call, in order

    def add(self, tool: str, title: str, params: dict, outcome: dict, rationale: str = "") -> Evidence:
        result = outcome.get("result") or {}
        ev = Evidence(
            id=f"E{len(self.items) + 1}",
            tool=tool,
            title=title,
            status=outcome.get("status", "ok"),
            params={k: v for k, v in params.items() if k != "rationale"},
            summary={k: v for k, v in result.items()
                     if k not in ("tables", "figures", "pdfs", "versions", "warnings")},
            tables={k: v for k, v in (result.get("tables") or {}).items() if isinstance(v, str)},
            figures={k: v for k, v in (result.get("figures") or {}).items() if isinstance(v, str)},
            error=outcome.get("error"),
            warnings=_as_list(result.get("warnings")) +
            [w for w in _as_list(outcome.get("r_warnings")) if not any(n in w for n in NOISY_R_WARNINGS)],
            rationale=rationale,
            elapsed_sec=float(outcome.get("elapsed_sec") or 0.0),
            versions=result.get("versions") or {},
        )
        self.items[ev.id] = ev
        self.trace.append({"evidence": ev.id, "tool": tool, "status": ev.status, "rationale": rationale,
                           "time": time.strftime("%H:%M:%S")})
        return ev

    def add_reference(self, pmid: str, meta: dict) -> str:
        key = f"PMID:{pmid}"
        self.references.setdefault(key, meta)
        return key

    def add_dataset(self, accession: str, meta: dict) -> str:
        self.datasets.setdefault(accession, meta)
        return accession

    def latest(self, tool: str, ok_only: bool = True) -> Optional[Evidence]:
        for ev in reversed(list(self.items.values())):
            if ev.tool == tool and (ev.status == "ok" or not ok_only):
                return ev
        return None

    def known_citation(self, token: str) -> bool:
        token = token.replace(" ", "")
        if token.startswith("PMID:"):
            return token in self.references
        if re.fullmatch(r"E\d+", token):
            return token in self.items
        return token in self.datasets

    def check_citations(self, text: str) -> List[str]:
        """Return citation tokens in `text` that do not match any registered evidence."""
        unknown = []
        for match in CITATION_RE.finditer(text or ""):
            for token in re.split(r"\s*[,;]\s*", match.group(1)):
                if not self.known_citation(token):
                    unknown.append(token.replace(" ", ""))
        return sorted(set(unknown))

    def compact(self, ev: Evidence, limit: int = 12000) -> str:
        """JSON handed back to the model for one tool call, bounded in size."""
        payload = {"evidence_id": ev.id, "status": ev.status}
        if ev.error:
            payload["error"] = ev.error
        payload.update(ev.summary)
        if ev.warnings:
            payload["warnings"] = ev.warnings
        if ev.tables:
            payload["tables"] = {k: Path(v).name for k, v in ev.tables.items()}
        text = json.dumps(payload, default=str)
        # Shrink the longest lists until the payload fits, so the model always gets valid JSON.
        keep = 20
        while len(text) > limit and keep >= 1:
            for key, value in payload.items():
                if isinstance(value, list) and len(value) > keep:
                    payload[key] = value[:keep] + [f"... {len(value) - keep} more (use inspect_table)"]
            text = json.dumps(payload, default=str)
            keep //= 2
        return text

    def save(self, path: str) -> None:
        data = {
            "evidence": {k: asdict(v) for k, v in self.items.items()},
            "references": self.references,
            "datasets": self.datasets,
            "trace": self.trace,
        }
        Path(path).write_text(json.dumps(data, indent=2, default=str))
