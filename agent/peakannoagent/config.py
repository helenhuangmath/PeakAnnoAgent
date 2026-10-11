"""Run configuration: the peak file, genome, and the user's known biology."""

from __future__ import annotations

import json
import os
import re
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import List, Optional

SUPPORTED_GENOMES = ("mm10", "mm39", "hg38", "hg19")
GENOME_ALIASES = {"mouse": "mm10", "human": "hg38", "grcm38": "mm10", "grch38": "hg38",
                  "grcm39": "mm39", "grch37": "hg19"}
ORGANISM = {"mm10": "Mus musculus", "mm39": "Mus musculus", "hg38": "Homo sapiens", "hg19": "Homo sapiens"}
DEFAULT_MODEL = "claude-opus-5-5"

# Histone marks that form broad domains: summit-centred motif scans are not meaningful.
BROAD_MARKS = ("h3k27me3", "h3k9me3", "h3k36me3", "h3k79me2", "h4k20me1")


def normalize_genome(genome: str) -> str:
    g = GENOME_ALIASES.get(genome.lower(), genome.lower())
    if g not in SUPPORTED_GENOMES:
        raise ValueError(f"Unsupported genome '{genome}'. Choose one of {', '.join(SUPPORTED_GENOMES)}.")
    return g


def split_list(value) -> List[str]:
    if not value:
        return []
    if isinstance(value, (list, tuple)):
        return [str(v).strip() for v in value if str(v).strip()]
    return [v.strip() for v in re.split(r"[,;\s]+", str(value)) if v.strip()]


@dataclass
class BiologyContext:
    """What the user already knows. Free text plus structured hints the agent can act on."""

    description: str = ""
    assay: str = ""                      # e.g. "TF ChIP-seq", "ATAC-seq", "H3K27ac CUT&Tag"
    target: str = ""                     # ChIPped TF or histone mark
    cell_type: str = ""
    condition: str = ""                  # treatment / comparison, e.g. "LCMV clone 13 d30 vs Arm"
    expected_genes: List[str] = field(default_factory=list)
    expected_tfs: List[str] = field(default_factory=list)
    custom_gmt: Optional[str] = None     # user signatures (GMT) to test alongside MSigDB
    rnaseq: Optional[str] = None         # differential-expression table to integrate

    @property
    def is_broad_mark(self) -> bool:
        text = f"{self.target} {self.assay}".lower()
        return any(m in text for m in BROAD_MARKS)

    @property
    def is_tf_assay(self) -> bool:
        text = f"{self.assay} {self.target}".lower()
        return bool(self.target) and not self.is_broad_mark and not re.search(r"h[2-4]k|h2a|atac|dnase", text)

    def is_empty(self) -> bool:
        return not any([self.description, self.assay, self.target, self.cell_type, self.condition,
                        self.expected_genes, self.expected_tfs])


@dataclass
class RunConfig:
    peaks: str
    genome: str
    out_dir: str
    prefix: str = "run"
    biology: BiologyContext = field(default_factory=BiologyContext)
    blacklist: Optional[str] = None
    offline: bool = False                # skip ChIP-Atlas / literature / GEO
    use_llm: bool = True
    model: str = DEFAULT_MODEL
    effort: str = "high"
    max_turns: int = 40
    rscript: str = field(default_factory=lambda: os.environ.get("PEAKANNOAGENT_RSCRIPT", "Rscript"))
    chipatlas_wait_min: float = 25.0
    started_at: str = field(default_factory=lambda: datetime.now().isoformat(timespec="seconds"))

    def __post_init__(self):
        self.genome = normalize_genome(self.genome)
        self.peaks = str(Path(self.peaks).expanduser().resolve())
        self.out_dir = str(Path(self.out_dir).expanduser().resolve())
        if not Path(self.peaks).is_file():
            raise FileNotFoundError(f"Peak file not found: {self.peaks}")
        for attr in ("custom_gmt", "rnaseq"):
            value = getattr(self.biology, attr)
            if value:
                path = Path(value).expanduser().resolve()
                if not path.is_file():
                    raise FileNotFoundError(f"{attr} file not found: {path}")
                setattr(self.biology, attr, str(path))

    @property
    def organism(self) -> str:
        return ORGANISM[self.genome]

    def to_dict(self) -> dict:
        return asdict(self)


def load_biology_file(path: str) -> dict:
    """A biology file is JSON (structured fields) or any text/markdown (free description)."""
    text = Path(path).expanduser().read_text()
    if path.lower().endswith(".json"):
        data = json.loads(text)
        if not isinstance(data, dict):
            raise ValueError("Biology JSON must be an object.")
        return data
    return {"description": text.strip()}


def build_biology(file_fields: dict, **overrides) -> BiologyContext:
    merged = dict(file_fields)
    for key, value in overrides.items():
        if value not in (None, "", []):
            merged[key] = value
    known = set(BiologyContext.__dataclass_fields__)
    unknown = set(merged) - known
    if unknown:
        raise ValueError(f"Unknown biology field(s): {', '.join(sorted(unknown))}")
    merged["expected_genes"] = split_list(merged.get("expected_genes"))
    merged["expected_tfs"] = split_list(merged.get("expected_tfs"))
    return BiologyContext(**merged)
