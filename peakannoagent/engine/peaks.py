"""Reading peak files. Coordinates are zero-based and half-open, as in BED, and stay that way."""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import List

import numpy as np


@dataclass
class Peaks:
    ids: List[str]
    chrom: np.ndarray          # object array of chromosome names
    start: np.ndarray          # int64, zero-based
    end: np.ndarray            # int64, exclusive

    def __len__(self) -> int:
        return len(self.ids)

    def subset(self, keep) -> "Peaks":
        keep = np.asarray(keep)
        idx = np.flatnonzero(keep) if keep.dtype == bool else keep
        return Peaks([self.ids[i] for i in idx], self.chrom[idx], self.start[idx], self.end[idx])

    def write_bed(self, path) -> str:
        with open(path, "w") as handle:
            for i in range(len(self)):
                handle.write(f"{self.chrom[i]}\t{self.start[i]}\t{self.end[i]}\t{self.ids[i]}\n")
        return str(path)


def _whole_number(text: str) -> int:
    try:
        return int(text)
    except ValueError:
        value = float(text)                      # raises ValueError for non-numbers
        if value != int(value):
            raise ValueError("Peak coordinates must be whole numbers.") from None
        return int(value)


def _is_number(text: str) -> bool:
    try:
        float(text)
        return True
    except ValueError:
        return False


def read_peak_file(path) -> Peaks:
    """BED, narrowPeak/broadPeak, or a tab-separated table with Chr/Start/End columns."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        lines = [ln.rstrip("\r\n") for ln in handle
                 if ln.strip() and not ln.startswith(("#", "track", "browser"))]
    if not lines:
        raise ValueError("The peak file contains no peaks.")
    rows = [ln.split("\t") if "\t" in ln else ln.split() for ln in lines]
    header = None
    first = rows[0]
    if len(first) >= 3 and not (_is_number(first[1]) and _is_number(first[2])):
        header, rows = [h.strip() for h in first], rows[1:]
    if any(len(r) < 3 for r in rows):
        raise ValueError("Peak file must contain at least three columns: Chr, Start, End.")
    c_chr, c_start, c_end, c_id = 0, 1, 2, (3 if rows and len(rows[0]) > 3 else None)
    if header:
        lower = [h.lower() for h in header]
        if all(k in lower for k in ("chr", "start", "end")):
            c_chr, c_start, c_end = lower.index("chr"), lower.index("start"), lower.index("end")
            c_id = next((lower.index(k) for k in ("peakid", "name", "peak") if k in lower), None)
    if not rows:
        raise ValueError("The peak file contains no peaks.")
    chrom = np.array([r[c_chr] for r in rows], dtype=object)
    start = np.array([_whole_number(r[c_start]) for r in rows], dtype=np.int64)
    end = np.array([_whole_number(r[c_end]) for r in rows], dtype=np.int64)
    ids = [r[c_id] if c_id is not None and c_id < len(r) else "" for r in rows]
    if not all(ids) or len(set(ids)) != len(ids) or all(i in (".", "") for i in ids):
        ids = [f"peak_{i + 1:05d}" for i in range(len(rows))]
    return Peaks(ids, chrom, start, end)


def check_intervals(peaks: Peaks) -> None:
    if len(peaks) and (np.any(peaks.start < 0) or np.any(peaks.end <= peaks.start)):
        raise ValueError("Peaks must use valid zero-based, half-open intervals: 0 <= Start < End.")
