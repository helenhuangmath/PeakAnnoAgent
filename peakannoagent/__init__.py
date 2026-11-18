"""PeakAnnoAgent: fast peak annotation, pathways, motifs and regulators, with an optional LLM agent.

One Python package. `peakannoagent.engine` holds the analysis (gene-model references, the fast
annotator, pathway enrichment, pathway networks, motif enrichment); the rest runs it as tools,
adds public-data and literature search, lets Claude plan and interpret, and writes the HTML report.

    from peakannoagent import build_reference, read_peak_file, annotate_peaks, run_annotation
    ref = build_reference("hg38", "gencode")
    table = annotate_peaks(read_peak_file("peaks.bed"), ref)
"""

__version__ = "1.1.0"

from .engine import (Engine, Peaks, Reference, annotate_peaks, build_reference, read_peak_file,  # noqa: E402
                     run_annotation, run_enrichment, run_peak_qc)

__all__ = ["Engine", "Peaks", "Reference", "annotate_peaks", "build_reference", "read_peak_file", "run_annotation",
           "run_enrichment", "run_peak_qc", "__version__"]
