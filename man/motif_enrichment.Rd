\name{run_motif_enrichment}
\alias{run_motif_enrichment}
\alias{run_homer_motifs}
\alias{run_monalisa_motifs}
\alias{read_homer_known_results}
\title{Known Motif Enrichment}
\description{Find enriched known transcription-factor motifs using HOMER or monaLisa with JASPAR matrices.}
\usage{
run_motif_enrichment(peaks, out_dir, prefix, species = "mm10",
  backend = "auto", size = 200, max_peaks = 5000,
  padj_cutoff = 0.05, top = 25, ...)
run_homer_motifs(peaks, out_dir, prefix, species = "mm10", size = 200,
  max_peaks = 5000, genome = NULL, threads = 4, denovo = FALSE)
run_monalisa_motifs(peaks, out_dir, prefix, species = "mm10", size = 200,
  max_peaks = 5000, bsgenome_pkg = NULL)
read_homer_known_results(path)
}
\arguments{
  \item{peaks}{A data.table of zero-based, half-open peaks with Chr, Start, End and PeakID.}
  \item{out_dir}{Output directory.}
  \item{prefix}{Output file prefix.}
  \item{species}{Supported genome build.}
  \item{backend}{One of auto, homer or monalisa. Auto uses HOMER when available.}
  \item{size}{Sequence window width, in base pairs.}
  \item{max_peaks}{Maximum number of peaks to scan.}
  \item{padj_cutoff}{Adjusted p-value threshold.}
  \item{top}{Maximum number of significant motifs in the summary.}
  \item{...}{Additional HOMER backend arguments.}
  \item{genome}{Optional HOMER genome identifier.}
  \item{threads}{HOMER thread count.}
  \item{denovo}{Whether HOMER should also discover de novo motifs.}
  \item{bsgenome_pkg}{Optional BSgenome package name.}
  \item{path}{Path to HOMER knownResults.txt.}
}
\value{The main wrapper returns summaries, output paths, warnings and versions. Backend
wrappers return a standardized motif table with backend metadata. The reader returns a data.table
of motif names, TFs, family, target/background percentages, enrichment and p-values.}
\details{monaLisa uses a GC-matched genomic background. Motif enrichment indicates compatible
sequence specificity and does not prove factor occupancy. Broad histone domains should not use this scan.}
