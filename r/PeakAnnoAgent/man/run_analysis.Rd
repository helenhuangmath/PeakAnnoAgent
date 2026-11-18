\name{run_analysis}
\alias{run_analysis}
\title{Run the Full Analysis}
\description{Run \code{peakanno-agent run}: peak QC, annotation, pathway enrichment and networks,
motif enrichment, optional public-data and literature lookups, and the HTML report.}
\usage{
run_analysis(peaks, genome, out_dir, prefix = "run", llm = FALSE, offline = FALSE,
  reference = "gencode", gtf = NULL, blacklist = NULL, assay = NULL, target = NULL,
  cell_type = NULL, condition = NULL, biology = NULL, expected_genes = NULL,
  expected_tfs = NULL, gmt = NULL, rnaseq = NULL, palette = "chippeakanno",
  echo = interactive())
}
\arguments{
  \item{peaks, genome, out_dir, prefix, reference, gtf, palette}{As in \code{\link{annotate_peaks}}.}
  \item{llm}{Let Claude plan the analysis and interpret the results. Needs the Anthropic SDK and
    \code{ANTHROPIC_API_KEY}. The default is the fixed, deterministic plan.}
  \item{offline}{Skip ChIP-Atlas, literature and GEO lookups.}
  \item{blacklist}{Blacklist BED file; overlapping peaks are removed.}
  \item{assay, target, cell_type, condition, biology}{What is known about the experiment.}
  \item{expected_genes, expected_tfs}{Character vectors of positive controls.}
  \item{gmt}{GMT file of your own gene signatures.}
  \item{rnaseq}{Differential-expression table to integrate.}
  \item{echo}{Show the progress log.}
}
\value{A list with \code{status} (\code{"success"} or \code{"partial"}), \code{reasons},
\code{report} (path of the HTML report), \code{annotation} and \code{enrichment} data frames,
the paths of all \code{tables} and \code{figures}, the \code{manifest} and \code{out_dir}.}
\examples{
\dontrun{
res <- run_analysis("peaks.bed", genome = "mm10", out_dir = "results/tex",
  assay = "ATAC-seq", cell_type = "CD8 T cells")
browseURL(res$report)
}
}
