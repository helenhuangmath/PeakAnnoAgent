\name{run_chippeakanno_pipeline}
\alias{run_chippeakanno_pipeline}
\title{Annotate Peaks and Derive Gene Sets}
\description{Run ChIPpeakAnno nearest-gene and promoter annotation, genomic region summaries,
TSS profiles and gene-set exports with a genome-wide gene universe.}
\usage{
run_chippeakanno_pipeline(peaks, out_dir, prefix, species = "mm10",
  txdb_pkg = NULL, orgdb_pkg = NULL, promoter_upstream = 2000,
  promoter_downstream = 500, distal_max_distance = 100000)
}
\arguments{
  \item{peaks}{A data.table with Chr, Start, End and PeakID. Coordinates are zero-based, half-open.}
  \item{out_dir}{Output directory.}
  \item{prefix}{Output file prefix.}
  \item{species}{Supported genome build, such as mm10 or hg38.}
  \item{txdb_pkg}{Optional TxDb package name overriding the genome default.}
  \item{orgdb_pkg}{Optional organism annotation package name.}
  \item{promoter_upstream, promoter_downstream}{Promoter extent in base pairs relative to the TSS.}
  \item{distal_max_distance}{Maximum absolute TSS distance for distal peak-to-gene links, in base pairs.}
}
\value{A named list of summaries, tables, figures, gene-set and universe paths, warnings and package versions.}
\details{The exported peak-gene table retains BED coordinates. Bioconductor range objects
are converted to one-based closed coordinates internally. Nearest-gene links are associations,
and do not establish enhancer targets.}
