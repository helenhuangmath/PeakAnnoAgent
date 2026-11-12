\name{run_peak_qc}
\alias{run_peak_qc}
\title{Check and Clean Genomic Peak Intervals}
\description{Report peak widths, duplicates, chromosome distribution and blacklist overlaps,
and export a cleaned BED file.}
\usage{run_peak_qc(peaks, out_dir, prefix, blacklist = NULL)}
\arguments{
  \item{peaks}{A data.table with Chr, Start, End and unique PeakID columns. Coordinates are zero-based and half-open.}
  \item{out_dir}{Directory in which to write tables and figures.}
  \item{prefix}{Output file prefix.}
  \item{blacklist}{Optional path to a zero-based, half-open BED blacklist. Touching intervals do not overlap.}
}
\value{A named list with counts, width summaries, warnings, cleaned BED path, tables and figure paths.}
\details{Invalid coordinates, duplicate intervals, nonstandard contigs and blacklisted intervals
are removed. An error is raised if no valid peaks remain.}
