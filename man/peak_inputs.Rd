\name{read_peak_file}
\alias{read_peak_file}
\alias{peaks_to_granges}
\title{Read BED Peaks and Convert Coordinate Systems}
\description{Read genomic peak intervals without altering their BED coordinates,
then explicitly convert them at the Bioconductor boundary.}
\usage{
read_peak_file(path)
peaks_to_granges(peaks)
}
\arguments{
  \item{path}{BED file or headered TSV containing Chr, Start and End. Both use
  zero-based, half-open intervals. A BED name column is used as PeakID when valid.}
  \item{peaks}{A table containing Chr, Start and End, with optional PeakID.}
}
\value{The reader returns a data.table retaining zero-based, half-open coordinates.
The converter returns a GRanges with one-based, closed coordinates and PeakID names.}
\details{The BED interval chr1:0-1 becomes the single R base chr1:1-1. Intervals
chr1:0-1 and chr1:1-2 do not overlap. The input table is not modified by conversion.
Fractional coordinates and integer overflow are rejected on read. Invalid intervals
are removed by run_peak_qc; peaks_to_granges rejects any remaining invalid intervals.
Peak-gene tables and BED exports retain the input convention. Raw backend range
tables returned by ChIPseeker, ChIPpeakAnno and annotatr retain their native R convention.}
