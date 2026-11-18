\name{annotate_peaks}
\alias{annotate_peaks}
\title{Annotate Peaks}
\description{Annotate genomic peaks against a choice of gene models by running
\code{peakanno-agent annotate}, and return its tables as data frames.}
\usage{
annotate_peaks(peaks, genome, out_dir = tempfile("peakanno_"), prefix = "run",
  reference = "gencode", gtf = NULL, promoter_upstream = 2000,
  promoter_downstream = 500, distal_max_distance = 100000,
  gene_types = c("auto", "protein_coding", "all"), figures = TRUE,
  palette = "chippeakanno")
}
\arguments{
  \item{peaks}{Path of a BED file, or a data frame with \code{Chr}, \code{Start}, \code{End}
    (or chromosome, start and end as its first three columns) and optionally \code{PeakID}.
    Coordinates are zero-based and half-open, as in BED.}
  \item{genome}{Genome build of the peaks: \code{"mm10"}, \code{"mm39"}, \code{"hg38"} or \code{"hg19"}.}
  \item{out_dir}{Folder for the tables and figures.}
  \item{prefix}{Prefix of the output files.}
  \item{reference}{Gene models: \code{"gencode"} (default), \code{"gencode:<release>"} such as
    \code{"gencode:v44"}, \code{"refseq"} or \code{"ucsc"}.}
  \item{gtf}{Path of a GTF file to annotate against instead (GENCODE, Ensembl, RefSeq or custom).}
  \item{promoter_upstream, promoter_downstream}{Promoter window around a TSS, in base pairs.}
  \item{distal_max_distance}{Peaks further than this from any TSS are not linked to a gene.}
  \item{gene_types}{Genes to annotate against: protein-coding (\code{"auto"}, the default) or all.}
  \item{figures}{Draw the figures.}
  \item{palette}{Figure colours: \code{"chippeakanno"} (ChIPpeakAnno's default colours, the default),
    \code{"npg"} (Nature Publishing Group), \code{"spectral"}, or a Matplotlib
    colormap name such as \code{"Set2"}, \code{"Dark2"}, \code{"viridis"} or \code{"RdYlBu_r"}.}
}
\value{A list with \code{annotation} (one row per peak: feature, nearest gene, distance to TSS,
promoter genes, region class), \code{feature_distribution}, \code{region_distribution},
\code{tss_distance_bins}, \code{gene_sets}, \code{universe}, \code{figures} (paths of the PDF
figures), \code{summary} (the run summary) and \code{out_dir}.}
\details{The first use of a reference downloads the gene models and builds a cache; later runs
take about a second. The Python program must be installed: \code{pip install peakannoagent}.}
\seealso{\code{\link{run_analysis}}, \code{\link{paa_command}}}
\examples{
\dontrun{
res <- annotate_peaks("peaks.bed", genome = "hg38", reference = "gencode")
head(res$annotation)
res$feature_distribution
}
}
