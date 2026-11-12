\name{agent_tool_names}
\alias{agent_tool_names}
\alias{run_agent_tool}
\alias{supported_genomes}
\title{Discover and Run Analysis Tools}
\description{Discover supported genome builds and execute an analysis tool using a named parameter list.}
\usage{
agent_tool_names()
supported_genomes()
run_agent_tool(name, params)
}
\arguments{
  \item{name}{A tool name returned by \code{agent_tool_names()}.}
  \item{params}{Named list. Analysis tools require \code{peaks} (BED or TSV path),
  \code{out_dir}, \code{prefix}, and \code{genome}. Pathway enrichment instead takes
  \code{gene_sets}, \code{universe}, and \code{methods}. Motif analysis accepts
  \code{backend}, \code{size}, and \code{max_peaks}. See the command-line reference
  in the repository for complete tool parameters.}
}
\value{Discovery functions return character vectors. \code{run_agent_tool} returns a named
list of summary statistics, output paths, warnings and software versions; tool errors raise an R error.}
\details{All input BED and Chr/Start/End TSV intervals use zero-based, half-open
coordinates. Conversion to one-based closed ranges occurs at \code{peaks_to_granges}.
Bioconductor dependencies are loaded only when the selected tool needs them.}
\examples{
agent_tool_names()
supported_genomes()
}
