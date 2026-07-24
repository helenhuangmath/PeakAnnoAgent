\name{PeakAnnoAgent}
\alias{PeakAnnoAgent}
\alias{PeakAnnoAgent-package}
\alias{arg_value}
\alias{default_enrichr_libraries}
\alias{export_annotation_web_inputs}
\alias{export_gene_lists}
\alias{export_metascape_input}
\alias{infer_species_defaults}
\alias{known_methods}
\alias{parse_cli_args}
\alias{peaks_to_granges}
\alias{plot_annotation_summary}
\alias{plot_enrichment_dotplot}
\alias{read_gene_table}
\alias{read_peak_file}
\alias{run_annotatr_annotation}
\alias{run_chipeakanno_annotation}
\alias{run_chipseeker_annotation}
\alias{run_clusterprofiler_go}
\alias{run_clusterprofiler_kegg}
\alias{run_david_enrichment}
\alias{run_enrichr}
\alias{run_gprofiler2}
\alias{run_homer_annotation}
\alias{run_msigdb_enrichment}
\alias{run_nearest_txdb_annotation}
\alias{run_reactomepa}
\alias{run_rgreat_annotation}
\alias{run_webgestaltr}
\alias{split_csv}
\alias{validate_methods}
\docType{package}
\title{PeakAnnoAgent: Peak Annotation and Pathway Analysis Agent}
\description{
PeakAnnoAgent provides reusable R functions and a command-line runner for
genomic peak annotation, pathway enrichment, and publication-ready figures.
}
\details{
The package is organized around optional method backends. Lightweight input,
export, plotting, and validation helpers can run with minimal dependencies.
Specialized annotation and enrichment methods require their corresponding
Bioconductor, command-line, or API-backed tools.

The command-line entrypoint is \code{scripts/peak_agent.R}. Public functions
are grouped into annotation wrappers, enrichment wrappers, figure builders, and
input/output helpers.
}
\author{
Hua Huang
}
\seealso{
\code{read_peak_file}, \code{read_gene_table}, \code{plot_annotation_summary},
\code{plot_enrichment_dotplot}
}
\keyword{package}
