# Methods Supported by the Framework

This repository separates analysis methods into interchangeable backends.

## Peak annotation

- `chipseeker`: local R/Bioconductor annotation with `ChIPseeker`, `TxDb.*`, and `org.*.eg.db`.
- `nearest`: lightweight nearest-gene assignment using `TxDb` genes and `GenomicRanges::distanceToNearest`.
- `chipeakanno`: annotation using `ChIPpeakAnno::annotatePeakInBatch`.
- `annotatr`: region annotation using the `annotatr` annotation catalog.
- `rgreat`: submits regions to GREAT through `rGREAT`.
- `homer`: calls HOMER `annotatePeaks.pl` from the shell. HOMER must already be installed and available in `PATH`.
- `webexport`: writes clean BED files for web tools.

## Pathway enrichment

- `clusterprofiler`: local GO enrichment with an explicit universe, recommended for publishable analyses when your tested background is known.
- `kegg`: KEGG ORA through `clusterProfiler::enrichKEGG`.
- `reactome`: Reactome pathway ORA through `ReactomePA::enrichPathway`.
- `msigdb`: MSigDB enrichment through `msigdbr` and `clusterProfiler::enricher`.
- `enrichr`: API-backed enrichment across GO, KEGG, Reactome, WikiPathways, MSigDB Hallmark, and TRRUST libraries.
- `gprofiler`: API-backed enrichment through `gprofiler2`.
- `webgestalt`: API-backed ORA through `WebGestaltR`.
- `metascape`: exports grouped gene lists in a Metascape-friendly table for upload or downstream API automation.
- `david`: uses `RDAVIDWebService` when installed and `DAVID_EMAIL` is configured.
- `export`: writes long and wide gene-list tables for manual upload tools.

## Recommended reporting standard

For manuscripts, run at least one local method with an explicit universe, such as `clusterprofiler`, and use API methods as complementary evidence. Save the exact input gene sets, background gene universe, package versions, and generated enrichment tables with every figure.

Figures follow Nature's artwork guidelines (`R/figures.R`). Every figure is written as a vector PDF with Helvetica text that stays editable, and with no Dingbats glyphs. Text is black, axes are thin and L-shaped, and there are no gridlines. ggplot figures use an 18 pt base font, and gene symbols are set in italics. Bar plots are dark grey. Significance colours run from light to dark red (#990000), with Penn Blue (#011F5B) for reference lines. A PNG copy of each figure is embedded in the HTML report, which links each one to its PDF.

## Agent analyses (agent/ + R engine)

- `peak_qc`: standard chromosomes, duplicates, width summary, optional blacklist removal; writes a cleaned BED.
- `annotate_peaks`: ChIPpeakAnno `annotatePeakInBatch` (nearest gene by TSS, plus all genes whose TSS window overlaps the peak), `assignChromosomeRegion`, `binOverFeature`; builds promoter / distal / all gene sets and a genome background.
- `enrich_pathways`: clusterProfiler ORA (GO, KEGG, ReactomePA, any MSigDB collection) against the genome background; user GMT signatures use a hypergeometric test against the same full background.
- `find_motifs`: HOMER `findMotifsGenome.pl` (known motifs) when installed, otherwise monaLisa with JASPAR CORE vertebrate PWMs vs a GC-matched random genomic background.
- `chipatlas_*`: ChIP-Atlas Enrichment Analysis (overlap with all public ChIP/ATAC experiments vs permuted background) and ChIP-Atlas target-gene tables.
- `integrate_rnaseq`, `check_expectations`: hypergeometric tests of DE genes vs peak-linked genes; positive-control recovery of expected genes and TFs.
- `search_literature`, `search_geo`: Europe PMC and NCBI GEO searches; retrieved PMIDs/accessions are the only citations the report accepts.
