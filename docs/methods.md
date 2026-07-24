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
