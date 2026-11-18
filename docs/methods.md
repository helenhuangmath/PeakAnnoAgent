# Methods

How each analysis step works. The Methods section of every report states the same, with the gene
models and gene-set versions used in that run.

## Peak QC

Peaks are read as zero-based, half-open intervals. QC counts duplicated intervals, peaks on
non-standard contigs and peaks that overlap a blacklist, reports the width distribution and peaks
per chromosome, and writes a cleaned peak set (standard chromosomes, no duplicates, no blacklisted
peaks) that later steps use.

## Gene models

| `reference` | Gene models |
|---|---|
| `gencode` (default) | GENCODE basic: hg38 v46, hg19 v46lift37, mm39 vM35, mm10 vM25 |
| `gencode:<release>` | a chosen GENCODE release |
| `refseq` | NCBI RefSeq (UCSC ncbiRefSeq GTF) |
| `ucsc` | UCSC knownGene (UCSC knownGene GTF) |
| `gtf = <file or URL>` | any GTF with `gene_id` and `transcript_id` |

The GTF is parsed once into sorted arrays of transcription start sites, exons, introns, UTRs (the
exonic parts outside each transcript's coding range) and transcript ends, and cached. Peaks are
annotated against protein-coding genes unless `gene_types` is `all`. When a GTF carries no gene
types, genes with a coding sequence are called protein-coding.

## Peak annotation

- **Feature.** Each peak gets the first that applies of: promoter (the peak overlaps the window
  from `promoter_upstream` bp upstream to `promoter_downstream` bp downstream of any TSS), 5' UTR,
  3' UTR, exon (first or other), intron (first or other), downstream (within 300 bp past a
  transcript end), distal intergenic.
- **Nearest gene.** The gene of the TSS nearest to the peak centre, with the signed distance
  (negative upstream of the gene).
- **Promoter genes.** Every gene with a TSS whose promoter window the peak overlaps.
- **Region class.** `promoter`; `distal` if the nearest TSS is within `distal_max_distance`
  (100 kb); otherwise `distal_far`, which is not linked to a gene.
- **Gene sets.** `promoter` (promoter genes), `distal` (nearest genes of distal peaks) and `all`
  (their union), with the annotated genes of the same gene models as the background.

Nearest-gene links are associations and do not establish enhancer targets.

## Pathway enrichment

Each gene set is tested for over-representation in each gene-set library with the hypergeometric
test, and p-values are adjusted with the Benjamini-Hochberg method within each library and gene
set. For the public libraries the background is the annotated genes that occur in the library, and
sets of 10-500 genes are tested; for your own GMT the background is all annotated genes and sets
of 3 genes or more are tested. Gene symbols are matched without regard to case.

| Method | Source |
|---|---|
| `go_bp`, `go_mf`, `go_cc`, `hallmark`, `reactome`, `wikipathways`, `tf_targets`, `cell_types` | MSigDB 2024.1, human and mouse collections |
| `immunesigdb` | MSigDB 2024.1 C7 ImmuneSigDB (human; matched to mouse genes by symbol) |
| `kegg` | MSigDB KEGG Medicus (human), Enrichr KEGG 2019 (mouse) |
| `custom` | your GMT file |

## Pathway networks

- **Enrichment map.** The 30 most significant terms are nodes. Two terms are linked when the
  Jaccard index of their hit genes is at least 0.25, keeping each term's three strongest links.
  Modules are found by greedy modularity clustering and named after their most significant term.
- **Gene-pathway network.** The lead term of each of the top five modules, linked to its genes;
  genes in two or more of those terms are reported as hubs.

## Motif enrichment

Known-motif enrichment uses HOMER `findMotifsGenome.pl` on 200 bp windows centred on the peaks
(at most 5,000 peaks), when HOMER and its genome package are installed. Otherwise the step is
skipped.

## Public data and literature

Co-binding with public ChIP-seq and ATAC-seq experiments uses ChIP-Atlas Enrichment Analysis
(random-permutation background, Fisher's exact test, Benjamini-Hochberg correction). Papers come
from Europe PMC and datasets from NCBI GEO.

## Reporting standard

Every tool result has an evidence ID. The report cites only those IDs and the PMIDs and accessions
retrieved in the run; any other citation is struck through. Failed and skipped steps are shown with
their reason.
