# Gene models and the annotation method

Two things decide how peaks are annotated: **which gene models** they are compared with, and
**which genes** of those models are used. Both are set per run, for `run` and `annotate` alike.

## Gene models

| `--reference` | Gene models | Download on first use |
|---|---|---|
| `gencode` (default) | GENCODE basic annotation: hg38 v46, hg19 v46lift37, mm39 vM35, mm10 vM25 | 20-40 MB |
| `gencode:<release>` | a GENCODE release you choose, for example `gencode:v44` or `gencode:vM25` | 20-40 MB |
| `refseq` | NCBI RefSeq, from the UCSC ncbiRefSeq GTF | 20-40 MB |
| `ucsc` | UCSC knownGene, from the UCSC knownGene GTF | 20-40 MB |

`--gtf FILE` annotates against any GTF file instead (GENCODE, Ensembl, RefSeq or your own), and
overrides `--reference`. Chromosome names without the `chr` prefix, as in Ensembl files, are
converted.

```bash
peakanno-agent annotate --peaks peaks.bed --genome hg38 --out results
peakanno-agent annotate --peaks peaks.bed --genome mm10 --out results --reference gencode:vM25
peakanno-agent annotate --peaks peaks.bed --genome hg38 --out results --reference refseq
peakanno-agent annotate --peaks peaks.bed --genome hg38 --out results --gtf Homo_sapiens.GRCh38.112.gtf.gz
```

```{important}
The gene models must be on the same genome build as the peaks. hg38 takes any current GENCODE human
release and mm39 takes mouse vM26 or later. mm10 is frozen at vM25, and hg19 needs the lifted
"GRCh37" files (`v46lift37`).
```

### Which genes are used

Peaks are annotated against **protein-coding genes** by default. GENCODE and Ensembl list many
pseudogenes and lncRNAs, which would otherwise become the nearest gene of many peaks and enter the
enrichment background. `--gene-types all` uses every gene. When a GTF has no gene types, a gene
with a coding sequence counts as protein-coding.

### The cache

A reference is built the first time it is used (about 6 s for mouse GENCODE, about 10 s for human)
and read from a cache afterwards. Downloaded files and built references are stored in
`~/.cache/peakannoagent`, or in the folder named by `PEAKANNOAGENT_CACHE`. A built reference is
5-10 MB. Delete the folder to rebuild everything.

## How peaks are annotated

Gene models are parsed once into sorted arrays of transcription start sites, exons, introns, UTRs
and transcript ends. Annotating is then a series of binary searches over those arrays, which takes
well under a second for 100,000 peaks.

- **Feature.** The first that applies of promoter, 5' UTR, 3' UTR, exon (first or other), intron
  (first or other), downstream (within 300 bp of a transcript end) and distal intergenic. This is
  the priority ChIPseeker uses.
- **Promoter.** A peak that overlaps the promoter window of any TSS, per transcript. A peak at an
  alternative promoter is therefore a promoter peak, and both genes of a bidirectional promoter are
  recorded. The window is `--promoter-upstream` and `--promoter-downstream` (2000 and 500 bp by
  default); windows wider than 1 kb are reported in 1 kb distance bins.
- **Nearest gene.** The gene of the TSS nearest to the peak centre. The distance is negative
  upstream of the gene.
- **Region class.** `promoter`, `distal`, or `distal_far` for peaks further than
  `--distal-max-distance` (100 kb) from any TSS, which are not linked to a gene for enrichment.

### Agreement with ChIPseeker and ChIPpeakAnno

The method follows the conventions of the two Bioconductor packages. Before the R engine was
removed, the same algorithm was compared with them on the bundled example using UCSC gene models:
it picked the same nearest gene as ChIPseeker for 99.7% of peaks and the same feature type for 90%.
The feature differences are definitions: "downstream" means within 300 bp of a transcript end,
transcripts with no gene ID are left out, and a peak inside any gene's promoter window is a promoter
peak even when another gene's TSS is nearer.

## Figures

| Figure | Shows | Equivalent in |
|---|---|---|
| `feature_pie`, `feature_bar` | Share of peaks per feature type | `ChIPseeker::plotAnnoPie`, `plotAnnoBar` |
| `tss_distance_bins` | Peaks by distance to the nearest TSS, upstream and downstream | `ChIPseeker::plotDistToTSS` |
| `feature_upset` | Peaks by the combination of features they overlap | `ChIPseeker::upsetplot` |
| `chromosome_coverage` | Peak positions along chromosomes | `ChIPseeker::covplot` |
| `region_distribution` | Share of peaks per genomic region | `ChIPpeakAnno::assignChromosomeRegion` |
| `tss_profile` | Peak count in bins around TSSs | `ChIPpeakAnno::binOverFeature` |
| `tss_distance` | Histogram of distance to the nearest TSS | |

```{image} ../figures/feature_pie.png
:alt: Feature-type pie chart
:width: 620px
```
