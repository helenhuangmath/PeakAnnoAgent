# annotate

Annotate peaks and nothing else: no agent, no enrichment, no report. Once the gene models are
cached, it takes about a second.

```bash
peakanno-agent annotate --peaks peaks.bed --genome hg38 --out results/anno --prefix k562
```

It prints the feature breakdown and the path of a JSON summary:

```text
351 peaks annotated against GENCODE vM25 basic (mm10, 21,833 protein-coding genes) in 1.1 s
   20.23%  Promoter (<=1kb)
    1.99%  Promoter (1-2kb)
    1.71%  3' UTR
    4.27%  Other Exon
    7.41%  1st Intron
   17.09%  Other Intron
    0.28%  Downstream (<=300bp)
   47.01%  Distal Intergenic
results/anno/tex_annotation.json
```

## Options

| Option | Default | Description |
|---|---|---|
| `--peaks` | required | BED file, or a table with `Chr`, `Start`, `End` |
| `--genome` | required | `mm10`, `mm39`, `hg38` or `hg19` |
| `--out` | required | Output directory |
| `--prefix` | `run` | Prefix of every output file |
| `--reference` | `gencode` | `gencode`, `gencode:<release>`, `refseq` or `ucsc` (see [Gene models](../execution.md)) |
| `--gtf` | none | Annotate against this GTF file instead |
| `--palette` | `chippeakanno` | Figure colours: `chippeakanno` (ChIPpeakAnno's default colours, default), `npg` (Nature Publishing Group), `spectral`, or a Matplotlib colormap name such as `Set2`, `Dark2`, `tab10`, `viridis` or `RdYlBu_r` |
| `--promoter-upstream` | 2000 | Promoter window upstream of a TSS, bp |
| `--promoter-downstream` | 500 | Promoter window downstream of a TSS, bp |
| `--distal-max-distance` | 100000 | Peaks further than this from any TSS are not linked to a gene |
| `--gene-types` | `auto` | `auto` and `protein_coding`: protein-coding genes; `all`: every gene |
| `--no-figures` | off | Write the tables only |
| `--no-download` | off | Use cached gene models only; never download |

## Output

| File | Content |
|---|---|
| `<prefix>_peak_genes.tsv` | One row per peak: feature, features overlapped, nearest gene (symbol, ID, type, transcript), distance to TSS, promoter genes, region class |
| `<prefix>_feature_distribution.tsv` | Peaks per feature type |
| `<prefix>_region_distribution.tsv` | Peaks per genomic region |
| `<prefix>_tss_distance_bins.tsv` | Peaks per distance-to-TSS bin, upstream and downstream |
| `<prefix>_peak_gene_sets.tsv`, `<prefix>_gene_universe.tsv` | Gene sets (promoter, distal, all) and the background, for enrichment |
| `<prefix>_annotation.json` | The run summary: reference used, percentages, timings, warnings and file paths |
| `<prefix>_<figure>.pdf` | Eight figures as vector PDFs (see [Figures](../execution.md#figures)) |

The per-peak table looks like this:

```text
PeakID         Chr    Start     End       feature           nearest_symbol  distance_to_tss  promoter_symbols  region_class
tex_peak_0217  chr11  82045467  82045817  Other Intron      Taf15           18515            NA                distal
tex_peak_0218  chr11  83529612  83530012  Promoter (<=1kb)  Ccl5            28               Ccl5              promoter
tex_peak_0219  chr11  83581085  83581435  Distal Intergenic Ccl9            7875             NA                distal
```
