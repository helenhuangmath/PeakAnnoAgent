# Tutorial 2: ENCODE examples with a known answer

Three real ENCODE peak sets come with the repository. Each has a well-known expected result, so they are
a good way to check that your installation works and to see how the assay changes the analysis.

| Example | Genome | Peaks | ENCODE | Expected result |
|---|---|---|---|---|
| CTCF ChIP-seq, K562 | hg38 | 1000 strongest | ENCSR000BPJ, file ENCFF582SNT (IDR thresholded) | CTCF motif dominates; co-binding with cohesin (RAD21, SMC3) |
| H3K27ac ChIP-seq, human CD8 T cells | hg38 | 1000 strongest | ENCSR835OJV, file ENCFF283CAQ (replicated) | T-cell identity genes and pathways; RUNX, ETS, TCF/LEF, T-bet/Eomes regulators |
| H3K4me3 ChIP-seq, mouse BMDM | mm10 | 800 strongest | ENCSR916IRY, file ENCFF813FPW (replicated) | Nearly every peak on a promoter; macrophage identity genes |

Each example has a peak file and a biology file in `examples/data/`.

## CTCF in K562: a transcription factor

```bash
peakanno-agent run \
  --peaks examples/data/example_ctcf_k562_hg38.bed --genome hg38 \
  --biology-file examples/data/example_ctcf_k562_biology.json \
  --out results/ctcf --prefix ctcf
```

What to look for:

- **Motifs.** CTCF (and its paralogue CTCFL/BORIS) should be the top known motif, far ahead of the rest.
  Because the biology file sets `"target": "CTCF"`, the positive-control check confirms that the ChIPped
  factor is recovered.
- **ChIP-Atlas.** Public CTCF, RAD21 and SMC3 experiments should top the co-binding table.
- **Annotation.** Most peaks are distal or intergenic, as expected for insulator sites.

## H3K27ac in human CD8 T cells: an active-chromatin mark

```bash
peakanno-agent run \
  --peaks examples/data/example_h3k27ac_cd8_hg38.bed --genome hg38 \
  --biology-file examples/data/example_h3k27ac_cd8_biology.json \
  --out results/h3k27ac --prefix h3k27ac
```

What to look for:

- **Pathways.** T-cell activation and differentiation modules in the enrichment map, and CD8 T-cell
  signatures in ImmuneSigDB, which the agent adds because the cell type is an immune cell.
- **Positive controls.** CD8A, CD8B, CD3E, LCK and the other expected genes should carry peaks.
- **Regulators.** Lineage factors (RUNX3, ETS1, TCF7, EOMES) among the motif and ChIP-Atlas results.

## H3K4me3 in mouse macrophages: a promoter mark

```bash
peakanno-agent run \
  --peaks examples/data/example_h3k4me3_bmdm_mm10.bed --genome mm10 \
  --biology-file examples/data/example_h3k4me3_bmdm_biology.json \
  --out results/h3k4me3 --prefix h3k4me3
```

What to look for:

- **Annotation.** Almost every peak should be in the promoter class. If it is not, the promoter window or
  the genome is wrong.
- **Motif scanning is skipped or flagged.** Broad histone marks do not give meaningful known-motif
  results, so the agent does not rely on them.
- **Pathways.** Macrophage identity (Csf1r, Cd68, Lyz2, Spi1) together with housekeeping functions.

## Run them from the web app

All three appear on the home page of the [web app](web_app.md) under *Examples*. **Run** starts the
analysis in one click. **Edit first** fills the form, so you can change the settings before running.
