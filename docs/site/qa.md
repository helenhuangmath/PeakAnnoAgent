# Common Q & A

## Running

### Do I need an API key?

No. `--no-llm` runs a fixed analysis plan and writes a factual report with no model involved. An
`ANTHROPIC_API_KEY` is needed only for agent mode (Claude plans and interprets) and for the chat
in the web app.

### Can it run without internet access?

Partly. `--offline` skips ChIP-Atlas, the literature search and the GEO search. Gene models and
gene-set libraries are still downloaded the first time they are needed, then cached. For a machine
with no internet at all, copy the cache folder from another machine, or pass your own files with
`--gtf` and `--gmt`.

### How long does a run take?

On the bundled example (351 peaks), without the model: about a second for annotation and a few
seconds for pathway enrichment and the networks, once the reference data are cached. The very
first run adds the downloads and about 6-10 s to build the gene-model reference. With network
access, ChIP-Atlas adds 5-30 minutes because the job runs on the DDBJ supercomputer; the agent
submits it early and collects the result later.

### Which assays does it support?

Any assay that yields peaks: transcription-factor and histone ChIP-seq, ATAC-seq, CUT&RUN,
CUT&Tag, DNase-seq. Telling it the assay (`--assay`) changes the analysis: broad histone marks skip
motif scanning, and enhancer assays get a distal-only motif scan.

## Annotation

### Which gene models should I use?

Use the ones the rest of your project uses, so gene names and counts agree with your RNA-seq. The
default is the GENCODE basic annotation; choose a release with `--reference gencode:<release>`, or
use `refseq`, `ucsc` or your own GTF. See [Gene models and the annotation method](execution.md).

### Will I get the same answer as ChIPseeker or ChIPpeakAnno?

Nearly. The method follows their conventions. When the same algorithm was compared with ChIPseeker
on the bundled example with UCSC gene models, it picked the same nearest gene for 99.7% of peaks
and the same feature type for 90%. The feature differences are definitions: "downstream" means
within 300 bp of a transcript end, a peak inside any gene's promoter window is a promoter peak even
when another gene's TSS is nearer, and transcripts with no gene ID are left out.

### Why is the nearest gene never a pseudogene or lncRNA?

Peaks are annotated against protein-coding genes by default. Use `--gene-types all` (or
`gene_types="all"` in Python and R) to annotate against every gene.

### Some peaks have no nearest gene. Why?

They lie on a sequence that the gene models do not contain, usually an unplaced contig or a
chromosome name in another style. The run reports how many. Peak QC removes non-standard contigs
before annotation.

### Does a nearest-gene link mean the peak regulates that gene?

No. It is an association. Distal peaks often regulate a gene that is not the nearest one; for
those, the ChIP-Atlas and RNA-seq evidence is more informative.

### Where is the reference cache, and how do I clear it?

In `~/.cache/peakannoagent`, or in the folder named by `PEAKANNOAGENT_CACHE`; `peakanno-agent check`
prints the location.
Delete the folder to clear it; references are rebuilt on the next run.

## Results

### Can I trust the hypotheses the agent writes?

Treat them as hypotheses. Each one cites the evidence IDs it rests on, and every citation is
checked against the results of the run; citations that could not be verified are listed in the
manifest. The statistics come from the analysis engine, not from the model.

### Are the results reproducible?

The deterministic mode gives the same tables and figures for the same input, package versions and
gene models; pathway-module clustering is seeded. The manifest records the input checksum, the
parameters and the package versions. In agent mode the analysis choices and the text can differ
between runs.

### How do I get the figures for a paper?

Every figure is written as a vector PDF with editable text under `results/` in the output folder.
The report embeds a picture of each figure and links to the PDFs. Change the colours with
`--palette` (default `chippeakanno`; also `npg`, `spectral` or a Matplotlib colormap name).

## Problems

### `peakanno-agent check` says a required package is missing

The engine needs NumPy and Matplotlib. Reinstall with `pip install -e .` in the environment that
`peakanno-agent` runs from.

### Motif enrichment is always skipped

It needs [HOMER](http://homer.ucsd.edu/homer/) with the genome package for your build
(`perl configureHomer.pl -install mm10`), and `findMotifsGenome.pl` on the `PATH`.

### A gene-model or gene-set download fails

Download the file yourself and pass it with `--gtf` (gene models) or `--gmt` (gene sets). The error
message names the URL that was tried.

### Where do I report a bug?

On [GitHub issues](https://github.com/helenhuangmath/PeakAnnoAgent/issues), with the version, the command, the genome and the relevant part
of the log. A small public or synthetic input helps.
