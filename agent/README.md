# peakannoagent (Python)

The agent layer of PeakAnnoAgent. Give it a peak file, a genome, and what you already know about the
biology. It then:

1. QCs and cleans the peaks
2. annotates them with ChIPpeakAnno
3. tests pathways and gene sets
4. finds enriched motifs
5. finds co-binding TFs in public ChIP-seq (ChIP-Atlas)
6. checks your positive controls
7. integrates RNA-seq
8. searches the literature and GEO

Claude chooses the steps and settings and interprets the results. The output is an HTML report in which
every claim links to the result it came from.

The R package in the repository root is the analysis engine. This package runs it as tools through
`scripts/peak_tools.R`.

## Install

```bash
# R engine (once)
Rscript scripts/install_dependencies.R

# Python agent
pip install -e "agent[agent]"     # Anthropic SDK for agent mode
pip install -e "agent[mcp]"       # optional: MCP server (Python >= 3.10)
```

The core package needs only the Python standard library (Python 3.9+). Agent mode needs the Anthropic
SDK, which requires Python 3.10+, and Claude credentials (`ANTHROPIC_API_KEY`, or `ant auth login`).

If the Bioconductor packages live in a different R installation, point the agent at it:

```bash
export PEAKANNOAGENT_RSCRIPT=/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/bin/Rscript
peakanno-agent check --genome mm10     # R packages, HOMER, credentials, network
```

## Run

```bash
peakanno-agent run \
  --peaks examples/data/example_tex_peaks_mm10.bed --genome mm10 --out results/tex --prefix tex \
  --biology-file examples/data/example_tex_biology.json \
  --gmt examples/data/example_tex_signatures.gmt \
  --rnaseq examples/data/example_tex_rnaseq_de.tsv
```

Open `results/tex/tex_report.html`.

**Known biology.** Give it with flags (`--assay`, `--target`, `--cell-type`, `--condition`,
`--expected-genes`, `--expected-tfs`, `--biology "free text"`), or with a JSON file of the same fields, or
with a markdown notes file. It changes how the run behaves:
- The promoter window adapts to the assay.
- Motif scans are skipped for broad histone marks.
- Enhancer assays get a distal-only motif scan.
- Gene-set collections are chosen to fit the cell type.
- Expected genes and TFs, and the ChIPped TF itself, become positive controls.
- Literature queries are built from your system.

**Modes**

| flag | what happens |
|---|---|
| (default) | Claude (`claude-opus-5-5`, adaptive thinking) plans and interprets |
| `--no-llm` | fixed analysis plan, automated factual summary, no API key needed |
| `--offline` | no ChIP-Atlas / Europe PMC / GEO calls |

ChIP-Atlas jobs run on the DDBJ supercomputer and usually take 5–30 minutes. The agent submits early and
collects later. `--chipatlas-wait` sets how long the fixed plan waits.

## Outputs

| file | content |
|---|---|
| `<prefix>_report.html` | Report. Includes: summary, hypotheses with next experiments, interpretation, every result with figures and table previews, an igv.js browser, references, an auto-generated Methods section with package versions, and the analysis log |
| `results/` | All TSV tables, plus Nature-style vector PDFs (Helvetica, editable text) and PNG previews of every figure |
| `<prefix>_evidence.json` | Every tool result, papers, datasets, and the step trace |
| `<prefix>_manifest.json` | Config, timing, LLM token usage, and any unverified citations |

## Guardrails

- **Statistics come only from the tools.** The model never computes them.
- **Citations are checked.** Every tool result gets an id (`E1`, `E2`, …). The report must cite these
  ids, PMIDs, or accessions that were actually retrieved in this run. Unknown citations are bounced back
  to the model once; any that remain are struck through in red in the report.
- **Untrusted text stays data.** User notes are wrapped in `<user_notes>` tags. Literature and GEO text
  are passed as data, not as instructions.
- **Failures degrade gracefully.** If the model call fails or is declined, the run falls back to the fixed
  plan. Each failed tool is recorded and reported, not fatal.

## Use from Claude Code (MCP)

```bash
claude mcp add peakanno -- python -m peakannoagent.mcp_server
```

Claude Code then becomes the agent. It calls `start_run`, then the analysis tools, then `write_report`.

## Tests

```bash
cd agent && python -m unittest discover -s tests     # Python (uses the R bridge if Rscript is available)
Rscript tests/run_tests.R                            # R engine, from the repo root
```
