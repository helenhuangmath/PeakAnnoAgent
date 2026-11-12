# peakannoagent (Python)

The agent layer of PeakAnnoAgent. Give it a peak file, a genome, and what you already know about the
biology. It then:

1. QCs and cleans the peaks
2. annotates them with ChIPpeakAnno
3. tests pathways and gene sets
4. draws clustered pathway networks (SimpNetworkViz) and runs a cross-tool pathway consensus (PathwayAgent), when installed
5. finds enriched motifs
6. finds co-binding TFs in public ChIP-seq (ChIP-Atlas)
7. checks your positive controls
8. integrates RNA-seq
9. searches the literature and GEO

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
pip install -e "agent[web]"       # web app (Flask)
pip install -e "agent[mcp]"       # optional: MCP server (Python >= 3.10)
```

The core package needs only the Python standard library (Python 3.10+). Agent mode needs the Anthropic
SDK, which requires Python 3.10+, and Claude credentials (`ANTHROPIC_API_KEY`, or `ant auth login`).

If the Bioconductor packages live in a different R installation, point the agent at it:

```bash
export PEAKANNOAGENT_RSCRIPT=/path/to/your/Rscript
peakanno-agent check --genome mm10     # R packages, HOMER, credentials, network
```

**Companion tools (optional).** `network_viz` calls [SimpNetworkViz](https://github.com/helenhuangmath/SimpNetworkViz)
(`snv` on `PATH`, or `PEAKANNOAGENT_SNV`) on this run's enrichment table and writes the clustered pathway
network, a gene–pathway network, cluster labels, an interactive HTML viewer and Cytoscape files.
`pathway_consensus` calls [PathwayAgent](https://github.com/helenhuangmath/PathwayAgent) (`pathwayagent` on
`PATH`, a checkout next to this repository with its `.venv`, or `PEAKANNOAGENT_PATHWAYAGENT`) on the promoter,
distal and all peak-linked gene lists with the genome background, and reports consensus themes, tiers and key
genes with their figures. Offline runs use PathwayAgent's local engines only (`--local-only`, preset `quick`);
online runs default to preset `standard`. Both steps are skipped, not failed, when the program is missing.

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

## Web app

```bash
peakanno-agent web --runs-dir runs --port 8766 --open
```

| page | what it does |
|---|---|
| `/` | New-run form (peaks, genome, known biology, GMT, RNA-seq, mode) and the list of runs. The *Examples* gallery offers four peak sets (one synthetic, three from ENCODE); each card runs the example in one click (`use_example=<key>`) or fills the form with it (`/?example=<key>`). The *Test data* panel lists every bundled file (download at `/examples/<file>`). *Open an existing result folder* registers any `peakanno-agent run` output directory. |
| `/runs/<id>` | Live log while the run works; afterwards the run in parts with a sidebar: Overview, Summary, Hypotheses, Interpretation, Peaks & annotation, Pathways, Motifs, Regulators, Expression & controls, Literature & data, Search results, Ask Claude, Methods, Analysis log. `#pathways` or `?tab=search&q=Tox` opens a part directly. Rendered by `peakannoagent.webview` from the evidence store; for older runs without a manifest the narrative is recovered from the report. |
| `/api/runs/<id>/search?q=Tox&kind=gene` | JSON search over the run index. Kinds: gene, peak, pathway, module, motif, regulator, expression, control, paper, dataset, hypothesis, section, result. |
| `/api/runs/<id>/gene/<symbol>` | Everything the run knows about one gene. |
| `/api/runs/<id>/chat` | `GET` history, `POST {"message": …}` to ask, `DELETE` to clear. The transcript is saved as `<prefix>_chat.json` next to the run. |

The chat assistant (`peakannoagent.chat`) sees the run configuration, the known biology, the report
narrative and a compact view of every result, and has read-only tools: `search_run`, `gene_lookup`,
`inspect_table`, `search_literature`, `search_geo` (the last two are hidden for offline runs). It is held
to the same rule as the report: every claim cites an evidence id, PMID or accession from the run or from
a tool call in the conversation, and unknown citations are flagged in the answer.

Runs started from the app execute in a background thread of the server process; their status and log are
in `status.json` in the run directory, so the list survives a restart. The app is meant for a lab machine
or a laptop, not the public internet: it has no authentication and serves files from the run directories.

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

## Version 1.0.2 release checks

Python 3.10 or later is required. Wheels and source distributions bundle the R engine and small examples; R and its analysis dependencies must be installed separately.

Run outcomes are `success`, `partial`, or `failed`, with CLI exits 0, 3, or 1 (invalid invocation: 2). Failed core analyses retain partial outputs. `peakanno-agent check --offline` exits nonzero when required local dependencies are missing.

See [the release guide](../docs/release_1.0.2.md) for coordinate conventions, reproducible dependencies, package checks and biological validation.
