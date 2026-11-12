# Tutorial 3: The web app

The web app does everything the command line does, from a browser. You can also search a finished run
and ask Claude questions about it.

## Start the app

```bash
pip install -e "agent[web]"
peakanno-agent web --port 8766 --open          # http://127.0.0.1:8766, runs are written to ./runs
```

| Option | Default | Description |
|---|---|---|
| `--runs-dir` | `./runs` | Where runs are written |
| `--host`, `--port` | `127.0.0.1`, `8766` | Address to listen on |
| `--model` | `claude-opus-5-5` | Model for agent runs and the chat |
| `--rscript` | `Rscript` on `PATH` | Rscript with the Bioconductor packages |
| `--open` | off | Open the browser |

:::{warning}
The app has no authentication and serves files from the run folders. Run it on a lab machine or a
laptop, not on the public internet.
:::

## 1. Submit a run

The home page is a three-step form: **your peaks**, **what you know**, **how to run**. You can upload
a BED file, or pick one of the four examples. **Run** starts an example in one click. **Edit first**
fills the form with it so you can change it.

```{image} ../../figures/web/home_tabs.png
:alt: Home page with Run your analysis, Examples and Runs tabs
```

Run one of the bundled example peak sets locally to explore a finished report.

The home page separates **Run your analysis**, **Examples** and **Runs** into tabs. Use Examples to start from bundled peak sets, or Runs to reopen previous results.

## 2. Watch the run

While a run is in progress, the page streams the log. It updates itself when the results are ready.

## 3. Read the results

The run page has a sidebar with one part per topic:

```{image} ../../figures/web/run_overview.png
:alt: Run overview
```

| Part | What it shows |
|---|---|
| **Overview** | Key numbers, the biology you gave, and the agent's plan with its reasons |
| **Summary · Interpretation** | The agent's reading and key results by topic. Every claim carries an evidence badge that opens the matching result |
| **Suggestions** | Testable hypotheses with confidence, evidence and the experiment that would test each |
| **Peaks & annotation · Pathways · Motifs · Regulators · Expression & controls · Literature & data** | One card per analysis step, with statistics, figures (PNG and PDF), tables and warnings |
| **Search results** | Search across genes, peaks, pathways, motifs, regulators, papers, datasets and hypotheses |
| **Ask Claude** | A chat about this run, with cited evidence |
| **Methods · Analysis log** | Auto-written methods with package versions, and the full tool trace |

```{image} ../../figures/web/run_hypotheses.png
:alt: Hypotheses
```

## 4. Explore the pathways

The Pathways part opens with an interactive explorer:

- **Top terms**: a clickable bar chart per collection and peak group.
- **Pathway network**: the enrichment map with its modules. You can drag nodes and click a pathway
  to list its genes.
- **Gene–pathway map**: top terms and their genes. Click a gene to see its pathways, peaks and
  expression.

Every view downloads as SVG, PNG or TSV. See [Tutorial 4](pathway_networks.md).

```{image} ../../figures/web/run_pathway_network.png
:alt: Pathway explorer
```

## 5. Search and ask

Type a gene, a peak ID, a pathway or a motif into the search box. The search results group everything
the run knows about it:

```{image} ../../figures/web/run_search.png
:alt: Search for Batf
```

**Ask Claude** answers questions about the run. It can search the run, look up genes, read any table, and
query Europe PMC and GEO. Every claim cites an evidence ID, PMID or accession, and unverified
citations are flagged.

```{image} ../../figures/web/run_chat.png
:alt: Ask Claude about the run
```

## 6. Download

**Full HTML report** opens the self-contained report. It has the same pathway explorer, with the data
embedded, so it works offline. **Download all** gives a zip of the report, the evidence file and every
table and figure.
