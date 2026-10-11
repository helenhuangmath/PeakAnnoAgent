# Agent Design

PeakAnnoAgent has two layers. Claude plans and interprets. A deterministic R engine computes.

```
User: peaks + genome + known biology (+ GMT signatures, RNA-seq DE table)
        │
   ┌────▼──────────────────────────────────────────────┐
   │ Agent loop (agent/peakannoagent/agent.py)         │
   │ Claude: plan → call tool → read result → adapt    │
   │ deterministic fallback: run_pipeline()            │
   └────┬──────────────────────────────────────────────┘
        │ tools (agent/peakannoagent/tools.py), each result → evidence id E1, E2, …
        ├── R engine via scripts/peak_tools.R (JSON in/out)
        │     peak_qc · annotate_peaks (ChIPpeakAnno) · enrich_pathways (clusterProfiler/MSigDB/custom)
        │     great_enrichment (rGREAT) · find_motifs (HOMER or monaLisa+JASPAR)
        ├── Python: integrate_rnaseq · check_expectations · inspect_table
        └── Network: chipatlas_submit/collect · chipatlas_target_genes · search_literature (Europe PMC)
                     · search_geo (NCBI)
        │
   report.py → <prefix>_report.html (citations verified against the evidence store)
```

## Why this split

- **Reproducible statistics.** Every number comes from R or Python code with recorded package versions.
  The model chooses parameters and interprets the results; it never does arithmetic.
- **One tool contract, three front ends.** The same tool registry backs:
  - the Claude API loop (`peakanno-agent run`)
  - the fixed pipeline (`--no-llm`)
  - an MCP server (`python -m peakannoagent.mcp_server`), so Claude Code or another MCP client can
    drive the analysis interactively
- **The engine stays a normal R package.** `scripts/peak_agent.R` still works on its own for HPC runs.

## Agent loop

- The loop uses the Claude Messages API: model `claude-opus-5-5`, adaptive thinking, effort `high`, and
  server-side refusal fallback.
- The system prompt is stable. Run-specific details and the biology go in the first user turn, so the
  prompt cache stays warm across turns.
- Every tool schema has a `rationale` field. The model's reason for each step is logged in the report's
  analysis log.
- Parallel tool calls run in order, and their results return in a single user message.
- `write_report` ends the run. Unknown citations are returned to the model once for correction.
- Turn budget: `--max-turns`, default 40. Near the limit the model is told to finish.

## Extending

To add a tool:

1. If it is R-backed, add a `tool_<name>(params)` function to `R/agent_tools.R`, and add the name to
   `agent_tool_names()`.
2. Add a `Tool(...)` entry with a handler in `agent/peakannoagent/tools.py`.
3. Optionally, give it a section in `report.TOOL_SECTIONS` and a sentence in `report.methods_text`.

The handler must return `_record(...)`, so the result gets an evidence id. The agent, the MCP server, and
the report then pick up the new tool automatically.

## Ideas not built yet

- Peak-to-gene linking beyond nearest gene: ABC / ENCODE-rE2G, Hi-C/TADs, eQTLs
- Differential binding (DiffBind/csaw) for multi-condition designs
- Variant layers (GWAS Catalog overlap, motifbreakR, AlphaGenome) for human peaks
- Footprinting (TOBIAS) for ATAC-seq when BAMs are available
- Shiny front end that calls the same CLI
