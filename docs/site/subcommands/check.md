# check

Check that everything a run needs is installed and reachable.

```bash
peakanno-agent check
peakanno-agent check --offline
```

## Options

| Option | Default | Description |
|---|---|---|
| `--genome` | `mm10` | Genome build to report on |
| `--offline` | off | Check local dependencies only, without contacting public services |

## What it reports

```text
Python 3.12.4
  ok       matplotlib 3.10.8
  ok       numpy 2.2.6
HOMER findMotifsGenome.pl: not found (motif enrichment is skipped)
Cache: /home/you/.cache/peakannoagent (2 gene-model reference(s) built)
SimpNetworkViz: not found (optional)
PathwayAgent: not found (optional)
anthropic SDK: not installed (needed for agent mode; --no-llm works without it)
ant CLI: not found; ANTHROPIC_API_KEY: not set
Europe PMC: reachable
NCBI E-utilities: reachable
ChIP-Atlas: reachable
```

The exit status is 0 when the required packages (NumPy and Matplotlib) are present and 1 otherwise.
HOMER, the companion tools and the API key are optional; steps that need them are skipped.
