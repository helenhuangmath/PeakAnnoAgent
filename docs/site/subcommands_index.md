# Subcommands

`peakanno-agent` has five subcommands. `run` does the whole analysis and `annotate` does the
annotation only; the others start the web app, check the installation and list the agent's tools.

```text
peakanno-agent run      --peaks PEAKS --genome GENOME --out OUT [options]
peakanno-agent annotate --peaks PEAKS --genome GENOME --out OUT [options]
peakanno-agent web      [--runs-dir runs] [--port 8766] [--open]
peakanno-agent check    [--genome mm10] [--offline]
peakanno-agent tools
```

| Subcommand | Description |
|---|---|
| [`run`](peakannoagent/tool_usage.md) | Annotate peaks, find pathways, motifs and regulators, write the HTML report |
| [`annotate`](subcommands/annotate.md) | Annotate peaks only: tables and figures, in about a second |
| [`web`](subcommands/web.md) | Start the local web app |
| [`check`](subcommands/check.md) | Check packages, optional tools, the cache, credentials and network access |
| [`tools`](subcommands/tools.md) | List the analysis tools the agent can call |

```{toctree}
:maxdepth: 1

peakannoagent/tool_usage
subcommands/annotate
subcommands/web
subcommands/check
subcommands/tools
```
