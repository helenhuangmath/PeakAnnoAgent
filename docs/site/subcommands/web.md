# web

Start the web app. It is a local page served from your own computer, not a hosted service.

```bash
peakanno-agent web --open
```

Then open <http://127.0.0.1:8766>. The page stays available while the command is running.

## Options

| Option | Default | Description |
|---|---|---|
| `--runs-dir` | `runs` | Folder in which each run gets its own directory |
| `--host` | `127.0.0.1` | Address to listen on |
| `--port` | `8766` | Port |
| `--model` | `claude-opus-5-5` | Default Claude model for agent runs and the chat |
| `--open` | off | Open the browser when the server starts |

## What the page does

- **Submit.** Upload a peak file or pick a bundled example, choose the genome, and describe the
  biology. Under **Gene annotation** you choose the reference gene models (GENCODE, RefSeq, UCSC
  or your own GTF).
- **Watch.** A live log shows each tool call and its status.
- **Read.** One card per analysis step, with figures, tables and downloads; the pathway explorer
  is interactive.
- **Search and ask.** Search genes, peaks, pathways and motifs across the run, and ask Claude
  questions that are answered with cited evidence.

Agent mode and the chat need `ANTHROPIC_API_KEY`. The deterministic mode needs nothing.
[Tutorial 3](../tutorials/web_app.md) walks through the page.
