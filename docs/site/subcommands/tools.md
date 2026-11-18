# tools

List the analysis tools the agent can call, with their descriptions.

```bash
peakanno-agent tools
```

Tools marked `[network]` contact a public service and are skipped with `--offline`. In agent mode
Claude chooses which tools to run and with which settings; with `--no-llm` a fixed plan runs them
in order.

The same list, with every argument, is on the [Agent tools](../reference/agent_tools.md) page.
