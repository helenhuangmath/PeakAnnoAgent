"""MCP server exposing the same tools, so Claude Code (or any MCP client) can be the agent.

    claude mcp add peakanno -- python -m peakannoagent.mcp_server

The client calls start_run first, then the analysis tools, then write_report.
Requires `pip install 'peakannoagent[mcp]'` (MCP Python SDK 2.x).
"""

from __future__ import annotations

import asyncio
import json
from pathlib import Path

import mcp_types as types
from mcp.server.lowlevel import Server
from mcp.server.stdio import stdio_server

from .config import RunConfig, build_biology
from .report import build_report
from .tools import TOOLS, WRITE_REPORT_SCHEMA, RunContext, run_tool

START_RUN = {
    "name": "start_run",
    "description": "Start an analysis session: peak file, genome, output directory and the user's known biology. "
                   "Must be called before any other PeakAnnoAgent tool.",
    "input_schema": {
        "type": "object",
        "properties": {
            "peaks": {"type": "string"}, "genome": {"type": "string"}, "out_dir": {"type": "string"},
            "prefix": {"type": "string"}, "offline": {"type": "boolean"},
            "description": {"type": "string"}, "assay": {"type": "string"}, "target": {"type": "string"},
            "cell_type": {"type": "string"}, "condition": {"type": "string"},
            "expected_genes": {"type": "array", "items": {"type": "string"}},
            "expected_tfs": {"type": "array", "items": {"type": "string"}},
            "custom_gmt": {"type": "string"}, "rnaseq": {"type": "string"}, "blacklist": {"type": "string"},
            "rscript": {"type": "string"},
        },
        "required": ["peaks", "genome", "out_dir"],
    },
}

_session: dict = {}


def _schemas():
    return [START_RUN] + [t.schema() for t in TOOLS.values()] + [WRITE_REPORT_SCHEMA]


def _call(name: str, args: dict) -> dict:
    if name == "start_run":
        bio_keys = ("description", "assay", "target", "cell_type", "condition", "expected_genes", "expected_tfs",
                    "custom_gmt", "rnaseq")
        biology = build_biology({}, **{k: args.get(k) for k in bio_keys})
        extra = {"rscript": args["rscript"]} if args.get("rscript") else {}
        config = RunConfig(peaks=args["peaks"], genome=args["genome"], out_dir=args["out_dir"],
                           prefix=args.get("prefix", "run"), biology=biology, blacklist=args.get("blacklist"),
                           offline=bool(args.get("offline", False)), use_llm=True, model="mcp-client", **extra)
        _session["ctx"] = RunContext(config=config, log=lambda m: None)
        return {"status": "ok", "out_dir": config.out_dir, "genome": config.genome}
    ctx = _session.get("ctx")
    if ctx is None:
        return {"status": "error", "error": "Call start_run first."}
    if name == "write_report":
        path = str(Path(ctx.config.out_dir) / f"{ctx.config.prefix}_report.html")
        rep = build_report(ctx.evidence, ctx.config, args, path, state=ctx.state)
        ctx.evidence.save(str(Path(ctx.config.out_dir) / f"{ctx.config.prefix}_evidence.json"))
        return {"status": "ok", "report": rep["path"], "unverified_citations": rep["unverified_citations"]}
    ev = run_tool(ctx, name, args)
    return json.loads(ctx.evidence.compact(ev))


async def _list_tools(ctx, params) -> types.ListToolsResult:
    return types.ListToolsResult(tools=[types.Tool(name=s["name"], description=s["description"],
                                                   input_schema=s["input_schema"]) for s in _schemas()])


async def _call_tool(ctx, params) -> types.CallToolResult:
    try:
        result = await asyncio.to_thread(_call, params.name, dict(params.arguments or {}))
        is_error = result.get("status") == "error"
    except (KeyError, ValueError, FileNotFoundError) as err:
        result, is_error = {"status": "error", "error": str(err)}, True
    return types.CallToolResult(content=[types.TextContent(type="text", text=json.dumps(result, default=str))],
                                is_error=is_error)


server = Server("peakannoagent", on_list_tools=_list_tools, on_call_tool=_call_tool)


async def _main() -> None:
    async with stdio_server() as (read_stream, write_stream):
        await server.run(read_stream, write_stream, server.create_initialization_options())


def main() -> None:
    asyncio.run(_main())


if __name__ == "__main__":
    main()
