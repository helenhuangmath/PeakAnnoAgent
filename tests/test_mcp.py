"""Exercise the optional MCP SDK through its in-process protocol client."""

import asyncio
import json
import tempfile
import unittest
from pathlib import Path
from importlib.util import find_spec

from helpers import EXAMPLE_PEAKS, FakeBridge

if find_spec("mcp") is not None:
    from mcp import Client
    from peakannoagent import mcp_server
else:
    Client = None


@unittest.skipIf(Client is None, "MCP extra not installed")
class MCPTest(unittest.TestCase):
    def test_tool_discovery_execution_and_failed_report_status(self):
        async def scenario(directory):
            async with Client(mcp_server.server) as client:
                tools = await client.list_tools()
                self.assertIn("start_run", [tool.name for tool in tools.tools])
                result = await client.call_tool("start_run", {"peaks": str(EXAMPLE_PEAKS), "genome": "mm10",
                                                              "out_dir": directory, "offline": True})
                self.assertFalse(result.is_error)
                mcp_server._session["ctx"].bridge = FakeBridge(directory)
                result = await client.call_tool("peak_qc", {})
                self.assertFalse(result.is_error)
                result = await client.call_tool("write_report", {"title": "Incomplete analysis",
                    "executive_summary": "QC only [E1].", "sections": [], "hypotheses": [], "caveats": []})
                self.assertTrue(result.is_error)
                self.assertEqual(json.loads(result.content[0].text)["status"], "failed")
                self.assertEqual(json.loads(Path(directory, "run_manifest.json").read_text())["status"], "failed")

        with tempfile.TemporaryDirectory() as directory:
            asyncio.run(scenario(directory))
