# Sphinx configuration for the PeakAnnoAgent documentation site.
# Build locally:  pip install -r docs/site/requirements.txt
#                 sphinx-build -b html docs/site docs/site/_build/html
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

project = "PeakAnnoAgent"
author = "Hua Huang"
copyright = "2026, Hua Huang"
_desc = (ROOT / "DESCRIPTION").read_text()
release = re.search(r"^Version:\s*(\S+)", _desc, re.M).group(1)
version = release

extensions = ["myst_parser"]
myst_enable_extensions = ["colon_fence", "deflist", "attrs_inline"]
myst_heading_anchors = 3
source_suffix = {".md": "markdown"}
master_doc = "index"
exclude_patterns = ["_build", "requirements.txt"]
suppress_warnings = ["myst.xref_missing", "myst.header"]

html_theme = "sphinx_rtd_theme"
html_title = f"PeakAnnoAgent {release}"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "navigation_depth": 3,
    "collapse_navigation": False,
    "style_nav_header_background": "#2c3e50",
}
html_context = {
    "display_github": True,
    "github_user": "helenhuangmath",
    "github_repo": "PeakAnnoAgent",
    "github_version": "main",
    "conf_py_path": "/docs/site/",
}


# Regenerate reference/agent_tools.md from the agent's tool definitions, so the page follows
# agent/peakannoagent/tools.py. The committed copy is kept if the package cannot be imported.
def _write_agent_tools():
    import sys

    sys.path.insert(0, str(ROOT / "agent"))
    try:
        from peakannoagent.tools import TOOLS
    except Exception as exc:  # pragma: no cover - docs build only
        print(f"[conf.py] agent tools page not regenerated: {exc}")
        return
    out = [
        "# Agent tools\n",
        "These are the tools the agent calls. In agent mode, Claude chooses which ones to run and with which",
        "settings. With `--no-llm`, a fixed plan runs them in order. Tools marked *network* contact a",
        "public service and are skipped with `--offline`. `peakanno-agent tools` prints the same list.\n",
        "Every tool runs in the R engine or calls a public API. Its result is stored with an evidence ID",
        "(`E1`, `E2`, …), which the report cites.\n",
    ]
    for t in TOOLS.values():
        out.append(f"## `{t.name}`" + (" *(network)*" if t.needs_network else "") + "\n")
        out.append(t.description + "\n")
        if t.properties:
            out.append("| Argument | Type | Description |\n|---|---|---|")
            for k, v in t.properties.items():
                typ = v.get("type", "")
                if typ == "array" and "items" in v:
                    typ = "array of " + v["items"].get("type", "")
                    if "enum" in v["items"]:
                        typ += ": " + ", ".join(f"`{x}`" for x in v["items"]["enum"])
                if "enum" in v:
                    typ += ": " + ", ".join(f"`{x}`" for x in v["enum"])
                req = " (required)" if k in t.required else ""
                out.append(f"| `{k}`{req} | {typ} | {v.get('description', '').replace('|', '/')} |")
            out.append("")
    text = "\n".join(out) + "\n"
    page = Path(__file__).parent / "reference" / "agent_tools.md"
    if not page.exists() or page.read_text() != text:
        page.write_text(text)


_write_agent_tools()
