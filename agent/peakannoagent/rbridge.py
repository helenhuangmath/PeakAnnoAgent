"""Runs the R analysis engine (scripts/peak_tools.R) and returns its JSON outcome."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


def repo_root() -> Path:
    env = os.environ.get("PEAKANNOAGENT_HOME")
    if env:
        return Path(env).expanduser().resolve()
    return Path(__file__).resolve().parents[2]


class RBridge:
    def __init__(self, rscript: str = "Rscript", root: Path = None, timeout: int = 4 * 3600,
                 log_dir: str = None):
        self.rscript = rscript
        self.root = Path(root) if root else repo_root()
        self.script = self.root / "scripts" / "peak_tools.R"
        self.timeout = timeout
        self.log_dir = Path(log_dir) if log_dir else None

    def available(self) -> bool:
        return self.script.is_file() and (shutil.which(self.rscript) or Path(self.rscript).is_file())

    def run(self, tool: str, params: dict) -> dict:
        if not self.available():
            return {"status": "error", "error": f"Rscript '{self.rscript}' or {self.script} not found."}
        with tempfile.TemporaryDirectory(prefix="peakanno_") as tmp:
            params_path = Path(tmp) / "params.json"
            result_path = Path(tmp) / "result.json"
            params_path.write_text(json.dumps(params))
            try:
                proc = subprocess.run(
                    [self.rscript, str(self.script), tool, str(params_path), str(result_path)],
                    capture_output=True, text=True, timeout=self.timeout, cwd=str(self.root),
                )
            except subprocess.TimeoutExpired:
                return {"status": "error", "error": f"R tool '{tool}' timed out after {self.timeout}s."}
            if self.log_dir:
                self.log_dir.mkdir(parents=True, exist_ok=True)
                (self.log_dir / f"{tool}.log").write_text(proc.stdout + "\n--- stderr ---\n" + proc.stderr)
            if not result_path.is_file():
                tail = "\n".join(proc.stderr.strip().splitlines()[-15:])
                return {"status": "error", "error": f"R tool '{tool}' crashed (exit {proc.returncode}):\n{tail}"}
            return json.loads(result_path.read_text())
