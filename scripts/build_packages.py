"""Build R and Python distributions without repository dotfiles or testing files."""

import argparse
import hashlib
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="dist", help="Distribution output directory")
    parser.add_argument("--r", default=os.environ.get("PEAKANNOAGENT_R", "R"), help="R executable")
    parser.add_argument("--no-isolation", action="store_true", help="Use the installed Python build dependencies")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    out = Path(args.out).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    description = (root / "DESCRIPTION").read_text()
    version = re.search(r"^Version:\s*(\S+)", description, re.MULTILINE).group(1)

    # Stage only R package inputs: no ignore file is needed in the repository.
    with tempfile.TemporaryDirectory(prefix="peakanno-package-") as temporary:
        stage = Path(temporary) / "PeakAnnoAgent"
        stage.mkdir()
        for name in ("DESCRIPTION", "NAMESPACE", "LICENSE", "LICENSE.md", "README.md"):
            shutil.copy2(root / name, stage / name)
        for name in ("R", "man", "inst", "vignettes", "examples/data"):
            shutil.copytree(root / name, stage / name)
        subprocess.run([args.r, "CMD", "build", str(stage)], cwd=out, check=True)

    command = [sys.executable, "-m", "build", "--outdir", str(out)]
    if args.no_isolation:
        command.append("--no-isolation")
    subprocess.run(command + [str(root / "agent")], check=True)
    names = [f"PeakAnnoAgent_{version}.tar.gz", f"peakannoagent-{version}.tar.gz",
             f"peakannoagent-{version}-py3-none-any.whl"]
    checksums = []
    for name in names:
        digest = hashlib.sha256((out / name).read_bytes()).hexdigest()
        checksums.append(f"{digest}  {name}\n")
    (out / "SHA256SUMS").write_text("".join(checksums))
    print(f"Built R package, Python source/wheel and checksums in {out}")


if __name__ == "__main__":
    main()
