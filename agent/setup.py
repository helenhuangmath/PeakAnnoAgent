"""Include the R engine in wheels and source distributions built from agent/."""

from pathlib import Path
from shutil import copy2, copytree

from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.command.sdist import sdist


def bundle_engine(destination):
    root = Path(__file__).resolve().parent.parent
    if not (root / "R").is_dir():
        # A source distribution already carries the engine as package data.
        if not (Path(__file__).resolve().parent / "peakannoagent" / "_engine" / "R").is_dir():
            raise RuntimeError("R engine missing: build from the repository or the complete source distribution")
        return
    for directory in ("R", "examples/data"):
        copytree(root / directory, destination / directory, dirs_exist_ok=True)
    (destination / "scripts").mkdir(parents=True, exist_ok=True)
    for name in ("install_dependencies.R", "peak_agent.R", "peak_tools.R"):
        copy2(root / "scripts" / name, destination / "scripts" / name)
    for name in ("DESCRIPTION", "LICENSE", "LICENSE.md", "NAMESPACE"):
        copy2(root / name, destination / name)


class BuildWithEngine(build_py):
    def run(self):
        super().run()
        bundle_engine(Path(self.build_lib) / "peakannoagent" / "_engine")


class SourceWithEngine(sdist):
    def make_release_tree(self, base_dir, files):
        super().make_release_tree(base_dir, files)
        bundle_engine(Path(base_dir) / "peakannoagent" / "_engine")


setup(cmdclass={"build_py": BuildWithEngine, "sdist": SourceWithEngine})
