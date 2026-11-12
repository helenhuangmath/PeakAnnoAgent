"""Input fingerprints and runtime versions for reproducible run manifests."""

import hashlib
import platform
from pathlib import Path


def input_provenance(config):
    paths = {"peaks": config.peaks, "blacklist": config.blacklist,
             "custom_gmt": config.biology.custom_gmt, "rnaseq": config.biology.rnaseq}
    inputs = {}
    for name, value in paths.items():
        if value:
            path = Path(value)
            digest = hashlib.sha256()
            with path.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            inputs[name] = {"path": str(path), "bytes": path.stat().st_size, "sha256": digest.hexdigest()}
    return inputs


def runtime_versions(store):
    return {"python": platform.python_version(), "platform": platform.platform(),
            "tools": {ev.id: ev.versions for ev in store.items.values() if ev.versions}}
