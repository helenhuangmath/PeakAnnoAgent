#!/usr/bin/env bash
# Fill sha256 in meta.yaml from the published GitHub release tarball.
# Run after the release tag (v<version>) exists on GitHub:  packaging/bioconda/update_sha256.sh
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
recipe="$here/meta.yaml"
version="$(sed -n 's/^{% set version = "\(.*\)" %}$/\1/p' "$recipe")"
url="https://github.com/helenhuangmath/PeakAnnoAgent/archive/refs/tags/v${version}.tar.gz"
tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT
if ! curl -fsSL "$url" -o "$tmp"; then
  echo "Could not download $url - publish the v${version} release on GitHub first." >&2
  exit 1
fi
sha="$(shasum -a 256 "$tmp" | cut -d' ' -f1)"
sed -i.bak -E "s/^(  sha256: ).*/\1${sha}/" "$recipe" && rm -f "$recipe.bak"
echo "v${version}: sha256 ${sha} written to $recipe"
