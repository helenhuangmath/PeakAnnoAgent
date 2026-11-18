# PeakAnnoAgent 1.0.2 — update-network-plots

Release version: **1.0.2**. Git tag: **v1.0.2**. Release title: **update-network-plots**.

Prepared November 12, 2026. This document distinguishes implemented fixes, local
verification, and checks requiring external infrastructure. Preparing these files
does not publish a GitHub release or upload packages to a package index.

## Local validation results

| Check | Result | Scope |
| --- | --- | --- |
| Python suite | 62 tests, 1 skipped, no failures | Python 3.12; includes real R QC bridge, scripted model clients, web regression and MCP protocol checks |
| R source suite | Passed | R 4.3.2; includes GenomicRanges boundary tests |
| R package build | Passed | Vignette rendered |
| Installed R package check | 0 errors, 0 warnings, 1 note | Optional suggested packages unavailable: annotatr, gprofiler2, JASPAR2024, rGREAT, RDAVIDWebService, WebGestaltR |
| Python wheel/source build | Passed | R engine, examples and web assets bundled |
| Installed-wheel smoke | Passed | Isolated environment, outside checkout; web home renders |
| Public endpoint connectivity | Passed | Europe PMC query, NCBI E-utilities and ChIP-Atlas status; system CA bundle configured |
| mm10 ATAC workflow | Success | 351 peaks, custom GMT, RNA-seq, pathways and all/distal motifs |
| hg38 CTCF workflow | Success | 1,000 peaks; CTCF highest-ranked motif |

Artifacts and logs are written under `release-artifacts/1.0.2/` (ignored by Git).
The validation R installation used JASPAR2020 as the supported fallback.
The Docker daemon was unavailable, and no live Anthropic credentials were set;
Linux container execution and live model behavior remain external checks.
The isolated Python installation had no default CA bundle. Endpoint checks
passed with `SSL_CERT_FILE=/etc/ssl/cert.pem`; certificate verification remained enabled.

## Completion contract

New manifests include `version`, `schema_version`, `status`, `required_steps`,
`missing_steps`, `unsuccessful_steps`, `reasons`, input hashes, runtime versions
and the individual tool trace.

| Outcome | Meaning | CLI exit | Web label |
| --- | --- | --- | --- |
| success | Required analyses completed; no unresolved tool errors | 0 | done |
| partial | QC/annotation completed, but another required step or attempted tool failed, or requested interpretation was unavailable | 3 | partial |
| failed | QC or annotation failed or was not run | 1 | failed |
| invalid invocation | Missing inputs, invalid parameters or a missing requested Python extra | 2 | error, where applicable |

Required analyses are QC, annotation and pathways, plus motifs for assays other
than broad histone domains and RNA-seq integration when supplied. Unavailable
optional tools are recorded as skipped. A legitimate zero-hit enrichment or
motif result is valid; it is distinguished from a dependency failure. Failure of
one pathway method is recorded even when other methods produce results.

Historical runs lacking a completion status are not automatically certified as
successful when reopened. Interrupted jobs retain that designation unless a new
manifest records a completed outcome. Partial outputs remain downloadable.

## Coordinates and compatibility

BED files and headered `Chr/Start/End` TSV inputs both use **zero-based, half-open**
coordinates. `chr1 0 1` is a single base. Conversion to one-based, closed R ranges
happens once in `peaks_to_granges`; blacklist overlap uses equivalent closed
copies. Cleaned BED and agent peak-gene tables retain BED coordinates. Raw
backend range tables retain their native R convention, as documented in R help.

Re-run analyses produced before 1.0.2 before using quantitative annotation or
overlap results. Existing figures and historical demo reports illustrate the
interface; their presence does not establish validation against corrected code.

Python requires 3.10 or later. The light R package requires R 4.1 or later;
the complete container recipe uses Bioconductor 3.20 with its corresponding R
version. Older Bioconductor installations may need the supported JASPAR2020/2022
fallback rather than JASPAR2024. The web app supports a single local user and
binds to localhost by default. Restart recovery preserves results; it does not
resume interrupted analyses. Multi-user hosted deployment is outside this
release's validation scope.

## Reproduce the checks

Use a dedicated Python environment. `constraints-release.txt` records resolved
validation dependencies, including transitive dependencies. Install the project
explicitly with these constraints:

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -c constraints-release.txt build "./agent[web,agent,mcp]"
.venv/bin/python -m unittest discover -s agent/tests
.venv/bin/python -m build --outdir release-artifacts/1.0.2 agent
.venv/bin/python -m pip install --force-reinstall --no-deps release-artifacts/1.0.2/*.whl
.venv/bin/python scripts/smoke_install.py
```

The default suite uses scripted model clients. Enable real R bridge checks with
`PEAKANNOAGENT_R_INTEGRATION=1` and set `PEAKANNOAGENT_RSCRIPT` to the intended
Rscript executable. Set `PEAKANNOAGENT_LIVE_TESTS=1` for the separately gated
SimpNetworkViz test. No model credentials are needed for these checks. An MCP
protocol test runs when the MCP extra is installed.

```bash
Rscript scripts/run_tests.R
R CMD build .
_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual PeakAnnoAgent_1.0.2.tar.gz
```

Installed-package tests use bundled fixtures and the package namespace. The
source runner also returns a failing status on failed expectations. If an
existing R installation has the `waldo::compare_proxy` namespace error, install
a current waldo in a dedicated test library; local verification used waldo 0.6.2
without changing the system R library.

Complete biological validation includes the 351-peak synthetic mm10 ATAC set,
its custom signatures and RNA-seq table, and the 1,000-peak ENCODE hg38 CTCF set.
It requires successful core steps, checks peak counts, and requires CTCF to be
the highest-ranked motif in the CTCF data:

```bash
python scripts/validate_biology.py --out release-artifacts/1.0.2/biology --rscript /path/to/Rscript
```

Or use the [official Bioconductor container](https://bioconductor.org/help/docker/):

```bash
docker build -t peakannoagent:1.0.2 .
docker run --rm --entrypoint /opt/peakanno-venv/bin/python \
  -v "$PWD/release-artifacts/1.0.2:/validation" peakannoagent:1.0.2 \
  /opt/PeakAnnoAgent/scripts/validate_biology.py --out /validation/biology
```

The container records installed R and Python packages. Its Bioconductor release
and Python constraints are fixed; R package downloads and base-image rebuilds
can still change. For exact long-term reproduction, retain the built image
digest and package lists alongside the release artifacts.

## External release checks

- Push reviewed changes and require the Python and R GitHub Actions jobs to pass.
- Run the manual biological workflow job to validate the Linux container.
- Validate a live Claude run with `ANTHROPIC_API_KEY` configured, using
  `scripts/validate_biology.py --live-llm`; this incurs API usage.
- Verify live ChIP-Atlas, literature/GEO and advertised optional backends separately.
- Attach the R tarball, Python wheel/source distribution, checksums and validation records to the reviewed release.

Run tests and analyses locally with the Python environment and R installation
selected in the commands above.
