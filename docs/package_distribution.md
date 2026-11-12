# Installing and publishing packages

PeakAnnoAgent already has two installable packages: the R analysis engine and
the Python agent. The v1.0.2 GitHub release contains their built distributions.

## Install v1.0.2 now

Install the Python wheel directly from the public release:

```bash
python -m pip install "peakannoagent[web,agent] @ https://github.com/helenhuangmath/PeakAnnoAgent/releases/download/v1.0.2/peakannoagent-1.0.2-py3-none-any.whl"
```

Install the R engine in R:

```r
install.packages("https://github.com/helenhuangmath/PeakAnnoAgent/releases/download/v1.0.2/PeakAnnoAgent_1.0.2.tar.gz",
                 repos = NULL, type = "source")
```

Install the required R/Bioconductor dependencies first; installing the Python
package does not install R, Bioconductor or genome databases.

## Build the packages locally

```bash
python -m pip install -c constraints-release.txt build
python scripts/build_packages.py --out dist
```

Set `--r /path/to/R` for a specific R installation. Pass `--no-isolation` when
using an already prepared Python build environment. The build stages only R
package files and excludes testing files from the distributions. It does not
depend on any repository dotfiles. Output is an R source package, a Python
wheel/source package, and SHA-256 checksums.

## Publish the Python package to PyPI

1. Create a PyPI account with two-factor authentication and confirm that the
   `peakannoagent` project name is available or controlled by the maintainer.
2. Build the distributions and inspect their metadata with `python -m twine check`.
3. Configure a project-scoped PyPI API token in the local publishing environment.
4. Upload only the Python distributions:

```bash
python -m pip install twine
python -m twine check dist/peakannoagent-1.0.2.tar.gz dist/peakannoagent-1.0.2-py3-none-any.whl
python -m twine upload dist/peakannoagent-1.0.2.tar.gz dist/peakannoagent-1.0.2-py3-none-any.whl
```

After publication users can install `pip install "peakannoagent[web,agent]"`.
Do not upload the R tarball to PyPI. Published filenames cannot be reused for
different package contents; use a new version for later fixes. Trusted
Publishing is an alternative if hosted publication workflows are added later.

See the [Python packaging recommendations](https://packaging.python.org/en/latest/guides/tool-recommendations/).

## Publish the R package to CRAN

The GitHub R tarball is installable now. For `install.packages("PeakAnnoAgent")`
through CRAN, first complete platform checks and resolve warnings or significant
notes from `R CMD check --as-cran`. Review optional dependency availability,
function help, examples, and external-service behavior, then submit the source
tarball through CRAN's submission form. Repository removal of dotfiles does not
prevent submission: build the clean R tarball with the staging script above.

See [CRAN's submission checklist](https://cran.r-project.org/web/packages/submission_checklist.html)
and [repository policy](https://cran.r-project.org/web/packages/policies.html).

PyPI publication requires maintainer credentials; CRAN additionally reviews
submissions. A GitHub release does not register a package with either service.
