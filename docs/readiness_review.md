# Repository readiness review

Reviewed: 2026-11-09. Scope: current working tree, including existing uncommitted changes.

This review includes local web, search, chat and demo additions that are not included
in this documentation/cleanup pull request. Findings and test counts describe that
reviewed working tree, rather than validation of this PR's base branch alone.

Assessment: substantial research prototype with implemented R tools, Python orchestration,
CLI, web interface, search, chat, reports, and demo outputs. It is not ready for a reliable
release or publication without addressing the correctness and packaging issues below.
An implementation-completion percentage would obscure these blockers.

## Findings, in priority order

1. **High: BED coordinates are not converted to R coordinates.**
   `R/utils.R:75` preserves the BED start, and `R/annotation.R:1` passes that start directly
   to `IRanges`. A BED interval `chr1 0 1` becomes an R range with start 0 and width 2,
   although the BED interval has width 1. This affects annotation and overlap boundaries.
   Define the coordinate convention for BED and TSV inputs explicitly, convert BED starts
   to one-based coordinates at the R boundary, preserve BED semantics for exports, and
   test one-base intervals and touching, non-overlapping intervals. Also review the
   inclusive `foverlaps` call in `R/qc.R:8` for BED blacklist boundary handling.

2. **High: failed core analysis can still return success.**
   `agent/peakannoagent/cli.py:cmd_run` returns 0 after `execute`, and
   `agent/peakannoagent/web.py:214` marks the job done whenever `execute` returns.
   A real offline run in the active R environment completed QC, failed annotation,
   GREAT and motifs because packages were missing, wrote an HTML report, and exited 0.
   Keep partial results, but distinguish successful, partial and failed analyses in
   the manifest, CLI exit status and web status. Define which steps are essential.

3. **High: web job identifiers collide.**
   `agent/peakannoagent/web.py:204` constructs identifiers from the second and a name slug.
   Two submissions with the same name in the same second use the same directory and
   replace the same registry entry. Reproduced by fixing the clock and creating two jobs.
   Use unique identifiers and exclusive directory creation before accepting uploads.

4. **Release blocker: installed-package R tests fail.**
   Source-level tests pass, but `R CMD check` fails. Tests resolve fixture files relative
   to the checkout layout, which is unavailable in the check test directory. The standalone
   runner also invokes tests outside the package namespace context, so internal helpers
   such as `enrichment_nodes` are unavailable. Make fixtures accessible in installed-package
   tests and use the standard testthat package runner; keep the standalone runner out of
   automatic package test discovery. The check environment additionally reports a
   `waldo`/`compare_proxy` error, which needs separate dependency-environment investigation.

5. **Release preparation remains incomplete.**
   The check reports 15 undocumented exports and undefined global bindings (`groups`,
   `i`, `j`, `node`). `DESCRIPTION` still contains a placeholder maintainer email.
   No CI workflow, dependency lockfile or container definition was found. Python bridge
   discovery assumes a source checkout; validate a built-wheel installation separately.
   Vignette warnings occurred with vignette building deliberately disabled and should
   be re-evaluated in a full release build.

## Validation collected

| Check | Result | Limits |
| --- | --- | --- |
| Python unittest discovery | 43 tests passed | Includes scripted LLM clients, R QC bridge, Flask test client and live SimpNetworkViz integration; does not validate live model calls or full biological backends |
| `Rscript tests/run_tests.R` | Passed | Active Rscript was R 4.1.3; repeated locale warnings |
| `R CMD build --no-build-vignettes --no-manual` | Passed | Build used framework R 4.3.2; no vignette rendering |
| `R CMD check --no-manual --no-vignettes`, force suggests disabled | 1 error, 3 warnings, 2 notes | Framework R 4.3.2; some suggested dependencies unavailable; package installs and loads |
| Real offline deterministic run | Report generated, core steps failed, exit 0 | Demonstrates failure-status issue; not a successful biological workflow |
| BED-to-GRanges boundary probe | BED width 1 became R width 2, start 0 | Reproduces coordinate conversion issue |
| Web job collision probe | Same identifier and directory | Two same-name jobs at a fixed timestamp |

The active Rscript had `data.table` and `jsonlite`, but lacked the main annotation,
pathway and motif packages. No new full end-to-end biological run or live LLM/API
validation was completed. Committed-looking demo artifacts in the working tree are
useful examples, but do not substitute for reproducible release validation; many
web/demo files are currently untracked.

## Claude file cleanup

No `CLAUDE.md`, `CLAUDE.local.md`, `.claude` directory or other Claude-named files were
found in this checkout. The initial review added and verified ignore rules for these
local instruction/settings files. At the user's subsequent request, this PR removes
the tracked `.gitignore` file entirely, superseding those ignore-rule additions.
Claude/Anthropic references in source and documentation implement the existing
application integration and were retained.
Existing user changes were preserved; no files were deleted during the review.

## PARCC execution requirement

The user specified that all test and analysis jobs must be submitted through `sbatch`.
The validation above ran before that instruction. No further test or analysis jobs
were started after it. Neither `sbatch` nor `squeue` was available in this tool session,
so these results must not be represented as PARCC batch-job validation. Future checks
must run through Slurm in the intended PARCC environment.

## Recommended completion order

1. Correct coordinate handling and add boundary regression tests.
2. Distinguish failed/partial analyses from successful runs and prevent web job collisions.
3. Repair installed-package tests, document exports and replace placeholder metadata.
4. Add reproducible dependencies and CI; validate wheel and R package installation.
5. Submit full mm10/hg38 example workflows through PARCC `sbatch`, recording job IDs,
   environment versions, outputs and essential-step success; validate live LLM mode separately.
