# PeakAnnoAgent v1.0.2

Version: **1.0.2**. Tag: **v1.0.2**. Release title: **v1.0.2**.

Released November 12, 2026.

Pathway networks now lay out modules separately, use consistent module colours and readable labels,
and export node positions. The release also corrects BED coordinate conversion and blacklist
boundary overlaps, records success/partial/failed outcomes consistently, isolates web jobs, and
bundles the R engine with Python distributions. See [the changelog](../CHANGELOG.md).

BED and Chr/Start/End TSV inputs use zero-based, half-open coordinates. Re-run pre-1.0.2 analyses
before using quantitative annotation or overlap results.

Local validation completed before publication: Python source and installed-wheel suites passed
(62 tests, one optional integration skipped); R package checks had zero errors and zero warnings,
with one note for unavailable optional suggested packages. Complete mm10 ATAC and hg38 CTCF
workflows succeeded, and CTCF was the top recovered motif. Public Europe PMC, NCBI and ChIP-Atlas
endpoints were reachable with a configured system CA bundle. Testing files and local validation
outputs are not included in this release.

Live model behavior and Docker execution have not been verified in this environment.
Python 3.10+ and the documented R dependencies are required. Use `peakanno-agent check --offline`
to check the selected local installation. Resolved Python dependencies are recorded in
`constraints-release.txt`; the container recipe uses Bioconductor 3.20. The web app supports local,
single-user operation.
