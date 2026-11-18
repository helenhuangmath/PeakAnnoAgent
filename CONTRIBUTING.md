# Contributing

Report reproducible bugs through GitHub issues with the version, command, genome,
coordinate convention and relevant log. A small public or synthetic input is preferable.

Changes should preserve the JSON tool contract and BED coordinate convention. Add
regression coverage for correctness fixes. Tests should use scripted model clients
and local fixtures; network services and large biological backends belong in
explicit integration checks.

Run the Python suite (`python -m unittest discover -s tests`) before submitting a pull
request. If you change the R interface, also run `R CMD check` on `r/PeakAnnoAgent`.
For packaging changes, build the source distribution and wheel and run the
installed-wheel smoke check away from the source checkout.

Participation in this project is governed by the
[Code of Conduct](https://github.com/helenhuangmath/PeakAnnoAgent/blob/main/CODE_OF_CONDUCT.md).
