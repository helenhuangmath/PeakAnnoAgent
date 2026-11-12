FROM bioconductor/bioconductor_docker:RELEASE_3_20

RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-venv python3-pip && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/PeakAnnoAgent
COPY DESCRIPTION NAMESPACE LICENSE LICENSE.md README.md constraints-release.txt ./
COPY R/ R/
COPY man/ man/
COPY inst/ inst/
COPY vignettes/ vignettes/
COPY examples/data/ examples/data/
COPY scripts/install_dependencies.R scripts/peak_agent.R scripts/peak_tools.R scripts/
COPY agent/pyproject.toml agent/setup.py agent/MANIFEST.in agent/README.md agent/
COPY agent/peakannoagent/ agent/peakannoagent/
RUN Rscript scripts/install_dependencies.R && R CMD INSTALL .
RUN python3 -m venv /opt/peakanno-venv && /opt/peakanno-venv/bin/pip install -c constraints-release.txt "./agent[web,agent,mcp]"
RUN Rscript -e 'writeLines(capture.output(sessionInfo()), "/opt/PeakAnnoAgent/R-session.txt"); write.table(installed.packages()[, c("Package", "Version")], "/opt/PeakAnnoAgent/R-packages.tsv", sep="\t", row.names=FALSE, quote=FALSE)'
RUN /opt/peakanno-venv/bin/pip freeze > /opt/PeakAnnoAgent/python-packages.txt
ENV PATH="/opt/peakanno-venv/bin:${PATH}"
WORKDIR /work
EXPOSE 8766
ENTRYPOINT ["peakanno-agent"]
CMD ["--help"]
