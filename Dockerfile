FROM python:3.12-slim

WORKDIR /opt/PeakAnnoAgent
COPY . .
RUN pip install --no-cache-dir -c constraints-release.txt ".[web,agent,mcp]" \
    && pip freeze > /opt/PeakAnnoAgent/python-packages.txt
# Gene models and gene-set libraries are downloaded on first use and kept here; mount a volume
# on /cache to keep them between runs.
ENV PEAKANNOAGENT_CACHE=/cache MPLCONFIGDIR=/tmp/matplotlib
RUN mkdir -p /cache
WORKDIR /work
EXPOSE 8766
ENTRYPOINT ["peakanno-agent"]
CMD ["--help"]
