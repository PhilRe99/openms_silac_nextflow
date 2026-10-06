FROM openms-silac/sdrf-pipelines:dea8ba9

# Nextflow report/trace/timeline collection runs ps inside the task container.
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
