FROM quay.io/biocontainers/sdrf-pipelines@sha256:b259135d4bf03dd76461eec481b755f745acb51c1f9d23a4709fbcfa0c999977

ENV XDG_CACHE_HOME=/opt/sdrf-cache \
    PYTHONDONTWRITEBYTECODE=1

RUN python -c \
    "from sdrf_pipelines.ols.ols import download_ontology_cache; download_ontology_cache()" \
    && chmod -R a+rX /opt/sdrf-cache