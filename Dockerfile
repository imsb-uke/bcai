FROM mambaorg/micromamba:1.5.10

WORKDIR /workspace

USER root
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        nodejs \
        npm \
        wget \
        git \
        libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# Create the conda env
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba env create -f /tmp/environment.yml -y && micromamba clean --all --yes

ENV MAMBA_DEFAULT_ENV=bcai
ENV CONDA_DEFAULT_ENV=bcai
ENV PATH=/opt/conda/envs/bcai/bin:$PATH
# Unbuffer stdout/stderr so print()-based logging shows up in `docker logs`
# immediately instead of sitting in a pipe buffer.
ENV PYTHONUNBUFFERED=1

# Add FastAPI + auth deps
COPY --chown=$MAMBA_USER:$MAMBA_USER api/requirements.txt /tmp/api-requirements.txt
RUN pip install --no-cache-dir -r /tmp/api-requirements.txt

# Copy application code
COPY --chown=$MAMBA_USER:$MAMBA_USER api/           /workspace/api/
COPY --chown=$MAMBA_USER:$MAMBA_USER mcp/           /workspace/mcp/
COPY --chown=$MAMBA_USER:$MAMBA_USER config/        /workspace/config/
COPY --chown=$MAMBA_USER:$MAMBA_USER docs/          /workspace/docs/
COPY --chown=$MAMBA_USER:$MAMBA_USER case_studies/  /workspace/case_studies/

EXPOSE 8000

# data/ is for drug libraries; user_data/ holds per-user dirs + SQLite DB.
# Both are volume-mounted at runtime; mkdir ensures they exist before the mount.
RUN mkdir -p /workspace/data /workspace/user_data && \
    chown $MAMBA_USER:$MAMBA_USER /workspace/data /workspace/user_data

COPY --chown=$MAMBA_USER:$MAMBA_USER entrypoint.sh /workspace/entrypoint.sh
ENTRYPOINT ["/workspace/entrypoint.sh"]
