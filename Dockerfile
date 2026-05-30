FROM mambaorg/micromamba:1.5.10

WORKDIR /workspace

# System dependencies
USER root
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        nodejs \
        npm \
        wget \
        tmux \
        git \
        libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

# Copy environment file
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# Create conda env with micromamba
RUN micromamba env create -f /tmp/environment.yml -y && \
    micromamba clean --all --yes

# Make bcai the default env
ENV MAMBA_DEFAULT_ENV=bcai
ENV CONDA_DEFAULT_ENV=bcai
ENV PATH=/opt/conda/envs/bcai/bin:$PATH

# Auto-activate bcai in every bash (also inside tmux panes)
RUN echo 'eval "$(micromamba shell hook -s bash)"' >> ~/.bashrc && \
    echo 'micromamba activate bcai' >> ~/.bashrc

# Use bash -lc for subsequent RUN instructions
SHELL ["bash", "-lc"]

# Optional sanity check (use micromamba run to exec inside bcai)
RUN echo "Conda envs:" && micromamba env list && \
    micromamba run -n bcai python -c "import sys; print('Python', sys.version)" && \
    micromamba run -n bcai python -c "import torch; print('Torch', __import__('torch').__version__)"

EXPOSE 8501

# Add entrypoint
COPY entrypoint.sh /workspace/entrypoint.sh
ENTRYPOINT ["/workspace/entrypoint.sh"]
	
# Default command
CMD ["bash"]
