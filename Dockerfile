FROM docker.io/condaforge/miniforge3:latest

RUN mkdir /tmp/medi /tmp/medi/bin

COPY medi.yml Makefile patches/*.patch /tmp/medi

RUN mamba env create -n medi -f /tmp/medi/medi.yml && \
    . ${CONDA_PREFIX}/etc/profile.d/conda.sh && conda activate medi && \
    cd /tmp/medi && make report && mv /tmp/medi/bin/kraken2-report /bin && \
    patch ${CONDA_PREFIX}/envs/medi/share/kraken2-2.1.3-4/libexec/build_kraken2_db.sh /tmp/medi/build.patch && \
    patch ${CONDA_PREFIX}/envs/medi/share/kraken2-2.1.3-4/libexec/download_genomic_library.sh /tmp/medi/download_genomic.patch && \
    conda clean --tarballs --index-cache --packages --yes && \
    find ${CONDA_PREFIX} -follow -type f -name '*.a' -delete && \
    find ${CONDA_PREFIX} -follow -type f -name '*.pyc' -delete && \
    conda clean --force-pkgs-dirs --all --yes  && \
    rm -rf /tmp/medi

ENTRYPOINT ["mamba", "run", "-n", "medi"]
