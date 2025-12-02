FROM conda-forge/miniforge3:latest

RUN mkdir /tmp/medi

COPY medi.yml Makefile patches /tmp/medi

RUN mamba env create -n medi -f /tmp/medi.yml && \
    patch ${CONDA_DIR}/envs/medi/share/kraken2-2.1.3-4/libexec/build_kraken2_db.sh /tmp/patches/build.patch && \
    patch ${CONDA_DIR}/envs/medi/share/kraken2-2.1.3-4/libexec/download_genomic_library.sh /tmp/patches/download_genomic.patch && \
    conda activate medi && make && \
    conda clean --tarballs --index-cache --packages --yes && \
    find ${CONDA_DIR} -follow -type f -name '*.a' -delete && \
    find ${CONDA_DIR} -follow -type f -name '*.pyc' -delete && \
    conda clean --force-pkgs-dirs --all --yes  && \
    echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate medi" >> /etc/skel/.bashrc && \
    echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate medi" >> ~/.bashrc && \
    rm -rf /tmp/medi
