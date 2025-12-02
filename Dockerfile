FROM condaforge/miniforge3:latest

RUN mkdir /tmp/medi /tmp/medi/bin

COPY medi.yml Makefile patches /tmp/medi

RUN mamba env create -n medi -f /tmp/medi/medi.yml && \
    . ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate medi && \
    cd /tmp/medi && make report && mv /tmp/medi/bin/kraken2-report /bin && \
    patch ${CONDA_DIR}/envs/medi/share/kraken2-2.1.3-4/libexec/build_kraken2_db.sh /tmp/medi/patches/build.patch && \
    patch ${CONDA_DIR}/envs/medi/share/kraken2-2.1.3-4/libexec/download_genomic_library.sh /tmp/medi/patches/download_genomic.patch && \
    conda clean --tarballs --index-cache --packages --yes && \
    find ${CONDA_DIR} -follow -type f -name '*.a' -delete && \
    find ${CONDA_DIR} -follow -type f -name '*.pyc' -delete && \
    conda clean --force-pkgs-dirs --all --yes  && \
    echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate medi" >> /etc/skel/.bashrc && \
    echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate medi" >> ~/.bashrc && \
    rm -rf /tmp/medi
