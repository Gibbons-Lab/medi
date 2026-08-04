#!/usr/bin/env bash

CONDA_DIR=$1
ENV_NAME="${2:-medi}"

echo "🔨 Building the report generator..." && \
. ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate ${ENV_NAME} && \
make report && \
echo "🔨 Patching kraken2..." && \
patch ${CONDA_DIR}/envs/${ENV_NAME}/medi/share/kraken2-2.1.3-4/libexec/build_kraken2_db.sh ./patches/build.patch && \
patch ${CONDA_DIR}/envs/${ENV_NAME}/share/kraken2-2.1.3-4/libexec/download_genomic_library.sh ./patches/download_genomic.patch