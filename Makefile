CXX := g++
K2DIR := "src/kraken2"
SRC := ${K2DIR}/src

all: report patch

repo:
	rm -rf ${K2DIR}
	git clone https://github.com/daydream-boost/kraken2 ${K2DIR}

report: repo
	${CXX} -O3 -std=c++11 \
		${SRC}/mmap_file.cc ${SRC}/reports.cc ${SRC}/taxonomy.cc \
		${SRC}/kraken2-report.cpp -o ./bin/kraken2-report

patch:
	patch ${CONDA_DIR}/envs/medi/share/kraken2-2.1.3-4/libexec/build_kraken2_db.sh ./patches/build.patch && \
    patch ${CONDA_DIR}/envs/medi/share/kraken2-2.1.3-4/libexec/download_genomic_library.sh ./patches/download_genomic.patch

.PHONY: clean
clean:
	rm -rf ${K2DIR}