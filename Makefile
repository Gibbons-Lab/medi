K2DIR := "src/kraken2"
SRC := ${K2DIR}/src

repo:
	rm -rf ${K2DIR}
	git clone https://github.com/daydream-boost/kraken2 ${K2DIR}

report: repo
	g++ -O3 -std=c++11 \
		${SRC}/mmap_file.cc ${SRC}/reports.cc ${SRC}/taxonomy.cc \
		${SRC}/kraken2-report.cpp -o ./bin/kraken2-report

clean:
	rm -rf ${K2DIR}