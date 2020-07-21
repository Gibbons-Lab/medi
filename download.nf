params.data_dir = "$baseDir/data"

Channel
    .fromPath("${params.data_dir}/PRJEB29065.tsv")
    .splitCsv(header: true, sep: "\t")
    .map{row -> tuple(row.submitted_ftp.split("/").last().replaceAll(".fastq.gz", ""), row.submitted_ftp)}
    .set{daily_files}

Channel
    .fromPath("${params.data_dir}/ibdmdb_healthy.csv")
    .splitCsv(header: true)
    .map{row -> tuple(row["External ID"], "https://ibdmdb.org/tunnel/static/HMP2/WGS/1818/${row['External ID']}.tar")}
    .set{ibdmdb_files}

process download_daily {
    cpus 1
    publishDir "${params.data_dir}/ssgx/raw/"

    input:
    tuple id, url from daily_files

    output:
    tuple id, path("${id}.fastq.gz") into daily_fastq

    """
    wget $url -O ${id}.fastq.gz
    """
}

process download_ibd {
    cpus 1
    publishDir "${params.data_dir}/mgx/raw"

    input:
    tuple id, url from ibdmdb_files

    output:
    tuple id, path("${id}_R1.fastq.gz"), path("${id}_R2.fastq.gz") into ibdmdb_fastq

    """
    wget $url -O ${id}.tar && \
    tar -xf ${id}.tar && rm ${id}.tar
    """
}

