#!/usr/bin/env nextflow

params.out_dir = "${baseDir}/data"
params.data_dir = "${baseDir}/data/raw"
params.kraken2_db = "${baseDir}/data/kraken2_db"
params.foods = "${baseDir}/data/dbs/food_matches.csv"
params.food_contents = "${baseDir}/data/dbs/food_contents.csv.gz"
params.single_end = true
params.trim_front = 5
params.min_length = 15
params.quality_threshold = 20
params.read_length = 150
params.threshold = 100

if (params.single_end) {
    Channel
        .fromPath("${params.data_dir}/*.fastq.gz")
        .map{row -> tuple(row.baseName.split("\\.fastq")[0], tuple(row))}
        .set{raw}
} else {
    Channel
        .fromFilePairs("${params.data_dir}/*_R{1,2}.fastq.gz")
        .ifEmpty { error "Cannot find any read files in ${params.data_dir}!" }
        .set{raw}
}

Channel
    .fromList(["D", "P", "G", "S"])
    .set{levels}

process preprocess {
    cpus 1
    publishDir "${params.out_dir}/preprocessed"
    input:
    set id, file(reads) from raw

    output:
    set id, file("${id}_filtered.fastq.gz"),
            file("${id}_fastp.json"),
            file("${id}.html") into processed_assembly, processed_align,
                                    processed_kraken

    script:
    if (params.single_end)
        """
        fastp -i ${reads[0]} -o ${id}_filtered.fastq.gz \
            --json ${id}_fastp.json --html ${id}.html \
            --trim_front1 ${params.trim_front} -l ${params.min_length} \
            -3 -M ${params.quality_threshold} -r
        """

    else
        """
        fastp -i ${reads[0]} -I ${reads[1]} \
            -o ${id}_filtered_R1.fastq.gz -O ${id}_filtered_R2.fastq.gz\
            --json ${id}_fastp.json --html ${id}.html \
            --trim_front1 ${params.trim_front} -l ${params.min_length} \
            -3 -M ${params.quality_threshold} -r
        """
}

process kraken {
    cpus 8
    publishDir "${params.out_dir}/kraken2"

    input:
    set id, file(reads), file(json), file(html) from processed_kraken

    output:
    set id, file("${id}.k2"), file("${id}.tsv") into kraken_reports

    script:
    if (params.single_end)
        """
        kraken2 --db ${params.kraken2_db} \
            --threads ${task.cpus} --gzip-compressed --output ${id}.k2 \
            --memory-mapping --report ${id}.tsv ${reads[0]}
        """

    else
        """
        kraken2 --db ${params.kraken2_db} --paired \
            --threads ${task.cpus} --gzip-compressed --output ${id}.k2 \
            --memory-mapping --report ${id}.tsv ${reads[0]} ${reads[1]}
        """
}

process count_taxa {
    cpus 4
    publishDir "${params.out_dir}/bracken"

    input:
    set id, file(kraken), file(report), lev from kraken_reports.combine(levels)

    output:
    set id, lev, file("${lev}_${id}.b2") into bracken_reports, reports_ready

    """
    mkdir ${lev} && cp ${report} ${lev}/${report} && \
        bracken -d ${params.kraken2_db} -i ${lev}/${report} \
        -l ${lev} -o ${lev}_${id}.b2 -r ${params.read_length} \
        -t ${params.threshold}
    """
}

process quantify {
    cpus 1

    input:
    set id, path(files) from bracken_reports
        .map{s -> tuple(s[0], s[2])}
        .groupTuple()

    output:
    path("${id}_food_content.csv") into food_contents
    path("${id}_food_abundance.csv") into food_abundances

    """
    Rscript $baseDir/scripts/quantify.R ${params.foods} ${params.food_contents} $files
    """
}


process merge_abundances {
    cpus 1
    publishDir "${params.out_dir}"

    input:
    path(files) from food_abundances.collect()

    output:
    path("food_abundances.csv") into merged_food_abundances

    """
    Rscript $baseDir/scripts/merge.R food_abundances.csv $files
    """
}

process merge_contents {
    cpus 1
    publishDir "${params.out_dir}"

    input:
    path(files) from food_contents.collect()

    output:
    path("food_contents.csv")

    """
    Rscript $baseDir/scripts/merge.R food_contents.csv $files
    """
}

process multiqc {
    publishDir "${params.out_dir}"

    input:
    path(foods) from merged_food_abundances

    output:
    path("multiqc_report.html")

    """
    multiqc ${params.out_dir}/preprocessed ${params.out_dir}/kraken2
    """
}
