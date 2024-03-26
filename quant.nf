#!/usr/bin/env nextflow
params.out_dir = "${launchDir}/data"
params.data_dir = "${launchDir}/data"
params.db = "data/medi_db"
params.foods = "${params.db}/food_matches.csv"
params.food_contents = "${params.db}/food_contents.csv.gz"
params.single_end = false
params.trim_front = 5
params.min_length = 50
params.quality_threshold = 20
params.read_length = 150
params.threshold = 10
params.confidence = 0.3
params.mapping = false

Channel
    .fromList(["D", "G", "S"])
    .set{levels}

process preprocess {
    cpus 4
    memory "8 GB"
    publishDir "${params.out_dir}/preprocessed"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id),
        path("${id}_filtered_R*.fastq.gz"),
        path("${id}_fastp.json"),
        path("${id}.html")

    script:
    if (params.single_end)
        """
        fastp -i ${reads[0]} -o ${id}_filtered_R1.fastq.gz \
            --json ${id}_fastp.json --html ${id}.html \
            --trim_front1 ${params.trim_front} -l ${params.min_length} \
            -3 -M ${params.quality_threshold} -r -w ${task.cpus}
        """

    else
        """
        fastp -i ${reads[0]} -I ${reads[1]} \
            -o ${id}_filtered_R1.fastq.gz -O ${id}_filtered_R2.fastq.gz\
            --json ${id}_fastp.json --html ${id}.html \
            --trim_front1 ${params.trim_front} -l ${params.min_length} \
            -3 -M ${params.quality_threshold} -r -w ${task.cpus}
        """
}

process kraken {
    cpus 8

    input:
    tuple val(id), path(reads), path(json), path(html)

    output:
    tuple val(id), path("${id}.k2"), path("${id}.tsv")

    script:
    if (params.single_end)
        """
        kraken2 --db ${params.db} \
            --confidence ${params.confidence} \
            --threads ${task.cpus} --gzip-compressed --output ${id}.k2 \
            --memory-mapping --report ${id}.tsv ${reads}
        """

    else
        """
        kraken2 --db ${params.db} --paired \
            --confidence ${params.confidence} \
            --threads ${task.cpus} --gzip-compressed --output ${id}.k2 \
            --memory-mapping --report ${id}.tsv  ${reads[0]} ${reads[1]}
        """
}

process architeuthis_filter {
    cpus 1
    publishDir "${params.out_dir}/kraken2", overwrite: true

    input:
    tuple val(id), path(k2), path(report)

    output:
    tuple val(id), path("${id}_filtered.k2"), path(report)

    """
    architeuthis mapping filter ${k2} \
        --data-dir ${params.db}/taxonomy \
        --min-consistency 0.95 --max-entropy 0.1 \
        --max-multiplicity 4 \
        --out ${id}_filtered.k2
    """


}

process summarize_mappings {
    cpus 1
    publishDir "${params.out_dir}/architeuthis"

    input:
    tuple val(id), path(k2), path(report)

    output:
    path("${id}_mapping.csv")

    """
    architeuthis mapping summary ${k2} --data-dir ${params.db}/taxonomy --out ${id}_mapping.csv
    """
}

process merge_mappings {
    cpus 1
    publishDir "${params.out_dir}", mode: "copy", overwrite: true

    input:
    path(mappings)

    output:
    path("mappings.csv")

    """
    architeuthis merge ${mappings} --out mappings.csv
    """
}

process count_taxa {
    cpus 4
    memory "16 GB"
    publishDir "${params.out_dir}/bracken", overwrite: true

    input:
    tuple val(id), path(kraken), path(report), val(lev)

    output:
    tuple val(id), val(lev), path("${lev}/${lev}_${id}.b2")

    """
    mkdir ${lev} && \
        sed 's/\\tR1\\t/\\tD\\t/g' ${report} > ${lev}/${report} && \
        bracken -d ${params.db} -i ${lev}/${report} \
        -l ${lev} -o ${lev}/${lev}_${id}.b2 -r ${params.read_length} \
        -t ${params.threshold} -w ${lev}/${id}_bracken.tsv
    """
}

process quantify {
    cpus 1
    memory "16 GB"
    publishDir "${params.out_dir}", mode: "copy", overwrite: true

    input:
    path(files)

    output:
    tuple path("food_abundance.csv"), path("food_content.csv")

    """
    quantify.R ${params.foods} ${params.food_contents} ${files}
    """
}

process merge_taxonomy {
    cpus 1
    memory "8 GB"

    input:
    tuple val(lev), path(reports)

    output:
    tuple val(lev), path("${lev}_merged.csv")

    """
    architeuthis merge ${reports} --out ${lev}_merged.csv
    """
}

process add_lineage {
    cpus 1
    memory "16 GB"
    publishDir "${params.out_dir}", mode: "copy", overwrite: true

    input:
    tuple val(lev), path(merged)

    output:
    path("${lev}_counts.csv")

    """
    architeuthis lineage ${merged} --data-dir ${params.db}/taxonomy --out ${lev}_counts.csv
    """
}


process multiqc {
    publishDir "${params.out_dir}", mode: "copy", overwrite: true

    input:
    path(report)

    output:
    path("multiqc_report.html")

    """
    multiqc ${params.out_dir}/preprocessed ${params.out_dir}/kraken2
    """
}

workflow {
    // find files
    if (params.single_end) {
        Channel
            .fromPath("${params.data_dir}/raw/*.fastq.gz")
            .map{row -> tuple(row.baseName.split("\\.fastq")[0], tuple(row))}
            .set{raw}
    } else {
        Channel
            .fromFilePairs([
                "${params.data_dir}/raw/*_R{1,2}_001.fastq.gz",
                "${params.data_dir}/raw/*_{1,2}.fastq.gz",
                "${params.data_dir}/raw/*_R{1,2}.fastq.gz"
            ])
            .ifEmpty { error "Cannot find any read files in ${params.data_dir}/raw!" }
            .set{raw}
    }

    // download taxa dbs
    //download_taxa_dbs()

    // quality filtering
    preprocess(raw)

    // quantify taxa abundances
    kraken(preprocess.out)
    architeuthis_filter(kraken.out)
    count_taxa(architeuthis_filter.out.combine(levels))
    count_taxa.out.map{s -> tuple(s[1], s[2])}
        .groupTuple()
        .set{merge_groups}
    merge_taxonomy(merge_groups)

    if (params.mapping) {
        // Get individual mappings
        summarize_mappings(architeuthis_filter.out)
        summarize_mappings.out.collect() | merge_mappings
    }

    // Add taxon lineages
    add_lineage(merge_taxonomy.out)

    // Quantify foods
    add_lineage.out.collect() | quantify

    // quality overview
    multiqc(merge_taxonomy.out.map{it[1]}.collect())
}
