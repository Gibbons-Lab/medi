#!/usr/bin/env nextflow
params.out_dir = "${launchDir}/data"
params.data_dir = "${launchDir}/data"
params.db = "${launchDir}/data/medi_db"
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
params.batchsize = 50
params.maxcpus = 24
params.dbmem = null

def db_size = null

Channel
    .fromList(["D", "G", "S"])
    .set{levels}

process preprocess {
    cpus 4
    memory "4 GB"
    publishDir "${params.out_dir}/preprocessed"
    time "1h"

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

process kraken_paired {
    cpus params.maxcpus
    memory db_size
    time 2.h + params.batchsize * 0.5.h

    input:
    tuple val(batch), val(ids), path(fwd_reads), path(rev_reads)

    output:
    path("*.k2")

    """
    #!/usr/bin/env python

    import sys
    import os
    from subprocess import run

    ids = "${ids.join(' ')}".split()
    fwd = "${fwd_reads}".split()
    rev = "${rev_reads}".split()

    assert len(ids) == len(fwd)
    assert len(ids) == len(rev)

    for i, idx in enumerate(ids):
        args = [
            "kraken2", "--db", "${params.db}", "--paired",
            "--confidence", "${params.confidence}",
            "--threads", "${task.cpus}", "--gzip-compressed",
            "--output", f"{idx}.k2", "--memory-mapping",
            fwd[i], rev[i]
        ]
        res = run(args)
        if res.returncode != 0:
            if os.path.exists(f"{idx}.k2"):
                os.remove(f"{idx}.k2")
            sys.exit(res.returncode)
    """
}

process kraken_single {
    cpus params.maxcpus
    memory db_size
    time 2.h + params.batchsize * 1.h

    input:
    tuple val(batch), val(ids), path(reads)

    output:
    path("*.k2")

    """
    #!/usr/bin/env python

    import sys
    import os
    from subprocess import run

    ids = "${ids.join(' ')}".split()
    fwd = "${reads}".split()

    assert len(ids) == len(fwd)

    for i, idx in enumerate(ids):
        args = [
            "kraken2", "--db", "${params.db}", "--paired",
            "--confidence", "${params.confidence}",
            "--threads", "${task.cpus}", "--gzip-compressed",
            "--output", f"{idx}.k2", "--memory-mapping", fwd[i]
        ]
        res = run(args)
        if res.returncode != 0:
            if os.path.exists(f"{idx}.k2"):
                os.remove(f"{idx}.k2")
            sys.exit(res.returncode)
    """
}

process architeuthis_filter {
    cpus 1
    publishDir "${params.out_dir}/kraken2", overwrite: true
    time 1.h
    memory "2 GB"

    input:
    tuple val(id), path(k2)

    output:
    tuple val(id), path("${id}_filtered.k2")

    """
    architeuthis mapping filter ${k2} \
        --data-dir ${params.db}/taxonomy \
        --min-consistency 0.95 --max-entropy 0.1 \
        --max-multiplicity 4 \
        --out ${id}_filtered.k2
    """
}

process kraken_report {
    cpus 1
    memory "200 MB"
    publishDir "${params.out_dir}/kraken2", overwrite: true
    time 30.m

    input:
    tuple val(id), path(k2)

    output:
    tuple val(id), path("*.tsv")

    """
    kraken2-report ${params.db}/taxo.k2d ${k2} ${id}.tsv
    """
}

process summarize_mappings {
    cpus 1
    publishDir "${params.out_dir}/architeuthis"
    time 1.h

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
    time 1.h

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
    memory "200 MB"
    publishDir "${params.out_dir}/bracken", overwrite: true
    time 1.h

    input:
    tuple val(id), path(report), val(lev)

    output:
    tuple val(id), val(lev), path("${lev}/${lev}_${id}.b2")

    """
    mkdir ${lev} && \
        fixk2report.R ${report} ${lev}/${report} && \
        bracken -d ${params.db} -i ${lev}/${report} \
        -l ${lev} -o ${lev}/${lev}_${id}.b2 -r ${params.read_length} \
        -t ${params.threshold} -w ${lev}/${id}_bracken.tsv
    """
}

process quantify {
    cpus 1
    memory "2 GB"
    publishDir "${params.out_dir}", mode: "copy", overwrite: true
    time 2.h

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
    memory "1 GB"
    time 2.h

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
    memory "4 GB"
    publishDir "${params.out_dir}", mode: "copy", overwrite: true
    time 2.h

    input:
    tuple val(lev), path(merged)

    output:
    path("${lev}_counts.csv")

    """
    architeuthis lineage ${merged} --data-dir ${params.db}/taxonomy --out ${lev}_counts.csv
    """
}


process multiqc {
    cpus 1
    memory "2 GB"
    publishDir "${params.out_dir}", mode: "copy", overwrite: true
    time 2.h

    input:
    path(report)

    output:
    path("multiqc_report.html")

    """
    multiqc ${params.out_dir}/preprocessed ${params.out_dir}/kraken2
    """
}

def batchify(ch, n, paired = true, batchsize = 10) {
    idx = Channel
        .from(0..(n-1))
        .map{it.intdiv(batchsize)}
    if (paired) {
        batched = idx.merge(ch)
            .map{tuple(it[0], it[1], it[2][0], it[2][1])}
            .groupTuple()
    } else {
        batched = idx.merge(ch).groupTuple()
    }
    return batched
}

workflow {
    // find files
    if (params.single_end) {
        Channel
            .fromPath("${params.data_dir}/raw/*.fastq.gz")
            .map{row -> tuple(row.baseName.split("\\.fastq")[0], tuple(row))}
            .set{raw}
        n = file("${params.data_dir}/raw/*.fastq.gz").size()
    } else {
        Channel
            .fromFilePairs([
                "${params.data_dir}/raw/*_R{1,2}_001.fastq.gz",
                "${params.data_dir}/raw/*_{1,2}.fastq.gz",
                "${params.data_dir}/raw/*_R{1,2}.fastq.gz"
            ])
            .ifEmpty { error "Cannot find any read files in ${params.data_dir}/raw!" }
            .set{raw}
        n = file("${params.data_dir}/raw/*.f*.gz").size() / 2
    }

    // Calculate db memory requirement
    if (params.dbmem) {
        db_size = MemoryUnit.of("${params.dbmem} GB")
    } else {
        db_size = MemoryUnit.of(file("${params.db}/hash.k2d").size()) + 6.GB
        log.info("Based on the hash size I am reserving ${db_size.toGiga()}GB of memory for Kraken2.")
    }

    // quality filtering
    preprocess(raw)


    // quantify taxa abundances
    batched = batchify(preprocess.out, n, !params.single_end, params.batchsize)
    if (params.single_end) {
        k2 = kraken_single(batched)
    } else {
        k2 = kraken_paired(batched)
    }

    k2.flatten().map{tuple it.baseName.split(".k2")[0], it} | architeuthis_filter | kraken_report
    count_taxa(kraken_report.out.combine(levels))
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
