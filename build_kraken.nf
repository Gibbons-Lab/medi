#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.additional_dbs = ["bacteria", "archaea", "human", "viral", "plasmid", "UniVec_Core"]
params.max_db_size = 500
params.confidence = 0.3
params.max_threads = 20

process setup_kraken_db {
    cpus 1
    memory "8 GB"

    output:
    path("medi_db")

    """
    kraken2-build --download-taxonomy --db medi_db
    """
}

process add_sequences {
    cpus 5
    memory "16 GB"

    input:
    path(fasta)
    path(db)

    output:
    path("$db")

    """
    gunzip -c $fasta > ${fasta.baseName} && \
    kraken2-build --add-to-library ${fasta.baseName} --db $db --threads ${task.cpus} && \
    rm ${fasta.baseName}
    """
}

process add_existing {
    cpus 4
    memory "16 GB"

    input:
    path(db)
    each group

    output:
    path("$db")

    script:
    if (group == "human")
        """
        kraken2-build --download-library $group --db $db --no-mask --threads ${task.cpus}
        """
    else
        """
        kraken2-build --download-library $group --db $db --threads ${task.cpus}
        """
}

process build_kraken_db {
    cpus params.max_threads
    memory "${params.max_db_size} GB"

    input:
    path(db)

    output:
    path("$db")

    """
    kraken2-build --build --db $db \
        --threads ${task.cpus} \
        --max-db-size ${(params.max_db_size as BigInteger) * (1000G**3)}
    """
}


process self_classify {
    cpus 1
    memory "${params.max_db_size} GB"

    input:
    path(k2)
    path(db)

    output:
    path(db)

    """
    kraken2 --db ${db} --threads ${task.cpus} \
        --confidence ${params.confidence} \
        --threads ${task.cpus} \
        --memory-mapping <( cat ${k2} ) > ${db}/database.kraken
    """
}

process build_bracken {
    cpus 20
    memory "64 GB"
    publishDir "$baseDir/data"

    input:
    path(db)

    output:
    path("$db")

    """
    bracken-build -d $db -t ${task.cpus} -k 35 -l 100 && \
    bracken-build -d $db -t ${task.cpus} -k 35 -l 150
    """
}

process library {
    cpus 1
    memory "4 GB"

    input:
    path(db)

    output:
    path("$db/library/*/*.f*a")

    """
    ls ${db}/library/*/*.f*a | wc -l
    """
}

process add_food_info {
    cpus 1
    memory "1 GB"
    publishDir "$baseDir/data"

    input:
    path(db)

    output:
    path("$db")

    """
    cp ${launchDir}/data/{food_matches.csv,food_contents.csv.gz} ${db}
    """
}

workflow {
    Channel.fromPath("${baseDir}/data/sequences/*.fna.gz").set{food_sequences}
    setup_kraken_db()
    add_existing(setup_kraken_db.out, params.additional_dbs)
    add_sequences(food_sequences, add_existing.out.last())
    build_kraken_db(add_sequences.out.last())

    library(build_kraken_db.out)

    self_classify(library.out, build_kraken_db.out)
    build_bracken(self_classify.out) | add_food_info
}
