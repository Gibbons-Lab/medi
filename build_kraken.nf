#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.additional_dbs = ["bacteria", "archaea", "human", "viral", "plasmid", "UniVec_Core"]
params.max_db_size = Math.round(500e9)
params.confidence = 0.3

process setup_kraken_db {
    cpus 1
    memory "8 GB"

    output:
    path("kraken2_db")

    """
    kraken2-build --download-taxonomy --db kraken2_db
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
    cpus 20
    memory "380 GB"

    input:
    path(db)

    output:
    path("$db")

    """
    kraken2-build --build --db $db --threads ${task.cpus} --max-db-size ${params.max_db_size}
    """
}

process classify {
    cpus 8
    memory "32 GB"

    input:
    path(fi)
    each path(db)

    output:
    path("${fi}.k2")

    """
    kraken2 --db ${db} \
            --confidence ${params.confidence} \
            --threads ${task.cpus} --output ${fi}.k2 \
            --memory-mapping ${fi}
    """
}

process merge {
    cpus 1
    memory "64 GB"

    input:
    path(k2)
    path(db)

    output:
    path(db)

    """
    cat ${k2} > ${db}/database.kraken
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
    path("$db/library/*/*.fna")

    """
    ls ${db}/library/*/*.fna | wc -l
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
    //Channel.fromPath("${baseDir}/data/sequences/*.fna.gz").set{food_sequences}
    //setup_kraken_db()
    //add_existing(setup_kraken_db.out, params.additional_dbs)
    //add_sequences(food_sequences, add_existing.out.last())
    //build_kraken_db(add_sequences.out.last())
    Channel.fromPath("${launchDir}/work/7b/5bd3edd1f1e5f1118bb0e62ab1d620/kraken2_db").set{build_kraken_db}

    library(build_kraken_db)

    classify(library.out.flatten(), build_kraken_db)
    merge(classify.out.collect(), build_kraken_db)
    build_bracken(merge.out) | add_food_info
}
