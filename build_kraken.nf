#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.additionalDbs = ["bacteria", "archaea", "human", "viral", "plasmid", "UniVec_Core"]
params.maxDbSize = 500
params.confidence = 0.3
params.threads = 20
params.rebuild = false
params.downloads = "${launchDir}/data"
params.out = "${launchDir}/data"
params.db = "${params.out}/medi_db"


workflow {
    if (!params.rebuild) {
        Channel.fromPath("${params.downloads}/sequences/*.fna.gz").set{food_sequences}
        setup_kraken_db()
        add_existing(setup_kraken_db.out, params.additionalDbs)
        add_sequences(food_sequences, add_existing.out.last())
        db = add_sequences.out.last()
    } else {
        db = Channel.fromPath(params.db)
    }

    build_kraken_db(db)
    self_classify(build_kraken_db.out)
    build_bracken(self_classify.out) | add_info
}


process setup_kraken_db {
    cpus 1
    memory "8 GB"

    output:
    path("medi_db")

    script:
    """
    kraken2-build --download-taxonomy --db medi_db
    """
}

process add_sequences {
    cpus 5
    memory "16 GB"
    publishDir params.out

    input:
    path(fasta)
    path(db)

    output:
    path("$db")

    script:
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
    cpus params.threads
    memory "${params.maxDbSize} GB"

    input:
    path(db)

    output:
    path("$db")

    script:
    """
    kraken2-build --build --db $db \
        --threads ${task.cpus} \
        --max-db-size ${(params.max_db_size as BigInteger) * (1000G**3)}
    """
}


process self_classify {
    cpus 1
    memory "${params.maxDbSize} GB"

    input:
    path(db)

    output:
    path(db)

    script:
    """
    kraken2 --db ${db} --threads ${task.cpus} \
        --confidence ${params.confidence} \
        --threads ${task.cpus} \
        --memory-mapping ${db}/library/*/*.f*a > ${db}/database.kraken
    """
}

process build_bracken {
    cpus 20
    memory "64 GB"
    publishDir params.out

    input:
    path(db)

    output:
    path("$db")

    script:
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

    script:
    """
    ls ${db}/library/*/*.f*a | wc -l
    """
}

process add_info {
    cpus 1
    memory "1 GB"
    publishDir params.out

    input:
    path(db)

    output:
    path("$db")

    script:
    """
    cp ${params.downloads}/dbs/{food_matches.csv,food_contents.csv.gz} ${db}
    cp ${params.downloads}/manifest.csv ${db}
    """
}
