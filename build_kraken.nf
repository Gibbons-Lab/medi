#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.additionalDbs = ["bacteria", "archaea", "human", "viral", "plasmid", "UniVec_Core"]
params.maxDbSize = "500 GB"
params.confidence = 0.3
params.threads = 20
params.rebuild = false
params.downloads = "${launchDir}/data"
params.out = "${launchDir}/data"
params.db = "${params.out}/medi_db"
params.additionalDecoys ="${params.out}/decoys"
params.useFtp = false
params.readLengths = [100, 150, 250]

/* Helper functions */

// Helper to calculate the required RAM for the Kraken2 database
def estimate_db_size(db, extra) {
    def db_size = null

    // Calculate db memory requirement
    db_size = MemoryUnit.of(db*.size().sum()) + extra
    log.info("Based on the hash size I am reserving ${db_size.toGiga()}GB of memory for Kraken2.")

    return db_size
}

workflow {
    println(params)

    rlens = channel.from(params.readLengths)

    if (!params.rebuild) {
        channel.fromPath("${params.out}/sequences/*.fna.gz").set{food_sequences}

        channel.fromPath("${params.additionalDecoys}/*.fna.gz")
        .set{decoy_sequences}

        sequences = decoy_sequences.concat(food_sequences)

        taxonomy = setup_kraken_db()
        add_existing(taxonomy, params.additionalDbs)
        add_sequences(sequences, taxonomy)
        lib = assemble_library(
            add_existing.out.collect(),
            add_sequences.out.collect()
        )
    } else {
        taxonomy = channel.of(file("${params.db}/taxonomy"))
        lib = channel.of(file("${params.out}/library"))
    }

    build_kraken_db(taxonomy, lib)
    self_classify(build_kraken_db.out, lib)
    build_bracken(build_kraken_db.out, self_classify.out, lib, rlens)

    add_info()
}


process setup_kraken_db {
    cpus 1
    memory "32 GB"
    time "12 h"

    output:
    path("medi_db/taxonomy")

    script:
    """
    kraken2-build --download-taxonomy --db medi_db
    """
}

process add_sequences {
    cpus 4
    memory "32 GB"
    time "48 h"

    input:
    path(fasta)
    path(taxonomy)

    output:
    path("medi_db/library/added/*.*")

    script:
    """
    mkdir medi_db && mv ${taxonomy} medi_db
    gunzip -c $fasta > ${fasta.baseName} && \
    kraken2-build --add-to-library ${fasta.baseName} --db medi_db --threads ${task.cpus} && \
    rm ${fasta.baseName}
    """
}

process add_existing {
    cpus 4
    memory "32 GB"
    time "48 h"

    input:
    path(taxonomy)
    each group

    output:
    path("medi_db/library/$group")

    script:
    if (group == "human")
        """
        mkdir medi_db && mv ${taxonomy} medi_db
        kraken2-build --download-library $group --db medi_db --no-mask --threads ${task.cpus} ${params.useFtp ? "--use-ftp" : ""}
        """
    else
        """
        mkdir medi_db && mv taxonomy medi_db
        kraken2-build --download-library $group --db medi_db --threads ${task.cpus} ${params.useFtp ? "--use-ftp" : ""}
        """
}

process assemble_library {
    cpus 1
    memory "8 GB"
    time "12h"
    publishDir params.out

    input:
    path(existing)
    path(sequences)

    output:
    path("library")

    script:
    """
    mkdir library && mkdir library/added
    mv ${existing} library
    mv ${sequences} library/added
    """
}

process build_kraken_db {
    cpus params.threads
    memory { MemoryUnit.of(params.maxDbSize) + 100.GB }
    cpus params.threads
    time "36 h"
    publishDir params.db

    input:
    path(taxonomy)
    path(library)

    output:
    tuple path("taxonomy"), path("*.k2d"), path("seqid2taxid.map")

    script:
    """
    alias find="find -L"  # Kraken2 bug see https://github.com/DerrickWood/kraken2/issues/67
    kraken2-build --build --db . \
        --threads ${task.cpus} \
        --max-db-size ${MemoryUnit.of(params.maxDbSize).toBytes()}
    """
}


process self_classify {
    cpus params.threads
    memory { estimate_db_size(bins, 250.GB) }
    time "1 d"

    input:
    tuple path(tax), path(bins), path(seqmap)
    path(library)

    output:
    path("database.kraken")

    script:
    """
    for f in ${library}/*/*.f*a; do
        echo "Classifying \$f..."
        kraken2 --db . --threads ${task.cpus} \
            --confidence ${params.confidence} \
            --threads ${task.cpus} \
            --memory-mapping \$f >> database.kraken
    done
    """
}

process build_bracken {
    cpus params.threads
    memory "64 GB"
    publishDir params.out
    time "12 h"
    publishDir params.db

    input:
    tuple path(taxonomy), path(bins), path(seqmap)
    path(self_class)
    path(lib)
    each rlen

    output:
    path("database*.kmer_distrib")

    script:
    """
    bracken-build -d . -t ${task.cpus} -k 35 -l ${rlen}
    """
}

process add_info {
    cpus 1
    memory "1 GB"
    publishDir params.db, mode: "copy", overwrite: true

    output:
    tuple path("manifest.csv"), path("food_matches.csv"), path("food_contents.csv.gz")

    script:
    """
    cp ${params.downloads}/dbs/{food_matches.csv,food_contents.csv.gz} .
    cp ${params.downloads}/manifest.csv .
    """
}
