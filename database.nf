#!/usr/bin/env nextflow

// Build the FOODB genome database

params.max_threads = 20

nextflow.enable.dsl = 2

foodb = "https://foodb.ca/public/system/downloads/foodb_2020_4_7_csv.tar.gz"
genbank_summary = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
taxdump = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"


process download {
    cpus 1
    publishDir "$baseDir/data/dbs"

    output:
    tuple path("foodb"), path("genbank_summary.tsv")

    """
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 $foodb -O foodb.tgz && \
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 $genbank_summary -O genbank_summary.tsv && \
    tar -xf foodb.tgz && \
    mv foodb_2020_04_07_csv foodb
    """
}

process get_taxids {
    cpus 1

    input:
    tuple path(foodb), path(gb_summary)

    output:
    tuple path("foodb"), path("taxids.tsv"), path("$gb_summary")

    """
    #!/usr/bin/env Rscript

    library(data.table)

    dt <- fread("$genbank_summary", sep="\t")[
        grepl("ftp.ncbi.nlm.nih.gov", ftp_path, fixed = TRUE)
    ]
    genbank <- dt[!is.na(taxid), .(taxid = as.character(unique(taxid)))]
    genbank[, "source" := "genbank"]
    dt <- fread("$foodb/Food.csv")
    foodb <- dt[!is.na(ncbi_taxonomy_id), .(taxid = ncbi_taxonomy_id)]
    foodb[, "source" := "foodb"]
    fwrite(rbind(genbank, foodb), "taxids.tsv", col.names=F, sep="\t")
    """
}

process download_taxa_dbs {
    cpus 1

    output:
    path("taxdump")

    """
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 \
        ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz && \
        mkdir taxdump && tar -xf taxdump.tar.gz --directory taxdump
    """
}

process get_lineage {
    cpus 1

    input:
    tuple path(foodb), path(taxids), path(gb_summary), path(taxadb)

    output:
    tuple path("$foodb"), path("lineage.txt"), path("lineage_ids.txt"), path("${gb_summary}")

    """
    taxonkit lineage --data-dir $taxadb -i 1 $taxids > raw.txt && \
    taxonkit reformat --data-dir $taxadb -i 3 raw.txt > lineage.txt && \
    taxonkit reformat --data-dir $taxadb -t -i 3 raw.txt > lineage_ids.txt
    """
}

process match_taxids {
    cpus 1
    publishDir "$baseDir/data"

    input:
    tuple path(foodb), path(lineage), path(lineage_ids), path(gb_summary)

    output:
    path("matches.csv")

    """
    match.R $lineage_ids $gb_summary
    """
}

process download_sequences {
    cpus 8
    memory "64 GB"

    publishDir "$baseDir/data/"

    input:
    path(matches)

    output:
    tuple path("sequences/*.fna.gz"), path("manifest.csv")

    """
    download.R $matches $task.cpus sequences
    """
}


process food_mappings {
    cpus 1
    memory "64 GB"
    publishDir "$baseDir/data/dbs"

    input:
    path(matches)

    output:
    tuple path("food_matches.csv"), path("food_contents.csv.gz")

    """
    food_mapping.R $baseDir/data/dbs/foodb $matches
    """
}

process sketch {
    cpus 2
    memory "4 GB"
    publishDir "${launchDir}/data/sketches"

    input:
    path(seq)

    output:
    path("*.sig")

    """
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000 ${seq}
    """
}

process ANI {
    cpus params.max_threads
    memory "64 GB"
    publishDir "${launchDir}/data", mode: "copy", overwite: true

    input:
    path(sigs)

    output:
    path("mash_ani.csv")

    """
    sourmash compare -k 21 --ani -p ${task.cpus} --csv mash_ani.csv ${sigs}
    """
}


workflow {
    download() | get_taxids
    download_taxa_dbs()
    get_lineage(get_taxids.out.combine(download_taxa_dbs.out))
        | match_taxids
        | download_sequences

    download_sequences.out.map{it[0]}.flatten().set{seqs}


    seqs | sketch
    ANI(sketch.out.collect())

    food_mappings(match_taxids.out)
}
