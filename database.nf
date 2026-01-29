#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.threads = 20
params.out = "${launchDir}/data"
params.additionalDecoys = null //"${params.out}/decoys.csv"
params.useRsync = true

workflow {
    def foodb = "https://foodb.ca/public/system/downloads/foodb_2020_4_7_csv.tar.gz"
    def genbank_summary = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
    def taxdump = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

    download_foodb_genbank(foodb, genbank_summary)
    curate_taxids(download_foodb_genbank.out, file("${projectDir}/data/missing_foodb_curated.csv"))
    get_taxids(download_foodb_genbank.out, curate_taxids.out)

    download_taxa_dbs(taxdump)
    get_lineage(get_taxids.out.combine(download_taxa_dbs.out))
        | match_taxids

    nuc = download_nucleotide(match_taxids.out)
    match_taxids.out
        .splitCsv(header: true)
        .filter{row -> row.db == "genbank"}
        .map{row -> row.id}
        .unique()
        .set{gb_ids}
    download_genbank(gb_ids.combine(match_taxids.out))
    download_genbank.out.collect() | merge_genbank
    merge_all(nuc.combine(merge_genbank.out))

    merge_all.out.map{it[1]}.flatten().set{seqs}

    seqs | sketch
    ANI(sketch.out.collect())

    println(params)
    // Add more decoy sequences
    if (params.additionalDecoys) {
        channel.fromPath("${params.additionalDecoys}").set{decoy_manifest}
        download_decoys(decoy_manifest)
    }

    food_mappings(download_foodb_genbank.out, match_taxids.out, curate_taxids.out)
}



process download_foodb_genbank {
    cpus 1
    memory "2 GB"
    publishDir "${params.out}/dbs"

    input:
    val foodb
    val genbank_summary

    output:
    tuple path("foodb"), path("genbank_summary.tsv")

    script:
    """
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 ${foodb} -O foodb.tgz && \
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 ${genbank_summary} -O genbank_summary.tsv && \
    tar -xf foodb.tgz && \
    mv foodb_*_csv foodb
    """
}

process curate_taxids {
    cpus 1
    memory "4 GB"
    time "30 m"
    publishDir params.out

    input:
    tuple path(foodb), path(gb_summary)
    path(curated)

    output:
    path("curated_foods.csv")

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    foods <- fread("foodb/Food.csv")
    foods[, "revised_taxonomy_id" := ncbi_taxonomy_id]
    setkey(foods, id)
    curated <- fread("${curated}")
    foods[J(curated[["id"]]), revised_taxonomy_id := curated[["revised_taxonomy_id"]]]

    fwrite(foods, "curated_foods.csv")
    """
}

process get_taxids {
    cpus 1
    memory "4 GB"
    time "1 h"

    input:
    tuple path(foodb), path(gb_summary)
    path(curated)

    output:
    tuple path("foodb"), path("taxids.tsv"), path("${gb_summary}")

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    dt <- fread("${gb_summary}", sep="\t")[
        grepl("ftp.ncbi.nlm.nih.gov", ftp_path, fixed = TRUE)
    ]
    genbank <- dt[!is.na(taxid), .(taxid = as.character(unique(taxid)))]
    genbank[, "source" := "genbank"]
    dt <- fread("${curated}")
    foodb <- dt[!is.na(ncbi_taxonomy_id), .(taxid = ncbi_taxonomy_id)]
    foodb[, "source" := "foodb"]
    fwrite(rbind(genbank, foodb), "taxids.tsv", col.names=F, sep="\\t")
    """
}

process download_taxa_dbs {
    cpus 1
    memory "500 MB"
    time "2 h"

    input:
    val(taxdump)

    output:
    path("taxdump")

    script:
    """
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 \
        ${taxdump} && \
        mkdir taxdump && tar -xf taxdump.tar.gz --directory taxdump
    """
}

process get_lineage {
    cpus 1
    memory "8 GB"
    time "1 h"

    input:
    tuple path(foodb), path(taxids), path(gb_summary), path(taxadb)

    output:
    tuple path("$foodb"), path("lineage.txt"), path("lineage_ids.txt"), path("${gb_summary}")

    script:
    """
    taxonkit lineage --data-dir $taxadb -i 1 $taxids > raw.txt && \
    taxonkit reformat -f "{K};{p};{c};{o};{f};{g};{s}" --data-dir $taxadb -i 3 raw.txt > lineage.txt && \
    taxonkit reformat -f "{K};{p};{c};{o};{f};{g};{s}" --data-dir $taxadb -t -i 3 raw.txt > lineage_ids.txt
    """
}

process match_taxids {
    cpus 1
    publishDir params.out, mode: "copy", overwrite: true
    memory "8 GB"
    time "1 h"

    input:
    tuple path(foodb), path(lineage), path(lineage_ids), path(gb_summary)

    output:
    path("matches.csv")

    script:
    """
    match.R $lineage_ids $gb_summary
    """
}

process download_nucleotide {
    cpus 2
    memory "32 GB"
    time "48h"

    input:
    path(matches)

    output:
    tuple path("nucleotide.csv"), path("sequences/*.fna.gz")

    script:
    """
    download.R $matches "nucleotide" sequences "all" ${params.useRsync}
    """
}

process download_genbank {
    cpus 2
    memory "32 GB"
    time "8h"
    errorStrategy "ignore"

    input:
    tuple val(id), path(matches)

    output:
    tuple path("sequences/*.fna.gz"), path("${id}.csv")

    script:
    """
    download.R $matches "genbank" sequences "${id}" ${params.useRsync}
    """
}

process download_decoys {
    cpus 2
    memory "32 GB"
    time "48h"
    publishDir params.out

    input:
    path(decoys)

    output:
    tuple path("decoys.csv"), path("decoys/*.fna.gz")

    script:
    """
    download.R $decoys decoys decoys "all" ${params.useRsync}
    """
}

process merge_genbank {
    cpus 1
    memory "4 GB"
    publishDir params.out
    time "1 h"

    input:
    path(files)

    output:
    tuple path("genbank.csv"), path("sequences/*.fna.gz")

    script:
    """
    mkdir sequences && mv *.fna.gz sequences/
    merge.R genbank.csv *.csv
    """
}

process merge_all {
    cpus 1
    memory "4 GB"
    publishDir params.out
    time "1 h"

    input:
    tuple path(nuc), path(nuc_seqs), path(gb), path(gb_seqs)

    output:
    tuple path("manifest.csv"), path("sequences/*.fna.gz")

    script:
    """
    mkdir sequences && mv ${nuc_seqs} ${gb_seqs} sequences/
    merge.R manifest.csv nucleotide.csv genbank.csv
    """
}


process food_mappings {
    cpus 1
    memory "64 GB"
    publishDir "${params.out}/dbs"
    time "30 m"

    input:
    tuple path(foodb), path(gb_summary)
    path(matches)
    path(curated)

    output:
    tuple path("food_matches.csv"), path("food_contents.csv.gz")

    script:
    """
    food_mapping.R ${foodb} ${matches} ${curated}
    """
}

process sketch {
    cpus 2
    memory "16 GB"
    publishDir "${params.out}/sketches"
    time "8 h"

    input:
    path(seq)

    output:
    path("*.sig")

    script:
    """
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000 ${seq}
    """
}

process ANI {
    cpus 8
    memory "96 GB"
    publishDir "${params.out}", mode: "copy", overwite: true
    time "8 h"

    input:
    path(sigs)

    output:
    path("mash_ani.csv")

    script:
    """
    sourmash compare -k 21 --ani -p ${task.cpus} --csv mash_ani.csv ${sigs}
    """
}
