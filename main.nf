// Build the FOODB genome database

foodb = "https://foodb.ca/public/system/downloads/foodb_2020_4_7_csv.tar.gz"
refseq_summary = "https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/release200.accession2geneid.gz"
genbank_summary = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
taxdump = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"


process download {
    cpus 1
    publishDir "$baseDir/data/dbs"

    output:
    tuple  path("foodb"), path("genbank_summary.tsv"), path("refseq_summary.tsv.gz") into dbs

    """
    wget $foodb -O foodb.tgz && \
    wget $refseq_summary -O refseq_summary.tsv.gz && \
    wget $genbank_summary -O genbank_summary.tsv && \
    tar -xf foodb.tgz && \
    mv foodb_2020_04_07_csv foodb
    """
}

process get_taxids {
    cpus 1

    input:
    tuple path(foodb), path(gb_summary), path(refseq_cat) from dbs

    output:
    tuple path("foodb"), path("taxids.tsv"), path("$gb_summary") into taxids

    """
    #!/usr/bin/env Rscript

    library(data.table)

    dt <- fread("zcat < $refseq_cat", sep="\t", header=F)
    refseq <- dt[!is.na(V1), .(taxid = as.character(unique(V1)))]
    refseq[, "source" := "refseq"]
    dt <- fread("$genbank_summary", sep="\t")
    genbank <- dt[!is.na(taxid), .(taxid = as.character(unique(taxid)))]
    genbank[, "source" := "genbank"]
    dt <- fread("$foodb/Food.csv")
    foodb <- dt[!is.na(ncbi_taxonomy_id), .(taxid = ncbi_taxonomy_id)]
    foodb[, "source" := "foodb"]
    fwrite(rbind(genbank, refseq, foodb), "taxids.tsv", col.names=F, sep="\t")
    """
}

process download_taxa_dbs {
    cpus 1

    output:
    path("taxdump") into taxa_db

    """
    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz && \
    mkdir taxdump && tar -xf taxdump.tar.gz --directory taxdump
    """
}

process get_lineage {
    cpus 1

    input:
    tuple path(foodb), path(taxids), path(gb_summary), path(taxadb) from taxids.combine(taxa_db)

    output:
    tuple path("$foodb"), path("lineage.txt"), path("lineage_ids.txt"), path("$gb_summary") into lineage

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
    tuple path(foodb), path(lineage), path(lineage_ids), path(gb_summary) from lineage

    output:
    path("matches.csv") into matches

    """
    Rscript $baseDir/scripts/match.R $lineage_ids $gb_summary
    """
}

process download_sequences {
    cpus 1
    publishDir "$baseDir/data/sequences"

    input:
    path(matches) from matches

    output:
    path("sequences/*.fna.gz") into sequences
    path("manifest.csv") into manifest

    """
    Rscript $baseDir/scripts/download.R $matches $task.cpus sequences
    """
}

process setup_kraken_db {
    publishDir "$baseDir/data"

    output:
    path("kraken2_db") into db

    """
    kraken2-build --download-taxonomy --db kraken2_db --use-ftp
    """
}

process add_sequences {
    cpus 5

    input:
    path(fasta) from sequences.flatMap()
    each path(db) from db

    output:
    tuple val("$fasta"), path("$db") into added_db

    """
    gunzip -c $fasta > ${fasta.baseName} && \
    kraken2-build --add-to-library ${fasta.baseName} --db $db --threads ${task.cpus} && \
    rm ${fasta.baseName}
    """
}

process build_kraken_db {
    cpus 20

    input:
    tuple val(fasta), path(db) from added_db.last()

    output:
    path("$db") into built_db

    """
    kraken2-build --build --db $db --threads ${task.cpus}
    """
}

process setup_bracken {
    cpus 20

    input:
    path(db) from built_db

    output:
    path("$db") into bracken_db

    """
    bracken-build -d $db -t ${task.cpus} -k 35 -l 100
    """
}
