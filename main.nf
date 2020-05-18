// Build the FOODB genome database

foodb = file(
    "https://foodb.ca/public/system/downloads/foodb_2020_4_7_csv.tar.gz")
refseq_catalog = file(
    "https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/RefSeq-release200.catalog.gz"
)
genbank_summary = file(
    "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
)


process download {

    output:
    path("food_taxids.txt"), path("genbank_taxids.txt"), path("refseq_taxids.txt"),
        path("foodb"), path("genbank_summary.tsv"), path("refseq_catalog.tsv.gz") into tax_dbs, dbs

    """
    wget $foodb -O foodb.tgz \
    wget $refseq_catalog -O refseq_catalog.tsv.gz \
    wget $genbank_summary -O genbank_summary.tsv \
    tar -xf foodb.tgz \
    zcat refseq_catalog.tsv.gz | cut -f 1 > refseq_taxids.txt \
    cut -f 6 genbank_summary.tsv | tail -n +3 > genbank_taxids.txt \
    cut -d "," -f 3 foodb/Foods.csv > food_taxids.txt
    """
}




