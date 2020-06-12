library(data.table)
library(reutils)
library(magrittr)
library(futile.logger)
library(Biostrings)
library(R.utils)

args <- commandArgs(trailingOnly = TRUE)

matches <- fread(args[1])
threads <- as.numeric(args[2])
out_folder <- args[3]

if (is.null(getOption("reutils.api.key"))) {
    rate <- 1
} else {
    rate <- 10
}

dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)

download_genome <- function(hit, out_dir="sequences") {
    hit <- copy(hit[1])
    id <- basename(hit$url)
    hit$url <- paste0(hit$url, "/", id, "_genomic.fna.gz")
    hit$filename <- file.path(out_dir, paste0(id, ".fna.gz"))
    flog.info("Downloading genome for assembly %s...", id)
    for (i in 1:10) {
        if (file.exists(hit$filename)) unlink(hit$filename)
        ret <- tryCatch(
            download.file(hit$url, hit$filename, quiet = TRUE),
            error = function(e) return(1),
            warning = function(e) return(1)
        )
        if (ret == 0) break
    }
    if (ret != 0) {
        flog.info("Failed downloading %s :(", hits$url)
        return(NULL)
    }
    fa <- readDNAStringSet(hit$filename)
    names(fa) <- paste0("kraken:taxid|", as.character(hit$matched_taxid),
                        " ", names(fa))
    writeXStringSet(fa, hit$filename, compress = "gzip")
    hit$num_records <- length(fa)
    hit$seqlength <- as.double(sum(width(fa)))
    return(hit)
}

download_sequences <- function(hits, taxid, out_dir="sequences") {
    hits <- copy(hits)
    filename <- file.path(out_dir, paste0(as.character(taxid), ".fna"))
    flog.info("Downloading sequences for taxon %s...", taxid)
    post <- epost(unique(hits$id), db = "nuccore")
    for (i in 1:10) {
        Sys.sleep(1/rate)
        if (file.exists(filename)) unlink(filename)
        fetch <- suppressMessages(
                efetch(post, db = "nuccore",
                       rettype = "fasta", retmode = "text")
        )
        if (length(getError(fetch)) == 1) {
            write(content(fetch), filename)
            break
        }
    }
    if (!file.exists(filename) || !grepl(">", content(fetch))) {
        flog.info("Failed downloading %s :(", taxid)
        return(NULL)
    }
    hit <- hits[1]
    hit$filename <- paste0(filename, ".gz")
    fa <- readDNAStringSet(filename)
    names(fa) <- paste0("kraken:taxid|", as.character(taxid),
                        " ", names(fa))
    writeXStringSet(fa, hit$filename, compress = "gzip")
    unlink(filename)

    hit$num_records <- length(fa)
    hit$seqlength <- as.double(sum(width(fa)))
    return(hit)
}

# Download additional contigs
if (any(matches$db == "nucleotide")) {
    contigs <- matches[
        db == "nucleotide",
        download_sequences(.SD, matched_taxid[1]),
        by = "matched_taxid"]
    flog.info("Downloaded contigs for %d additional taxa.", nrow(contigs))
}

# Download full genomes
genomes <- matches[db == "genbank", download_genome(.SD), by= "id"]
flog.info("Downloaded %d full genomes.", nrow(genomes))


manifest <- rbind(genomes, contigs)
manifest[, "orig_taxid" := NULL]
fwrite(manifest, "manifest.csv")

