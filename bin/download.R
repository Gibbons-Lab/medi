#!/usr/bin/env Rscript

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
    rate <- 0.9
} else {
    rate <- 9
}

dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)

ncbi_rsync <- function(url, out) {
    rsync_url <- gsub("https://", "rsync://", url)
    ret <- system2("rsync", c("--no-motd", rsync_url, out))
    Sys.chmod(out, "0755")
    return(ret)
}

download_genome <- function(hit, out_dir="sequences") {
    hit <- copy(hit[1])
    id <- basename(hit$url)
    hit$url <- paste0(hit$url, "/", id, "_genomic.fna.gz")
    hit$filename <- file.path(out_dir, paste0(id, ".fna.gz"))
    for (i in 0:7) {
        if (file.exists(hit$filename)) unlink(hit$filename)
        ret <- tryCatch(
            ncbi_rsync(hit$url, hit$filename),
            error = function(e) return(1),
            warning = function(e) return(1)
        )
        if (ret == 0) break
        Sys.sleep(2^i)
    }
    if (ret != 0) {
        flog.error("Failed downloading %s :(", hit$url)
        stop()
    }
    fa <- readDNAStringSet(hit$filename)
    short_names <- tstrsplit(names(fa), "\\s+")[[1]]
    names(fa) <- paste0(short_names, "_", 1:length(short_names),
                        "|kraken:taxid|", as.character(hit$matched_taxid),
                        " ", names(fa))
    writeXStringSet(fa, hit$filename, compress = "gzip")
    hit$num_records <- length(fa)
    hit$seqlength <- as.double(sum(width(fa)))
    flog.info("Downloaded genome for assembly %s...", id)
    return(hit)
}

download_sequences <- function(hits, taxid, out_dir="sequences") {
    hits <- copy(hits)
    filename <- file.path(out_dir, paste0(as.character(taxid), ".fna"))
    flog.info("Downloading sequences for taxon %s...", taxid)
    for (i in 0:7) {
        Sys.sleep(1/rate + 2^i)
        if (file.exists(filename)) unlink(filename)
        post <- epost(unique(hits$id), db = "nuccore")
        Sys.sleep(1/rate)
        fetch <- suppressMessages(
                efetch(post, db = "nuccore",
                       rettype = "fasta", retmode = "text")
        )
        if (length(getError(fetch)) == 1) {
            write(content(fetch), filename)
            if (file.exists(filename) && grepl(">", content(fetch))) {
                break
            }
        }
    }
    if (!file.exists(filename) || !grepl(">", content(fetch))) {
        flog.error("Failed downloading %s. UIDs=%s) :(", taxid, paste(unique(hits$id), collpase=", "))
        print(post)
        stop()
    }
    hit <- hits[1]
    hit$filename <- paste0(filename, ".gz")
    fa <- readDNAStringSet(filename)
    short_names <- tstrsplit(names(fa), "\\s+")[[1]]
    names(fa) <- paste0(short_names, "_", 1:length(short_names),
                        "|kraken:taxid|", as.character(taxid),
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
gb <- matches[db == "genbank"]
flog.info("Downloading %d genomes with %d threads.", gb[, uniqueN(id)], threads)
dls <- parallel::mclapply(
    gb[, unique(id)],
    function(i) download_genome(gb[id == i]),
    mc.cores=threads
)
genomes <- rbindlist(dls)
flog.info("Downloaded %d full genomes.", nrow(genomes))

manifest <- rbind(genomes, contigs)
manifest[, "orig_taxid" := NULL]
fwrite(manifest, "manifest.csv")

