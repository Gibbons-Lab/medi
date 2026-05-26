#!/usr/bin/env Rscript

library(data.table)
library(reutils)
library(magrittr)
library(futile.logger)
library(Biostrings)
library(R.utils)

MAX_SEQLENGTH = 2e9

args <- commandArgs(trailingOnly = TRUE)

matches <- fread(args[1])
group <- args[2]
out_folder <- args[3]
target_id <- args[4]


if (is.null(getOption("reutils.api.key"))) {
    rate <- 0.9
} else {
    rate <- 9
    api_key <- getOption("reutils.api.key")
}

dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)

ncbi_download <- function(accession, filename) {
    data_package <- paste0(accession, ".zip")
    ret <- system2(
        "datasets",
        c("download", "genome", "accession", accession,
          "--include", "genome",
          "--filename", paste0(accession, ".zip"))
    )

    if (ret != 0) return(ret)

    unzip(
        data_package,
        paste0("ncbi_dataset/data/", accession, "/", filename),
        junkpaths=TRUE
    )
    unlink(data_package)

    return(ret)
}

download_genome <- function(hit, out_dir="sequences") {
    hit <- copy(hit[1])
    id <- basename(hit$url)
    hit$filename <- file.path(out_dir, paste0(id, ".fna"))
    for (i in 0:7) {
        if (file.exists(hit$filename)) unlink(hit$filename)
        ret <- tryCatch(
            ncbi_download(hit$id, hit$filename),
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
    flog.info("Downloaded genome for assembly %s...", id)
    fa <- readDNAStringSet(hit$filename)
    short_names <- tstrsplit(names(fa), "\\s+")[[1]]
    names(fa) <- paste0(short_names, "_", 1:length(short_names),
                        "|kraken:taxid|", as.character(hit$matched_taxid),
                        " ", names(fa))
    writeXStringSet(fa, hit$filename, compress = "gzip")
    hit$num_records <- length(fa)
    hit$seqlength <- as.double(sum(width(fa)))

    return(hit)
}

fragmented_efetch <- function(hits, taxid, filename) {
    hits[, "group" := floor(cumsum(seqlength) / MAX_SEQLENGTH)]
    if (file.exists(filename)) unlink(filename)
    flog.info("Downloading sequences for taxon %s [fragmented into %d groups]...",
        taxid, hits[, uniqueN(group)])

    for (g in unique(hits$group)) {
        for (i in 0:7) {
            Sys.sleep(1/rate + 2^i)
            post <- epost(hits[group == g, id], db = "nuccore")
            Sys.sleep(1/rate)
            fetch <- suppressMessages(
                efetch(post, db = "nuccore",
                       rettype = "fasta", retmode = "text")
            )
            if (length(getError(fetch)) == 1) {
                write(content(fetch), filename, append=TRUE)
                if (file.exists(filename) && grepl(">", content(fetch))) {
                    done <- TRUE
                    break
                }
            }
        }
        if (i == 7) {
            if (file.exists(filename)) unlink(filename)
            done <- FALSE
            break
        }
    }

    return(done)
}

download_sequences <- function(hits, taxid, out_dir="sequences") {
    hits <- copy(hits) %>% unique(by="id")
    filename <- file.path(out_dir, paste0(as.character(taxid), ".fna"))
    done <- fragmented_efetch(hits, taxid, filename)
    if (!done) {
        flog.error("Failed downloading %s. UIDs=%s) :(", taxid, paste(unique(hits$id), collapse=", "))
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
if ((args[2] == "nucleotide") && (any(matches$db == "nucleotide"))) {
    contigs <- matches[
        db == "nucleotide",
        download_sequences(.SD, matched_taxid[1]),
        by = "matched_taxid"]
    flog.info("Downloaded contigs for %d additional taxa.", nrow(contigs))
    contigs[, "orig_taxid" := NULL]
    fwrite(contigs, "nucleotide.csv")
}

# Download full genomes
if (args[2] == "genbank") {
    target = args[4]
    gb <- matches[db == "genbank"]
    flog.info("Downloading genome %s with.", target)
    report <- download_genome(gb[id == target])
    flog.info(
        "Found %d records summing to %.3g Mbps.",
        report$num_records,
        report$seqlength / 1e6
    )

    report[, "orig_taxid" := NULL]
    fwrite(report, paste0(target, ".csv"))
}

if (args[2] == "decoys") {
    flog.info("Downloading %d additional decoys.", nrow(matches))
    decoys <- matches[, download_genome(.SD, out_folder), by = "id"]
    flog.info("Downloaded genomes for %d decoys summing to %.2g Mbps.",
                nrow(decoys), decoys[, sum(seqlength) / 1e6])
    fwrite(decoys, "decoys.csv")
}
