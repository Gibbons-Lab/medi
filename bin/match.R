#!/usr/bin/env Rscript

# helper functions to match based on phylogenetic ranks

library(data.table)
library(R.utils)
library(magrittr)
library(futile.logger)
library(reutils)

ASSEMBLY_SEARCH <- "txid%s[orgn]"
NT_SEARCH <- "txid%s[orgn] AND 10000:10000000000[SLEN] AND biomol_genomic[PROP]"
GB_CATEGORY_SCORES <- c(`na` = 0, `representative genome` = 10, `reference genome` = 20)
GB_ASSEMBLY_SCORES <- c(Contig = 0, Chromosome = 1, Scaffold = 1, `Complete Genome` = 2)


if (is.null(getOption("reutils.api.key"))) {
    rate <- 0.9
} else {
    rate <- 9
}

genbank_quality <- function(dt) {
    dt <- dt[order(-seq_rel_date)]
    dt[, "score" := 0]
    dt[, score := score + GB_CATEGORY_SCORES[refseq_category] + GB_ASSEMBLY_SCORES[assembly_level]]
    dt[is.na(score), "score" := 0]

    return(dt)
}

not_found <- function(res) {
    if (is.character(getError(res)$wrnmsg)) {
        return(getError(res)$wrnmsg == "No items found.")
    }
    return(FALSE)
}

find_taxon <- function(taxid, gb_taxa, gb_summary, col, db) {
    url <- NULL
    taxid <- as.character(taxid)[!is.na(taxid)]
    if (db == "genbank") {
        flog.info("Querying the assembly database for taxon %s...", taxid)
        tids <- gb_taxa[as.character(taxid), orig_taxid]
        if (all(is.na(tids))) {
            return(NULL)
        }
        matches <- gb_summary[tids]
        # If we have a complete genome we use the most recent one
        if (matches[, max(score)] > 1) {
            matches <- matches[score == max(score)][1]
        } else {
            matches <- matches[score == max(score)]
        }
        uids <- matches[, `#assembly_accession`]
        url <- matches[, ftp_path]
        refseq_category <- matches[, refseq_category]
        assembly_level <- matches[, assembly_level]
        seqlength <- matches$genome_size
        genome_type <- "full genome"
        name <- paste(matches$organism_name, matches$infraspecific_name)
    } else {
        r <- rate
        for (i in 0:7) {
            Sys.sleep(1 / rate + 2^i)
            tids <- as.character(taxid)
            flog.info("Querying the nt database for taxon %s...", taxid)
            ret <- suppressMessages(esearch(
                sprintf(NT_SEARCH, as.character(taxid)),
                retmax = 500,
                sort = "SLEN",
                db = "nuccore"
            ))
            Sys.sleep(1/rate)
            summ <- suppressMessages(
                esummary(ret, db="nuccore") %>% content("parsed") %>% data.table(fill=T)
            )
            summ[, "Slen" := as.integer(Slen)]
            if (ret$no_errors() || not_found(ret)) {
                break
            }
            if (i == 7) {
                flog.info("Querying failed for %s. Aborting.", taxid)
                stop()
            }
        }
        uids <- ret %>% uid()
        uids <- uids[!is.na(uids)]
        refseq_category <- "excluded"
        assembly_level <- "contig"
        seqlength <- summ$Slen
        genome_type <- summ$Genome
        name <- summ$Title
    }
    if (length(uids) == 0) {
        return(NULL)
    }
    return(data.table(
        id = uids, db = db, matched_taxid = taxid, url = url,
        refseq_category = refseq_category, assembly_level = assembly_level,
        seqlength=seqlength, name = name, genome_type = genome_type
    ))
}

ordered_match <- function(
    query, taxa, summary, db = "genbank", rank = "species") {
    coln <- paste0(rank, "_taxid")
    flog.info("Searching with rank `%s` in `%s`...", rank, db)
    query[["group"]] <- query[[coln]]
    taxa[["group"]] <- taxa[[coln]]
    taxa[, "group" := as.character(group)]
    setkey(taxa, "group")
    if (nrow(query) == 0) {
        return(NULL)
    }
    good <- query[, !is.na(group) & group != ""]
    hits <- query[
        good,
        find_taxon(group[1], taxa, summary, coln, db),
        by = "orig_taxid"
    ]
    if (nrow(hits) == 0) {
        return(NULL)
    }
    hits[["rank"]] <- rank

    return(hits)
}

RANKS <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
args <- commandArgs(trailingOnly = TRUE)

taxids <- fread(args[1],
    sep = "\t",
    col.names = c("orig_taxid", "source", "lineage", "names", "taxids")
)

gbs <- fread(args[2], sep = "\t")
# Some assemblies seem to have new existing release in never Genbank versions.
# Only keep the ones actually available.
gbs <- gbs[grepl("ftp.ncbi.nlm.nih.gov", ftp_path, fixed = TRUE)]
gbs[, "taxid" := as.character(taxid)]
gbs <- genbank_quality(gbs)
setkey(gbs, taxid)
taxids[, (RANKS) := tstrsplit(names, ";")]
taxids[, (paste0(RANKS, "_taxid")) := tstrsplit(taxids, ";")]
taxids[, "orig_taxid" := as.character(orig_taxid)]
food <- taxids[source == "foodb"]
gb_taxa <- taxids[source == "genbank"]

trials <- data.table(
    rank = rep(c("species", "genus"), 2),
    db = rep(c("genbank", "nucleotide"), each = 2)
)
flog.info(
    "Trying the following matching strategy:\n%s",
    capture.output(trials)
)
matches <- NULL
for (i in 1:nrow(trials)) {
    rank <- trials[i, rank]
    db <- trials[i, db]
    queries <- food[!orig_taxid %in% matches$orig_taxid]
    m <- ordered_match(queries, gb_taxa, gbs, db = db, rank = rank)
    matches <- rbind(matches, m, fill = TRUE, use.names = TRUE)
    gbs[`#assembly_accession` %chin% matches$id, score := score - 30]
}
matches <- unique(food[, c("orig_taxid", RANKS), with = FALSE])[
    matches,
    on = "orig_taxid"
]
fwrite(matches, "matches.csv")

flog.info(
    "Matched %d/%d taxa.",
    matches[, uniqueN(orig_taxid)], food[, uniqueN(orig_taxid)]
)
flog.info("Matches by strategy:\n")
print(matches[, uniqueN(id), by = c("db", "rank")])
