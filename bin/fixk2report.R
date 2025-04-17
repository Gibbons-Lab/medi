#!/usr/bin/env Rscript

# Fixes missing ranks in Kraken2 reports

library(data.table)
library(futile.logger)

domains <- c("Archaea", "Bacteria", "Eukaryota")

args <- commandArgs(trailingOnly = TRUE)
report <- fread(args[1], sep="\t", header=FALSE, strip.white=FALSE)

fix <- function(name) {
    name <- trimws(name)

    out <- "D1"
    if (name == "root") {
        out <- "R"
    } else if (name == "cellular organisms") {
        out <- "R1"
    } else if (name %in% domains) {
        out <- "D"
    }

    flog.info("Mapping %s -> rank: %s", name, out)
    return(out)
}

report[V4 == "" | (V4 %chin% as.character(0:9)), V4 := sapply(V6, fix)]
fwrite(report, args[2], sep="\t", quote=FALSE, col.names=FALSE)