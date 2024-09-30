#!/usr/bin/env Rscript

# Fixes missing ranks in Kraken2 reports

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
report <- fread(args[1], sep="\t", header=FALSE, strip.white=FALSE)

fix <- function(name) {
    name <- trimws(name)
    if (name == "root") {
        return("R")
    } else if (name == "cellular organisms") {
        return("R1")
    } else {
        return("D")
    }
}

report[V4 == "" | V4 == "1", V4 := sapply(V6, fix)]
report[V4 == "1", V4 := "D1"]
fwrite(report, args[2], sep="\t", quote=FALSE, col.names=FALSE)