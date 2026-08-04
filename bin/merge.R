#!/usr/bin/env Rscript

# merge a bunch of CSVs
# first argument is output file and the rest the files to merge

library(data.table)

read_medi = function(filepath) {
    tab <- fread(filepath)
    if ("id" %in% names(tab)) {
        tab[, "id" := as.character(id)]
    }

    return(tab)
}

args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
files <- args[-1]
print(files)

read <- lapply(files, read_medi)
fwrite(rbindlist(read, use.names=TRUE, fill=TRUE), out)
