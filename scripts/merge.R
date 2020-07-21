# merge a bunch of CSVs
# first argument is output file and the rest the files to merge

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
files <- args[-1]
print(files)

read <- lapply(files, fread)
fwrite(rbindlist(read), out)
