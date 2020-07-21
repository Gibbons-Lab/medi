# Quantify the food abundances and contents

library(data.table)

join <- function(x) {
    paste(unique(x), collapse = "|")
}

args <- commandArgs(trailingOnly = TRUE)
foods <- fread(args[1])
contents <- fread(args[2])
files <- args[-(1:2)]
species_file <- files[grepl("^S_", files)]
genus_file <- files[grepl("^G_", files)]
id <- gsub("^S_", "", strsplit(basename(species_file), "\\.b2")[[1]][1])
print(id)
abundances <- rbind(fread(species_file, sep="\t"),
                    fread(genus_file, sep="\t"))
abundances <- abundances[, .(
    taxon = name, matched_taxid = taxonomy_id, reads = new_est_reads,
    fraction = new_est_reads / sum(new_est_reads),
    sample_id = id)]
food_quant <- foods[abundances, on = "matched_taxid", nomatch = 0]
fwrite(food_quant, sprintf("%s_food_abundance.csv", id))

contents <- contents[abundances, on = "matched_taxid", nomatch = 0]
contents[, "reps" := .N, by = c("matched_taxid", "compound_id", "source_type")]
food_content <- contents[, .(
    abundance = sum(fraction * standard_content / reps),
    name = name[1],
    monomer_mass = mono_mass[1],
    kingdom = kingdom[1],
    superclass = superclass[1],
    class = class[1],
    subclass = subclass[1],
    CAS = CAS[1]
    ),
    by = c("compound_id", "source_type")
]
food_content[, "relative" := abundance / sum(abundance), by = "source_type"]
food_content[, "sample_id" := id]
fwrite(food_content, sprintf("%s_food_content.csv", id))
