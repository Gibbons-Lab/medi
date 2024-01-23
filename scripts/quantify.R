# Quantify the food abundances and contents

library(data.table)

join <- function(x) {
    paste(unique(x), collapse = "|")
}

args <- commandArgs(trailingOnly = TRUE)
foods <- fread(args[1])
contents <- fread(args[2])[standard_content > 0]
files <- args[-(1:2)]
species <- files[grepl("^S_", files)] |> fread()
species[, "taxrank" := "species"]
found_genera <- species[, tstrsplit(taxid_lineage, ";")[[6]] |> as.numeric()]
genus <- files[grepl("^G_", files)] |> fread()
genus <- genus[!taxonomy_id %chin% found_genera]
genus[, "taxrank" := "genus"]
abundances <- rbind(species, genus)
abundances <- abundances[, .(
    taxon = name, taxrank = taxrank,  matched_taxid = taxonomy_id,
    kraken_raw_reads = kraken_assigned_reads,
    reads = new_est_reads,
    sample_id = gsub("^([A-Z])_", "", sample_id, perl=T),
    lineage = lineage, taxid_lineage = taxid_lineage)]
abundances[, "relative_abundance" := reads / sum(reads), by="sample_id"]
abundances[, "total_reads" := sum(reads), by="sample_id"]
abundances[, "total_raw_reads" := sum(kraken_raw_reads), by="sample_id"]
abundances[
    ,
    "bacteria_reads" := sum(reads[grepl("k__Bacteria;", lineage)], na.rm = TRUE),
    by="sample_id"
]
abundances[
    ,
    "human_reads" := sum(reads[grepl("g__Homo;", lineage)], na.rm = TRUE),
    by="sample_id"
]


food_quant <- foods[abundances, on = "matched_taxid", nomatch = 0]
info <- abundances[
    taxrank == "species",
    .(sample_id, total_reads, total_raw_reads, bacteria_reads, human_reads)
] |> unique() |> setkey(sample_id)
missing <- info[, unique(sample_id[!sample_id %chin% food_quant$sample_id])]
empty <- lapply(
    missing,
    function(sid) data.table(
        taxon = NA, taxrank = "species", matched_taxid = NA, kraken_raw_reads = 0,
        reads = 0, relative_abundance = 0, sample_id = sid,
        lineage = NA, taxid_lineage = NA,
        total_reads = info[sid, total_reads],
        total_raw_reads = info[sid, total_raw_reads],
        bacteria_reads = info[sid, bacteria_reads],
        human_reads = info[sid, human_reads]
    )
)
summary(sapply(empty, nrow)) |> print()
food_quant <- rbind(food_quant, rbindlist(empty), fill=TRUE)

fwrite(food_quant, "food_abundance.csv")

contents <- contents[abundances, on = "matched_taxid", allow.cartesian = TRUE, nomatch = 0]
contents[, "reps" := .N, by = c("sample_id", "matched_taxid", "compound_id", "source_type", "unit")]
contents[, "relative_food_abundance" := reads / sum(reads[!duplicated(matched_taxid)], na.rm = T), by="sample_id"]
food_content <- contents[, .(
    abundance = sum(relative_food_abundance * standard_content / reps),
    name = name[1],
    monomer_mass = mono_mass[1],
    kingdom = kingdom[1],
    superclass = superclass[1],
    class = class[1],
    subclass = subclass[1],
    CAS = CAS[1]
    ),
    by = c("sample_id", "compound_id", "source_type", "unit")
]
food_content[unit == "mg/100g", "relative_abundance" := abundance / sum(abundance), by = c("source_type", "sample_id")]
empty_content <- data.table(sample_id = missing)
food_content <- rbind(food_content, empty_content, fill=TRUE)

fwrite(food_content, "food_content.csv")
