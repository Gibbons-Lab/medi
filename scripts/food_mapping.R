# Creates a mapping of Food content to taxonomy IDs

library(data.table)
args <- commandArgs(TRUE)

clean_unit <- function(x) {
    x <- gsub("/ ", "/", x)
    x <- gsub(" g", "g", x)
    x <- gsub("\\s\\w{2,}.+", "", x)
}

scales <- c(
    `mg/100g` = 1,
    `kcal/100g` = 1,
    `RE` = 1,
    `α-TE` = 1,
    `NE` = 1,
    `IU` = 1,
    `µg` = 1,
    `ug/g` = 0.1,
    `uM` = 0.1,
    `mg/kg` = 0.1,
    `ug/kg` = 0.0001,
    `ug/L` = 0.0001,
    `umol/g` = 1,
    `mg/l` = 0.1,
    `mg/g` = 0.01,
    `ppb` = 1,
    `ug/100g` = 0.001,
    `g/kg` = 100,
    `IU/100g` = 1
)

food <- fread(file.path(args[1], "Food.csv"))
content <- fread(file.path(args[1], "Content.csv"))
compounds <- fread(file.path(args[1], "Compound.csv"))
nutrients <- fread(file.path(args[1], "Nutrient.csv"))
matched <- unique(fread(args[2])[,
    .(orig_taxid, db, rank, matched_taxid, kingdom, phylum,
      class, order, family, genus, species)])

food <- food[!is.na(ncbi_taxonomy_id),
    .(food_id = id, wikipedia_id, food_group, food_subgroup,
      ncbi_taxonomy_id)]
content <- content[!is.na(standard_content), .(
    compound_id = source_id, food_id, orig_content, orig_min, orig_max,
    orig_unit, standard_content, preparation_type, source_type
)]
unit_map <- content[, unique(orig_unit) |> sapply(clean_unit)]
content[, "unit" := unit_map[orig_unit]]
content[, "scale" := scales[unit]]
content[, standard_content := standard_content * scale]
content[scale != 1, unit := "mg/100g"]
compounds <- compounds[, .(
    compound_id = id, name, description = annotation_quality,
    CAS = description, mono_mass = moldb_inchi, kingdom = kingdom,
    superclass = superklass, class = klass, subclass = subklass,
    source_type = "Compound"
)]
nutrients <- nutrients[, .(
    compound_id = id, name, description, source_type = "Nutrient"
)]
compounds <- rbind(compounds, nutrients, use.names = TRUE, fill = TRUE)
food_matches <- food[
    matched, on = c(ncbi_taxonomy_id = "orig_taxid"), nomatch = 0]
fwrite(food_matches, "food_matches.csv")
food_matches <- content[food_matches, on = "food_id", nomatch = 0]
food_matches <- compounds[
    food_matches, on = c("compound_id", "source_type"), nomatch = 0]
food_matches <- food_matches[
    !(name == "Energy" & unit != "kcal/100g")][
    !(unit == "mg/100g" & standard_content > 100000)
    ]
fwrite(food_matches, "food_contents.csv.gz")
