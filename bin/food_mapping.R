#!/usr/bin/env Rscript

# Creates a mapping of Food content to taxonomy IDs

library(data.table)
args <- commandArgs(TRUE)

clean_unit <- function(x) {
    x <- gsub("/ ", "/", x)
    x <- gsub(" g", "g", x)
    x <- gsub("\\s\\w{2,}.+", "", x)
    return(x)
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

energy <- function(ab) {
    # 1 - Fat, 2 - Proteins, 3 - Carbs, 4 - Fatty acids
    amounts <- c()
    amounts["fat"] <- 9 * max(
        ab[compound_id == 1, mean(standard_content, na.rm=T) / 1000 ],
        ab[compound_id == 4, mean(standard_content, na.rm=T) / 1000 ]
    )
    amounts["protein"] <- 4 * ab[compound_id == 2, mean(standard_content, na.rm=T) / 1000 ]
    amounts["carbs"] <- 4 * ab[compound_id == 3, mean(standard_content, na.rm=T) / 1000 ]
    if (length(amounts) < 3) {
        kcal <- NA
    } else {
        kcal <- sum(amounts)
    }
    res <- data.table(
        compound_id = 200000,
        orig_unit = "kcal/100g",
        source_type = "Nutrient",
        standard_content = kcal,
        orig_content = kcal,
        orig_min = kcal,
        orig_max = kcal
    )
}

food <- fread(args[3])
content <- fread(file.path(args[1], "Content.csv"))
compounds <- fread(file.path(args[1], "Compound.csv"))
nutrients <- fread(file.path(args[1], "Nutrient.csv"))
matched <- unique(fread(args[2])[,
    .(orig_taxid, db, rank, matched_taxid, kingdom, phylum,
      class, order, family, genus, species)])

food <- food[!is.na(revised_taxonomy_id),
    .(food_id = id, wikipedia_id, food_group, food_subgroup,
      revised_taxonomy_id)]
content <- content[!is.na(standard_content), .(
    compound_id = source_id, food_id, orig_content, orig_min, orig_max,
    orig_unit, standard_content, preparation_type, source_type
)]

# Add calculated energy
energies <- content[, energy(.SD), by=c("food_id", "preparation_type")]
content <- rbind(content, energies, fill=TRUE)

unit_map <- content[, unique(orig_unit) |> sapply(clean_unit)]
content[, "unit" := unit_map[orig_unit]]
content[, "scale" := scales[unit]]
content[, standard_content := standard_content * scale]
content[scale != 1, unit := "mg/100g"]

compounds <- compounds[, .(
    compound_id = id, name, description = annotation_quality,
    CAS = description, mono_mass = moldb_inchi, compound_kingdom = kingdom,
    compound_superclass = superklass, comound_class = klass, compound_subclass = subklass,
    source_type = "Compound"
)]
nutrients <- nutrients[, .(
    compound_id = id, name, description, source_type = "Nutrient"
)]
nutrients <- rbind(
    nutrients,
    data.table(
        compound_id = 200000,
        name = "Energy (calculated)",
        source_type = "Nutrient",
        description = "Energy content calculated from mean macronutrients amounts."
    )
)
compounds <- rbind(compounds, nutrients, use.names = TRUE, fill = TRUE)
food_matches <- food[
    matched, on = c(revised_taxonomy_id = "orig_taxid"), nomatch = 0]
fwrite(food_matches, "food_matches.csv")
food_matches <- content[food_matches, on = "food_id", nomatch = 0]
food_matches <- compounds[
    food_matches, on = c("compound_id", "source_type"), nomatch = 0]

# remove too high calory values
food_matches <- food_matches[!(unit == "kcal/100g" & orig_content > 900)]
# Correct the calories for which the standard_content contains incorrect values,
# but orig_content seems to be correct
food_matches[unit == "kcal/100g", standard_content := orig_content]

# Correct incorrect cholesterol units.
# There are two clear peaks in the histogram with one peak > 1000mg/100g
#  (probably ug -> mg error in FooDB)
food_matches[
    unit == "mg/100g" & name == "Cholesterol" & standard_content > 1000,
    standard_content := standard_content / 1000.0
]

# Only keep kcal for energy
# and remove obviously incorrect entries like 200g nutrient/100g food
food_matches <- food_matches[
    !(name == "Energy" & unit != "kcal/100g")][
    !(unit == "mg/100g" & standard_content > 100000)
    ]

fwrite(food_matches, "food_contents.csv.gz")
