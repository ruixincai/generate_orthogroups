log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")
library(data.table)
# or use library(dplyr)
# dplyr equivalent: read_csv
og_wide <- fread(snakemake@input['input_orthofinder'])
# convert to a *long* data.table
og_long <- melt(og_wide, id.vars = "Orthogroup", 
                variable.name = "species", 
                value.name = "all_proteins")
setkey(og_long, Orthogroup)
# split the all_proteins column
split_orthogroups <- og_long[,.(protein_id = unlist(strsplit(all_proteins, ", "))),
                             by = .(Orthogroup, species)]
# read the combined table
protein_to_gene <- fread(snakemake@input[["input_combinedtable"]])
og_with_tx <- merge(split_orthogroups,
                    protein_to_gene[, .(protein_id, transcript_id)],
                    by = "protein_id",
                    all.x = TRUE,
                    all.y = FALSE)
setkey(og_with_tx, Orthogroup, species, protein_id, transcript_id)
fwrite(og_with_tx, snakemake@output[["output_file"]])
# write a text file for each orthogroup (will be 13,982 text files)
# one line per transcript id
sessionInfo()