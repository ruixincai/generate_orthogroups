log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
# get the csv file names via path
file_names <- snakemake@input[['input_files']]
# names list names from file_names, replace ".csv" with ""
names(file_names) <- sub(".csv", "", basename(file_names))
# create an empty csv list including 3 file_names
csv_list <- lapply(file_names, fread)
# fill the csv list with all the data, name the first column "names_of_mosquitoes"
full_table <- rbindlist(csv_list, idcol = "species_of_mosquitoes")
fwrite(full_table, snakemake@output[["output_file"]])
sessionInfo()