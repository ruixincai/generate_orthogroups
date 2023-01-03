log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

# the outputs will be returned in the directory called "separate_ortho"
outputdirectory <- snakemake@output[['output_dir']]
# if the outputdirectory doesn't exist, create the outputdirectory
if(!dir.exists(outputdirectory)) {
  dir.create(outputdirectory)
}

# use the library(data.table) to read "ortho_match.csv"
ortho <- fread(snakemake@input[["input_file"]])

# x[i, j] to substitute the prefix "^rna-" to "" in the column transcript_id
# the new transcript_id will be generated in the new column called new_transcript_id
ortho[,new_transcript_id := gsub("^rna-", "", transcript_id)]
# only do these in ortho which named Orthogroup
setkey(ortho, Orthogroup)

# to get a function of x (Orthogroup)
# firstly create the name of file with "x.txt"
# fill the file with the transcript_ids (in the column called "new_transcript_id" in ortho) 
# without the column name
printortho <- function(x) {
  file_name <- paste0(outputdirectory, "/", x, ".txt")
  fwrite(ortho[x, .(new_transcript_id)], file_name, col.names = FALSE)}

# loop all the Orthogroup column
all_orthog <- ortho[, unique(Orthogroup)]
lapply(all_orthog, printortho)

sessionInfo()
