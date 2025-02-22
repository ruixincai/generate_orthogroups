log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(rtracklayer)

#files <- list.files(path="/Users/cairuixin/gff", pattern="*.gff", full.names=TRUE, recursive=FALSE)
#print(files)

# LookupParentGene <- function(cds_parent){
# mc[ID == cds_parent, unique(Parent)]}

# import and read the gff
gff <- import.gff(snakemake@input[['gff']])



# split gff file to metadata columns
mc <- as.data.table(mcols(gff))

# only left the "CDS" type
cds <- mc[as.character(type) == "CDS"]

# get the parent transcript for each CDS
cds[, parent_char := as.character(Parent)]

cds_out <- unique(cds[, .(gene = gene,
                          protein_id = protein_id,
                          transcript_id = parent_char)])

fwrite(cds_out, snakemake@output[['output_tables']])

sessionInfo()