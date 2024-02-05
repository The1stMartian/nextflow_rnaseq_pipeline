# Purpose is to annotate read counts files produced by DESeq2
# Script will do a left join (merge) of columns in a gtf file (annotations)
# with the DESeq2 analysis file. Note that both tables must have one column
# to use a key column for matching.

# User customize:
keyColumn = "gene_id"
inputFolder = "YOUR INPUT FOLDER PATH"
deFile = "YOUR DESEQ FILE (SHOULD BE IN THE INPUT FOLDER)"
gtfFile = "Your GTF FILE PATH"

#BiocManager::install("rtracklayer")
library(rtracklayer)
library(dplyr)

# Open annotations file
annotations = readGFF(gtfFile)

# Remove duplicate annotations
colnames(annotations)

# remove exon entries
annotations = annotations[annotations$type != "exon",] 

# remove duplicate key column entries 
annotations = annotations[!duplicated(annotations[keyColumn]),]

# remove unnecessary columns (user customize unwanted column names)
annotations = subset(annotations, select=-c(note, orig_transcript_id, product, 
  transcript_biotype, exon_number, orig_protein_id, transl_table, codons,
  partial, inference, exception, part, pseudo, transl_except, protein_id,
  source,	type, score, phase, transcript_id, gbkey))
head(annotations)

# create out file
annotatedOutFileName = paste0(inputFolder, "Annotated_", deFile)

# open deseq output file (input file here)
deseqFile = read.csv(paste0(inputFolder, deFile), header = TRUE)
head(deseqFile)
colnames(deseqFile)[1] = keyColumn

rc_annotated = left_join(deseqFile, annotations, by=keyColumn)
head(rc_annotated)

# Save annotated file
write.csv(rc_annotated, annotatedOutFileName, row.names = FALSE)
