#BiocManager::install("rtracklayer")
library(rtracklayer)
library(dplyr)

# Customize
inputFolder = "C:/nordman/deseq2/results_c1/"
deFile = "DEresults_comp1.csv"
annotations = readGFF("C:/nordman/genomefiles/dm6.gtf")

# Remove duplicate annotations
colnames(annotations)

# remove exon entries
annotations = annotations[annotations$type != "exon",] 

# remove duplicates
annotations = annotations[!duplicated(annotations$gene_id),]

# remove other rows
annotations = subset(annotations, select=-c(note, orig_transcript_id, product, 
  transcript_biotype, exon_number, orig_protein_id, transl_table, codons,
  partial, inference, exception, part, pseudo, transl_except, protein_id,
  source,	type, score, phase, transcript_id, gbkey))
head(annotations)

# Script
annotatedOutFileName = paste0(inputFolder, "Annotated_", deFile)

rc = read.csv(paste0(inputFolder, deFile), header = TRUE)
head(rc)
colnames(rc)[1] = "gene_id"
colnames(rc)

rc_annotated = left_join(rc, annotations, by="gene_id")
head(rc_annotated)

# Save annotated file
write.csv(rc_annotated, annotatedOutFileName, row.names = FALSE)
