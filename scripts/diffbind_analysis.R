#!/usr/bin/env Rscript


# Library having R packages <-- installed using script scripts/install_R_pkgs.R
user_lib <- "/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/R_libs"
.libPaths(user_lib)

# Load necessary library
library(DiffBind,lib.loc = user_lib)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_metadata <- args[1]        # First argument: metadata file (TSV) for peaks with SuspectList regions
input_metadata_noSL <- args[2]   # Second argument: metadata file (TSV) for peaks without SuspectList regions
output_pdf <- args[3]            # Third argument: output PDF file
output_pdf_noSL <- args[4]       # Fourth argument: output PDF file
# Read metadata 
metadata <- read.csv(input_metadata, sep = "\t")
metadata_noSL <- read.csv(input_metadata_noSL, sep = "\t")

# Perform DiffBind analysis
DBA <- dba(sampleSheet=metadata)

# Save plots to PDF
pdf(output_pdf)
dba.plotHeatmap(DBA, ColAttributes=c(DBA_FACTOR, DBA_REPLICATE, DBA_CONDITION))
dev.off()

# Process noSL metadata
DBA_noSL <- dba(sampleSheet=metadata_noSL)

pdf(output_pdf_noSL)
dba.plotHeatmap(DBA_noSL, ColAttributes=c(DBA_FACTOR, DBA_REPLICATE, DBA_CONDITION))
dev.off()