# Install required packages if not already installed
# Set the library path
lib_path <- "/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/R_libs"
dir.create(lib_path, showWarnings = FALSE, recursive = TRUE)

# Install required packages if not already installed
#install_if_missing <- function(p) {
#  if (!requireNamespace(p, quietly = TRUE)) {
#    tryCatch(
#      install.packages(p, lib = lib_path, repos = "http://cran.us.r-project.org", dependencies = TRUE),
#      error = function(e) {
#        message(paste("Failed to install package:", p))
#        message(e)
#      }
#    )
#  }
#}

# List of required packages
#required_packages <- c("pheatmap", "dplyr", "gridExtra", "grid", "ggplot2")

#lapply(required_packages, install_if_missing)

# Load necessary libraries
library(pheatmap, lib.loc = lib_path)
library(dplyr, lib.loc = lib_path)
library(gridExtra, lib.loc = lib_path)
library(grid, lib.loc = lib_path)
library(ggplot2, lib.loc = lib_path)
library(farver, lib.loc = lib_path)


sample_array <- read.table("/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/snakemake/SuspectPeak_Hunter/array_rep1.wo_IgG.wo_NegCtrls.tsv", header = TRUE)

target_color_list <- list(Targets=c(
  hoxa13b = "#FF0000",
  FLAG = "#FFFFFF",
  rbpja = "#00FF00",
  H3K4me1 = "#0000FF",
  H3K27ac = "#FFA500",
  H3K27me3 = "#FFFF00",
  H2AK119ub = "#800080",
  H3K4me3 = "#FFC0CB",
  H3K9me3 = "#00FFFF",
  pol2 = "#FF00FF",
  IgG = "#A52A2A",
  Prdm1a = "#000000",
  Tfap2a = "#808080",
  Nanog_HA = "#006400",
  H2AZ = "#ADD8E6",
  Ser2P = "#90EE90",
  pou5f3flag_FLAG = "#FFFFE0",
  pou5f3flag_HA = "#D3D3D3",
  pou5f3mut_HA = "#8B0000")#,
  #BioProject=c(PRJNA734348  = "#1B9E77",
  #CnR_Dora     = "#D95F02",
  #PRJNA719435  ="#7570B3",
  #PRJNA738523  ="#E7298A",
  #CnR_edlyn    ="#66A61E",
  #PRJNA899096  ="#E6AB02",
  #PRJNA1054088="#A6761D")
  , Target_Group=c("active_mark"=alpha("darkgreen",0.3), "inactive_mark"=alpha("darkred",0.3), "TF"=alpha("darkblue",0.3))
)




# Define the directory containing the files
#dir_path <- "/home/saisha/Desktop/unil_work/Aisha-Dora/snakemake/SuspectPeak_Hunter/05.Generating_SuspectLists/00.DownSampled_3000000.mergepeaks"
dir_path <- "/work/FAC/FBM/CIG/nvastenh/competition_model/Aisha-Dora/snakemake/SuspectPeak_Hunter/05.bootstrapping/05.Generating_SuspectLists.DS_10000000.mergepeaks"
# Get the list of files ending in "mergedPeaks.binary.tab"
file_list <- list.files(path = dir_path, pattern = "mergedPeaks\\.binary\\.tab$", full.names = TRUE)

# Function to process each file
process_file <- function(file) {
  f <- read.csv(file, sep="\t", header=T, check.names=FALSE)
  colnames(f) <- sub(".*\\/", "", colnames(f))
  colnames(f) <- gsub("\\.stringent.bed|'", "", colnames(f))
  f <- f %>% mutate(length = end - start)
  f <- f %>% filter(num >= 5, length > 20)
  return(f)
}

# Initialize a list to store the heatmaps
heatmaps <- list()

# Loop through each file in file_list
for (file in file_list) {
  # Process the file
  print(paste0("Processing:",basename(file)))
  f1 <- process_file(file)

  # Extract the matrix data (assuming the relevant columns are 6 to 14)
  data.mat <- f1[, 6:14]

  # Set row names to be region names
  rownames(data.mat) <- f1 %>%
    mutate(region = paste0(chrom, ":", start, "-", end)) %>%
    select(region) %>%
    as.list() %>%
    unlist()

  # Transpose the data matrix
  data.mat <- t(data.mat)

  # Assuming you have a sample_array data frame with Target and Stage information
  # (This part may need to be adapted based on your actual sample_array structure)
  Targets <- data.frame(Targets = sample_array$Target,
                        cell_stage = sample_array$Stage,
                        Target_Group=sample_array$Target_Group)
  rownames(Targets) <- sample_array$sample
  Targets <- Targets[rownames(data.mat), , drop = FALSE]

  # Generate the heatmap and store it in the list
   print(paste0("Generate Heatmap:",basename(file)))
  heatmaps[[basename(file)]] <- pheatmap(data.mat,
                                         cluster_rows = TRUE,
                                         cluster_cols = TRUE,
                                         annotation_colors = target_color_list,
                                         main = paste("Peaks regions shared in",gsub("\\.mergedPeaks\\.binary\\.tab", "", basename(file))),
                                         annotation_row = Targets,
                                         color = c("grey", "#27292b"), show_colnames = FALSE, treeheight_col = 0, legend = FALSE, silent = FALSE)
}

# Create the directory if it doesn't exist

dir.create("bootstrapping_pheatmaps.no_THR", showWarnings = FALSE)

# Loop through each heatmap and save as PNG
for (i in seq_along(heatmaps)) {
  png(file = paste0("bootstrapping_pheatmaps.no_THR/", names(heatmaps)[i], ".png"), width = 500, height = 350)
  grid.newpage()
  grid.draw(heatmaps[[i]]$gtable)
  dev.off()
}


