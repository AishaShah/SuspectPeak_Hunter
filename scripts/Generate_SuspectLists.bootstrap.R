## HELP FUNCTION
print_help_message <- function() {
  cat("
    Generating common/suspect peaks list 
    Use the following command-line arguments:
    --raw_peaks_file [path] : Path to the CSV file containing raw peaks.
    --filt_peaks_file [path] : Path to the CSV file containing filtered peaks.
    --chr_lengths_file [path] : Path to the file containing chromosome lengths.
    --metadata [path] : Path for metadata having sample names, targtes, stages,
    --num_sample [integer] : Total number of samples (default 32).
    --percentage_threshold [number] : Percentage threshold for blacklisting (default 20).
    --min_shared_region_len [integer] : Minimum length of shared region (default 10).
    --output_BL_bases [path] : Output path for the bases blacklisted.
    --output_BL_bases_filt [path] : Output path for the filtered bases blacklisted.
    --ignore_filtered_plots  [TRUE or FALSE] : Generate one or two blacklists.
  ")
}
## CHECK COMMAND LINE ARGUMENTS
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || args[1] == "--help") {
  print_help_message()
  quit("no")
}

# Parse arguments based on their prefix
options <- list(raw_peaks_file = NULL, filt_peaks_file = NULL, chr_lengths_file = NULL,
                metadata = NULL, num_sample = 32, percentage_threshold = 20, 
                min_shared_region_len = 10, output_BL_bases = NULL, output_BL_bases_filt = NULL,ignore_filtered_plots=TRUE)

for (arg in args) {
  key <- sub("^--", "", strsplit(arg, "=")[[1]][1])
  value <- strsplit(arg, "=")[[1]][2]
  options[[key]] <- value
}

#if (is.null(options$raw_peaks_file) || is.null(options$filt_peaks_file) || 
#    is.null(options$chr_lengths_file) || is.null(options$output_BL_bases) || 
#    is.null(options$output_BL_bases_filt)) {
#  stop("Missing necessary command line arguments. Run with '--help' for more information.")
#}

ignore_filtered_plots <- options$ignore_filtered_plots

if (is.null(options$raw_peaks_file)  || 
    is.null(options$chr_lengths_file) || is.null(options$output_BL_bases)) {
  stop("Missing necessary command line arguments. Run with '--help' for more information.")
}




options(repos = c(CRAN = "https://cloud.r-project.org"))

# Function to check if a package is installed
is_package_installed <- function(package) {is.element(package, installed.packages()[,1])}
# Check if required packages are installed, and install them if necessary
required_packages <- c("reshape2", "ggplot2","tidyr","dplyr","patchwork","ggh4x","ggVennDiagram","gridExtra")
missing_packages <- required_packages[!sapply(required_packages, is_package_installed)]
print(missing_packages)
if (length(missing_packages) > 0) {install.packages(missing_packages)}
library(dplyr)
library(tidyr)
library(reshape2) # for melt function
#library(googlesheets4)
library(ggplot2)
library(patchwork)
library(ggh4x)
library(ggVennDiagram)
library(gridExtra)



process_file <- function(file) {
  f <- read.csv(file, sep="\t", header=T, check.names=FALSE)
  colnames(f) <- sub(".*\\/", "", colnames(f))
  colnames(f) <- gsub("\\.stringent.bed|'", "", colnames(f))
  f <- f %>% mutate(length = end - start)
  #f <- f %>% filter(num >= 2, length > 10)
  return(f)
}

## LOAD DATA
raw.peaks.file <- process_file(options$raw_peaks_file)
if (ignore_filtered_plots!="TRUE") {filtered.peaks.file <- process_file(options$filt_peaks_file)}
chr_lengths <- read.table(options$chr_lengths_file)
metadata <-  read.csv(options$metadata, sep = "\t", header = TRUE, check.names = FALSE)

## load metadata

# Define the order of cell stages
stage_order <- c("Dome", "Shield", "15-18s", "6hpf", "24hpf", "high")
# Convert the Stage column to a factor with custom order
metadata$Stage <- factor(metadata$Stage, levels = stage_order)
# stage colors
stage_colors <- list(Stage = c(Dome = "#59b59a", Shield = "#df8948", "15-18s" = "#9895c4", "6hpf" = "#e863a7", "24hpf" = "#8ebb5c", high = "#e8be48"))


## chr info
chr_lengths <- chr_lengths %>%
  select(1:2) %>%  # Select only the first two columns
  rename(chrom = 1, Total_Bases = 2)  # Rename the columns

## Functions

# Function to process raw peaks data and summarize peaks shared across sample, targets and cell stages
# Usage: summarize_peak_data(raw_peaks_file, metadata_file)
#   raw_peaks_file: the raw peaks data file
#   metadata_file: the metadata file containing sample information
summarize_peak_data <- function(raw_peaks, metadata) {
  # Pivot data to long format, filter, join with metadata, and perform additional operations
  processed_data <- raw_peaks %>%
    pivot_longer(cols = -c(chrom, start, end, num, list, length), names_to = "SampleID", values_to = "Region_Found") %>% 
    filter(Region_Found == 1) %>% 
    left_join(metadata %>% select(SampleID = sample, Targets = Target, Stage = Stage, experiment_group = ctrl, Target_Group), by = "SampleID") %>% 
    mutate(region = paste0(chrom, ":", start, "-", end)) %>%
    group_by(region, chrom, start, end) %>% 
    summarize(num_samples = n_distinct(SampleID),
              num_targets = n_distinct(Targets), 
              num_stages = n_distinct(Stage))
  
  return(processed_data)
}
# Example usage:
# summarized_data <- summarize_peak_data(raw_peaks, metadata.csv)


get_shared_peaks_plot <- function(data, plot_column = "num_targets", x_title = "", title = "") {
  plot <- data %>% 
    group_by(.data[[plot_column]]) %>%
    summarize(freq = n()) %>%
    ungroup() %>% 
    mutate(total_obs = sum(freq),  # Calculate the total number of observations
           percentage = freq / total_obs * 100,  # Calculate the percentage
           Peak_Type = ifelse(.data[[plot_column]] == 1, "Unique", "Shared")) %>%
    ggplot(aes(x = .data[[plot_column]], y = freq, fill = Peak_Type)) + 
    geom_col() + 
    geom_text(aes(label = paste0(round(percentage, 3), " %")), vjust = 0.5, hjust = -0.03, angle = 90, size = 2.0) +  # Add text annotations
    scale_fill_manual(values = c("Unique" = "darkgreen", "Shared" = "#585858")) +
    labs(x = x_title, y = "Number of\nRegions Shared", title = title) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))  # Adjusting y-axis margins
  return(plot)
}


analyze_blacklisted_bases <- function(chr_lengths, 
                                      summary_data, 
                                      num_sample = 32, 
                                      min_shared_region_len = 10, 
                                      ignore_chrom = NULL) {
  
  summary_data <- summary_data %>% mutate(length = end - start)
  results <- data.frame()
  
  for (percentage_threshold in seq(0, 100, by = 10)) {
    min_num_samples <- as.integer((num_sample / 100) * percentage_threshold)
    BL_bases <- summary_data %>% filter(num_samples >= min_num_samples,
                                        length >= min_shared_region_len)
    
    BL_bases <- BL_bases %>% group_by(chrom) %>% summarize(bases_blacklisted = sum(length))
    BL_bases <- merge(chr_lengths, BL_bases, by = "chrom", all = TRUE)
    BL_bases <- BL_bases %>% mutate(bases_blacklisted = ifelse(is.na(bases_blacklisted), 0, bases_blacklisted),
                                    Total_Bases = ifelse(is.na(Total_Bases), 0, Total_Bases))
    BL_bases[BL_bases$chrom == "MT", ]$bases_blacklisted <- BL_bases[BL_bases$chrom == "MT", ]$Total_Bases
    BL_bases <- BL_bases %>% mutate(Sequence1 = ifelse(substr(chrom, 1, 1) == "K", "scaffolds", "Chromosome"),
                                    Sequence = ifelse(substr(chrom, 1, 1) == "K", "scaffolds", chrom))
    
    BL_bases$Sequence <- factor(BL_bases$Sequence, levels = c(1:25, "MT", "scaffolds"), order = TRUE)
    
    if (!is.null(ignore_chrom)) {
      BL_bases <- BL_bases %>% filter(!(chrom %in% ignore_chrom))}
    
    percentage <- BL_bases %>% 
      summarize(bases_blacklisted = sum(bases_blacklisted),
                Total_Bases       = sum(Total_Bases)) %>%
      mutate(   Percentage        = ((bases_blacklisted / Total_Bases) * 100))
    
    results <- rbind(results, data.frame(min_num_samples = min_num_samples, 
                                         percentage_threshold = percentage_threshold, 
                                         Percentage = percentage$Percentage))
  }
  
  # Plotting the percentage of samples against % of blacklisted genome
  plot <- ggplot(results, aes(x = Percentage, y = percentage_threshold)) +
    geom_point() + geom_line() + 
    geom_text(aes(label = paste0("n=", min_num_samples)), 
              size = 2.1, color = "black", 
              vjust = 0, hjust = -0.5) +
    labs(y = "Percentage of Samples", x = "% of Bases Blacklisted") +
    theme_grey()
  # Plotting the number of samples against % of blacklisted genome
  plot <- ggplot(results, aes(x = Percentage, y = min_num_samples)) +
    geom_point() + geom_line() + 
    labs(y = "Number of Samples", x = "% of Refrence Genome Blacklisted") +
    theme_grey()
  
  #return(list(plot=plot,results=results))
  return(plot)
}



## raw:
summary.raw<-summarize_peak_data(raw.peaks.file,metadata) 
c1p1 <- get_shared_peaks_plot(summary.raw,plot_column = "num_targets",x_title = "Number of Targets", title = "Peak Calling on\nunfiltered bams")
c1p2 <-  get_shared_peaks_plot(summary.raw,plot_column = "num_samples",x_title = "Number of Samples", title = "Peak Calling on\nunfiltered bams")
c1p3 <- get_shared_peaks_plot(summary.raw,plot_column = "num_stages",x_title = "Number of Stages", title = "Peak Calling on\nunfiltered bams")
column1 <- (c1p1 / 
              (c1p2+theme(plot.title=element_blank())) / 
              (c1p3+theme(plot.title=element_blank()))   
) #+ plot_layout(guides = "collect") & theme(legend.position = "bottom")



## filtered
if (ignore_filtered_plots!="TRUE") {
summary.filt<-summarize_peak_data(filtered.peaks.file,metadata) 
c2p1 <- get_shared_peaks_plot(summary.filt,plot_column = "num_targets",x_title = "Number of Targets", title = "Peak Calling on\nfiltered bams")
c2p2 <-  get_shared_peaks_plot(summary.filt,plot_column = "num_samples",x_title = "Number of Samples", title = "Peak Calling on\nfiltered bams")
c2p3 <- get_shared_peaks_plot(summary.filt,plot_column = "num_stages",x_title = "Number of Stages", title = "Peak Calling on\nfiltered bams")
column2 <- ((c2p1 +theme(axis.title.y = element_blank()))/ 
              (c2p2+theme(plot.title=element_blank(), axis.title.y = element_blank())) / 
              (c2p3+theme(plot.title=element_blank(), axis.title.y = element_blank()))   
) #+ plot_layout(guides = "collect") & theme(legend.position = "bottom")

### plot1
FigureA <- (column1 | column2) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
} else {
  FigureA<- column1  + plot_layout(guides = "collect") & theme(legend.position = "bottom", text=element_text(size=5)) 
}

### plot2
## chossing threshold for % of genome blacklisted:
#It is better to check out of total regions called as peaks, how many are in suspect list?? per sample?? per target?? or overall?
if (ignore_filtered_plots!="TRUE") {
FigureB<- (analyze_blacklisted_bases(chr_lengths, 
                                     summary.raw, 
                                     num_sample = as.integer(options$num_sample),
                                     min_shared_region_len=as.integer(options$min_shared_region_len))) | 
  (  analyze_blacklisted_bases(chr_lengths, 
                               summary.filt,
                               num_sample = as.integer(options$num_sample),
                               min_shared_region_len=as.integer(options$min_shared_region_len)))

} else {
  FigureB<- (analyze_blacklisted_bases(chr_lengths, 
                                       summary.raw, 
                                       num_sample = as.integer(options$num_sample),
                                       min_shared_region_len=as.integer(options$min_shared_region_len))) + theme(legend.position = "bottom", text=element_text(size=5)) 
}
#output_file <- "BL"
if (ignore_filtered_plots!="TRUE") {
  # Save Figure A as JPEG
  jpeg(paste(options$output_BL_bases_filt, "_figureA_Sharing.jpg", sep=""), width = 7 * 300, height = 7.5 * 300, units = "px", res = 300, quality = 100)
  print(FigureA)    
  dev.off()
  
  # Save Figure B with different dimensions as JPEG
  jpeg(paste(options$output_BL_bases_filt, "_figureB_prcntBL.jpg", sep=""), width = 6 * 300, height = 3 * 300, units = "px", res = 300, quality = 100)
  print(FigureB) 
  dev.off()
} else {
# Save Figure A as JPEG
jpeg(paste(options$output_BL_bases, "_figureA_Sharing.jpg", sep=""), width = (7 * 300)/3, height = (7.5 * 300)/1.5, units = "px", res = 300, quality = 100)
print(FigureA)  
dev.off()

# Save Figure B with different dimensions as JPEG
jpeg(paste(options$output_BL_bases, "_figureB_prcntBL.jpg", sep=""), width = (6 * 300)/3, height = (3 * 300)/1.5, units = "px", res = 300, quality = 100)
print(FigureB)  
dev.off()
}

## generating blacklists
summary.raw <- summary.raw %>% mutate(length = end - start)
min_num_samples <- as.integer((as.integer(options$num_sample) / 100) * as.integer(options$percentage_threshold))

BL_bases <- summary.raw %>% ungroup %>%
  filter(num_samples >= min_num_samples,
         length >= as.integer(options$min_shared_region_len)) %>% 
  mutate( INFO=paste0("nsamples:ntargets:nstages|",num_samples,":",num_targets,":",num_stages))  %>%
  select(chrom,start,end,INFO)


if (ignore_filtered_plots!="TRUE") {
  summary.filt <- summary.filt %>% mutate(length = end - start)
  BL_bases.filt <- summary.filt %>% ungroup %>%
    filter(num_samples >= min_num_samples,
           length >= as.integer(options$min_shared_region_len))  %>% 
    mutate( INFO=paste0("nsamples:ntargets:nstages|",num_samples,":",num_targets,":",num_stages))  %>%
    select(chrom,start,end,INFO)
}

## save the files and overlape using bedtools

if (ignore_filtered_plots!="TRUE") {
write.table(BL_bases.filt %>% ungroup(), #%>% select(chrom, start,end,INFO),
            file=paste(options$output_BL_bases_filt, ".bed", sep=""),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}

write.table(BL_bases %>% ungroup(),# %>% select(chrom, start,end,INFO), 
            file=paste(options$output_BL_bases,".bed", sep=""),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

