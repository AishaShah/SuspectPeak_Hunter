
library(ggplot2)
library(ggVennDiagram)
library(dplyr)
library(tidyr)


args <- commandArgs(trailingOnly = TRUE)
overlap <- read.csv(args[1],sep = "\t", header = T,check.names = FALSE)
min_len <- args[2]
output_path <- args[3]

overlap <- overlap %>% mutate(length=end-start)
overlap.short<-overlap  %>% filter(length>min_len)

x<- list(blacklist1=overlap.short %>% mutate(region=paste0(chrom,":",start,"-",end)) %>% filter(blacklist1==1) %>% select(region) %>% unlist(),
         blacklist2=overlap.short %>% mutate(region=paste0(chrom,":",start,"-",end)) %>% filter(blacklist2==1) %>% select(region) %>% unlist())

plot <- ggVennDiagram(x, category.names = c("BL 1","BL 2")) + 
  scale_fill_distiller(palette = "Reds", direction = 1)+ 
  labs(title = "Overlap between Blacklist 1 and Blacklist 2")

ggsave(paste0(output_path,"/venn_diagram.png"), plot, width = 4, height = 4, units = "in", dpi = 300)
