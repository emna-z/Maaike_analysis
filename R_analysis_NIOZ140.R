library(phyloseq)
library(ggplot2)
library(scales)
library(grid)
library(dplyr)
#####working_dir####
setwd(dir = "C:/Users/ezeghal/Documents/Maaike_analysis")
getwd()


#####data_import######
tax <- as.matrix(read.delim('taxTable_noSingletons.txt', row.names = 1, na.strings = "NA"))
tax <- tax_table(tax)
otu <- as.matrix(read.delim("otuTable_noSingletons_no_first_line.txt", row.names = 1))
otu <- otu_table(otu, taxa_are_rows = T)
map <- sample_data(read.delim('mapping_file_details.txt', row.names = 1))
physeq_object = merge_phyloseq(otu, tax, map)                 



####basic_info#####
ntaxa(physeq_object)
nsamples(physeq_object)  ###there's 71 samples in the OTU table and 74 in the mapping file any idea why? never mind, I found in the report that 3 of them had nothing whatsoever
sample_names(physeq_object)
taxa_names(physeq_object)
rank_names(physeq_object)
sample_sums(physeq_object)
taxa_sums(physeq_object)
min(sample_sums(physeq_object)) #it says 7 sequences here. In the report it says 8, not much of a difference but wondering how to deal wih a sample like that
max(sample_sums(physeq_object))
get_taxa_unique(physeq_object,"Phylum")
length(get_taxa_unique(physeq_object,"Phylum"))


#####subset & merge for plotting####
t1 <- subset_samples (physeq_object, timepoint=="T1")
t3 <- subset_samples (physeq_object, timepoint=="T3")
t6 <- subset_samples(physeq_object, timepoint=="T6")
#merge_samples(GlobalPatterns, group = factor(as.character(unlist(sample_data(GlobalPatterns)[,"SampleType"]))))

t3 <- transform_sample_counts(t3, function(x)  x/sum(x))
t3<- microbiome::transform(t3, "compositional")

pbar <- plot_bar(t3,x="Material", fill = "Domain")
pbar+ geom_bar(aes(color=Domain, fill=Domain), stat="identity", position="stack") + scale_fill_manual(values=wes_palette( name="Darjeeling1")) + 
  scale_color_manual(values=wes_palette( name="Darjeeling1"))+ theme_classic() +  facet_wrap (.~treatment)


######bar_plots_diversity#####
pseq <- microbiome::transform(physeq_object, "compositional")

p_bar <- plot_bar(pseq, fill = "Domain")
p_bar+geom_bar(aes(color=Domain, fill=Domain), stat="identity", position="stack") + theme_classic() + scale_fill_manual(values=wes_palette( name="Darjeeling1")) + 
  scale_color_manual(values=wes_palette( name="Darjeeling1"))+ theme_classic() +  facet_wrap (.~treatment) 

pbar <- plot_bar(merged_plastic, fill = "Phylum")
pbar+geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + theme_classic()


###############MIA & microbiome packages###################
library("tidyverse")
library("vegan")
library("phyloseq")
library("rmarkdown")
library("knitr")
library("BiocManager")
library("devtools")
library("DESeq2")
library("microbiome")
library("microbiomeutilities")
library("TreeSummarizedExperiment")

summarize_phyloseq(physeq_object)
alpha_tab <-microbiome::alpha(physeq_object, index = "all")
write.csv(alpha_tab, file = "./alpha_div_indexes.csv")
install.packages("ggrepel")
