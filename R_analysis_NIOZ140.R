library(phyloseq)
library(scales)
library(grid)
library(dplyr)
library(tidyverse)
library(vegan)
library(rmarkdown)
library(knitr)
library("BiocManager")
library("devtools")
library("DESeq2")
library("microbiome")
library("microbiomeutilities")
library("TreeSummarizedExperiment")


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
euk <- subset_taxa(physeq = physeq_object, Domain=="Eukaryota")
arch <- subset_taxa(physeq = pseq, Domain=="Archaea")


t3 <- transform_sample_counts(t3, function(x)  x/sum(x))
t3<- microbiome::transform(t3, "compositional")

CPCOLS <- c("#199442", "#ED1F1F", "#F5EE2C", "#B636D6", "#3D68E0", "#EBA53D", "#00688B", "#CDCD00", "#EE3A8C", "#00EE76", "#CD9B9B", "#00BFFF", "#FFF68F", "#FF7F50", "#68228B", "#ADFF2F", "#CD0000", "#0000FF", "#CD9B1D", "#FF34B3", "#BBFFFF", "#191970") 
C<- c("#14A821", "#E6DB45", "#EB2C2C", "#4BEE8", "#C66EE6")
pbar <- plot_bar(arch,x="Material", fill = "Family")
pbar+ geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + scale_fill_manual(values=CPCOLS) + 
  scale_color_manual(values= CPCOLS )+ theme_classic()+facet_wrap (.~timepoint)
######bar_plots_diversity#####
pseq <- microbiome::transform(physeq_object, "compositional")

p_bar <- plot_bar(pseq, fill = "Domain")
p_bar+geom_bar(aes(color=Domain, fill=Domain), stat="identity", position="stack") + theme_classic() + scale_fill_manual(values=wes_palette( name="Darjeeling1")) + 
  scale_color_manual(values=wes_palette( name="Darjeeling1"))+ theme_classic() +  facet_wrap (.~treatment) 

pbar <- plot_bar(merged_plastic, fill = "Phylum")
pbar+geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + theme_classic()


###############MIA & microbiome packages###################


summarize_phyloseq(physeq_object)
alpha_tab <-microbiome::alpha(physeq_object, index = "all")
write.csv(alpha_tab, file = "./alpha_div_indexes.csv")
install.packages("ggrepel")

################alpha_div_idexes################
a_div <- readxl::read_xlsx("alpha_div_indexes.xlsx", na = "NA")

plotGP <- plot_richness(physeq_object, "Material", "treatment", measures="Simpson")+facet_wrap()
plotGP + geom_boxplot(data=plotGP$data, aes(Material,value,color=NULL), alpha=0.1)


