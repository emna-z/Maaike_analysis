library(phyloseq)
library(grid)
library(tidyverse)
library(vegan)
library(rmarkdown)
library(knitr)
library("DESeq2")
library("microbiome")
library("microbiomeutilities")
library("TreeSummarizedExperiment")


#####data_import######
tax <- as.matrix(read.delim('./data/taxTable_noSingletons.txt', row.names = 1, na.strings = "NA"))
tax <- tax_table(tax)
otu <- as.matrix(read.delim("./data/otuTable_noSingletons_no_first_line.txt", row.names = 1))
otu <- otu_table(otu, taxa_are_rows = T)
map <- sample_data(read.delim("./data/mapping_file_details.txt", row.names = 1))
physeq_object = merge_phyloseq(otu, tax, map)                 



####basic_info#####
summarize_phyloseq(physeq_object)
ntaxa(physeq_object)
nsamples(physeq_object)  ###there's 71 samples in the OTU table and 74 in the mapping file any idea why? never mind, I found in the report that 3 of them had nothing whatsoever
sample_names(physeq_object)
taxa_names(physeq_object)
rank_names(physeq_object)
sample_sums(physeq_object)
taxa_sums(physeq_object)
min(sample_sums(physeq_object)) #it says 7 sequences here. In the report it says 8, not much of a difference but wondering how to deal wih a sample like that
physeq_object <-  subset_samples(physeq_object, sample_names(physeq_object)!="NIOZ140.90") #eliminate the sample with 7seq
#physeq_object <-  subset_samples(physeq_object,(sample_sums(physeq_object) >= 1000))
max(sample_sums(physeq_object))

##########gettingrid of wonky taxonomy assignments ###########
get_taxa_unique(physeq_object, "Domain") # unassigned in Domains
physeq_object <- subset_taxa(physeq_object, !is.na(Domain) & !Domain%in% c("", "Unassigned")) #let's eliminate those otus
get_taxa_unique(physeq_object, "Domain") # all good now
get_taxa_unique(physeq_object, "Phylum") # let's check the Phyla, there's "NA"
physeq_object <- subset_taxa(physeq_object, !is.na(Phylum) & !Phylum%in% c("NA"," NA" )) 
get_taxa_unique(physeq_object, "Phylum")
length(get_taxa_unique(physeq_object,"Phylum"))

#####subset & merge for plotting####
physeq_object <- subset_samples(physeq_object, timepoint %in% c("T1", "T6", "negative_c")) 
physeq_object <- prune_taxa(taxa_sums(physeq_object) > 1, physeq_object)

#t1 <- subset_samples (physeq_object, timepoint=="T1")

#t6 <- subset_samples(physeq_object, timepoint=="T6")


#merge_samples(GlobalPatterns, group = factor(as.character(unlist(sample_data(GlobalPatterns)[,"SampleType"]))))
euk <- subset_taxa(physeq_object, Domain=="Eukaryota")
arch <- subset_taxa(physeq_object, Domain=="Archaea")
bact <- subset_taxa(physeq_object, Domain=="Bacteria")

############### alpha div ###################
summarize_phyloseq(physeq_object)
alpha_tab <-microbiome::alpha(physeq_object, index = "all")
write.csv(alpha_tab, file = "./alpha_div/alpha_div_indexes_no_t3_no_8seq.csv")

microbiome::plot_taxa_prevalence(physeq_object, "Phylum", 1)+ theme(legend.position = "none") #prevalence

###### filter otus #########

#funtion to filter otu fraction < 0.01
filter_0.01 <- function(physeq, frac = 0.01){
  ## Estimate total abundance of OTUs
  total <- sum(phyloseq::taxa_sums(physeq))
  ## Remove OTUs
  res <- phyloseq::filter_taxa(physeq, function(x){ ( sum(x)/total ) > frac }, prune = TRUE)
  return(res)
}

physeq_object <- filter_0.01(physeq_object)

#condition <- function (x) {sum(x) >= 5}  
#taxaToKeep <- genefilter_sample(physeq_object, condition,  2 ) # min number of reads 5 & present in at least 2 samples
#physeq_object <- prune_taxa(taxaToKeep, physeq_object)



####### let's make bar plots #########
physeq_t <- transform_sample_counts(physeq_object, function(x)  x/sum(x))

p_bar <- plot_bar(physeq_object, x= "Material", fill = "Domain")
p_bar+geom_bar(aes(color=Domain, fill=Domain), stat="identity", position="stack") + theme_classic() +  facet_grid (treatment~timepoint) 


#colors#
#CPCOLS <- c("#199442", "#ED1F1F", "#F5EE2C", "#B636D6", "#3D68E0", "#EBA53D", "#00688B", "#CDCD00", "#EE3A8C", "#00EE76", "#CD9B9B", "#00BFFF", "#FFF68F", "#FF7F50", "#68228B", "#ADFF2F", "#CD0000", "#0000FF", "#CD9B1D", "#FF34B3", "#BBFFFF", "#191970") 
#("#14A821", "#E6DB45", "#EB2C2C", "#4BEE8", "#C66EE6")

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




################alpha_div_idexes################
a_div <- readxl::read_xlsx("alpha_div_indexes.xlsx", na = "NA")

plotGP <- plot_richness(physeq_object, "Material", "treatment", measures="Simpson")+facet_wrap()
plotGP + geom_boxplot(data=plotGP$data, aes(Material,value,color=NULL), alpha=0.1)


