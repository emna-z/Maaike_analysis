library(devtools)
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
#install.packages("ggpubr")
library(ggpubr)
library("wesanderson")

#####data_import######
tax <- as.matrix(read.delim('./data/taxTable_noSingletons.txt', row.names = 1, na.strings = "NA"))
tax <- tax_table(tax)
otu <- as.matrix(read.delim("./data/otuTable_noSingletons_no_first_line.txt", row.names = 1))
otu <- otu_table(otu, taxa_are_rows = T)
map <- sample_data(read.delim("./data/mapping_file_details.txt", row.names = 1, na.strings = c("NA", "")))
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

#####subset T3 & merge ####
sub1 <- subset_samples(physeq_object, timepoint %in% c("T1", "T6"))
sub2 <- subset_samples(physeq_object, surface=="negative_c")
physeq_object <- merge_phyloseq(sub1, sub2)

##########gettingrid of wonky taxonomy assignments ###########
get_taxa_unique(physeq_object, "Domain") # unassigned in Domains
physeq_object <- subset_taxa(physeq_object, !is.na(Domain) & !Domain%in% c("", "Unassigned")) #let's eliminate those otus
get_taxa_unique(physeq_object, "Domain") # all good now
get_taxa_unique(physeq_object, "Phylum") # let's check the Phyla, there's "NA"
physeq_object <- subset_taxa(physeq_object, !is.na(Phylum) & !Phylum%in% c("NA"," NA" )) 
get_taxa_unique(physeq_object, "Phylum")
length(get_taxa_unique(physeq_object,"Phylum"))
physeq_object <- subset_taxa(physeq_object, !Order%in% c(" Chloroplast")) 
physeq_object <- subset_taxa(physeq_object, !Family%in% c(" Mitochondria"))

physeq_object <- prune_taxa(taxa_sums(physeq_object) > 1, physeq_object) #no singletons

#t1 <- subset_samples (physeq_object, timepoint=="T1")

#t6 <- subset_samples(physeq_object, timepoint=="T6")


#merge_samples(GlobalPatterns, group = factor(as.character(unlist(sample_data(GlobalPatterns)[,"SampleType"]))))
#euk <- subset_taxa(physeq_object, Domain=="Eukaryota")
#arch <- subset_taxa(physeq_object, Domain=="Archaea")
#bact <- subset_taxa(physeq_object, Domain=="Bacteria")

############### alpha div ###################
summarize_phyloseq(physeq_object)
alpha_tab <-microbiome::alpha(physeq_object, index = "all")
write.csv(alpha_tab, file = "./alpha_div/alpha_div_indexes_no_t3_no_8seq.csv")
metad <- data.frame(physeq_object@sam_data) 
metad$Shannon <- alpha_tab$diversity_shannon 
metad$evenness_simpson <- alpha_tab$evenness_simpson 
metad <- filter(metad, timepoint %in% c("T1", "T6"))

p <- ggboxplot(metad, x = "Material", y = "evenness_simpson",
               color = "Material", palette =c("#5FB233FF" ,"#6A7F93FF" ,"#F57206FF" ,"#EB0F13FF", "#8F2F8BFF", "#1396DBFF"),
               add = "jitter", shape = "treatment", size = 1) + facet_wrap(~timepoint)
p


plot <- plot_richness(physeq_object, "Material", "treatment", measures="Simpson")+facet_wrap(~timepoint)
plot + geom_boxplot(data=plotGP$data, aes(Material,value,color=NULL), alpha=0.3)

microbiome::plot_taxa_prevalence(physeq_object, "Phylum")+ theme(legend.position = "none") #prevalence


#########merge samples per surface all replicates together##########
merged <- collapse_replicates(physeq_object, method = "sample", replicate_fields = c("description", "surface"))




####### let's make bar plots #########
#mean_PGroup = sapply(levels((as.factor(physeq_object_f@sam_data$description))),function(i){
 # rowMeans(otu_table(PGroup)[,SampleType==i])})
#install.packages("wesanderson")

physeq_t <- transform_sample_counts(merged, function(x)  x/sum(x)) #get relative abundance
physeq_t_no_control <- subset_samples(physeq_t, timepoint %in% c("T1", "T6")) #get ridof negative control for esthetic reasons :)

p_bar <- plot_bar(physeq_t_no_control, x= "Material", fill = "Domain")
p_bar+geom_bar(aes(color=Domain, fill=Domain), stat="identity", position="stack") + theme_classic() + scale_fill_manual(values=wes_palette( name="Darjeeling1")) + 
  scale_color_manual(values=wes_palette( name="Darjeeling1"))+ theme_classic() +  facet_grid (timepoint~treatment) 

p_bar <- plot_bar(physeq_t_no_control, x= "Material", fill = "Phylum")
p_bar+geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + theme_pubclean() +  facet_grid (treatment~timepoint) 


#colors#
CPCOLS <- c("#199442", "#ED1F1F", "#F5EE2C", "#B636D6", "#3D68E0", "#EBA53D", "#00688B", "#CDCD00", "#EE3A8C", "#00EE76", "#CD9B9B", "#00BFFF", "#FFF68F", "#FF7F50", "#68228B", "#ADFF2F", "#CD0000", "#0000FF", "#CD9B1D", "#FF34B3", "#BBFFFF", "#191970") 
c("#14A821", "#E6DB45", "#EB2C2C", "#4BEE8", "#C66EE6")

###### filter otus #########


#condition <- function (x) {sum(x) >= 5}  
#taxaToKeep <- genefilter_sample(physeq_object, condition,  2 ) # min number of reads 5 & present in at least 2 samples
#physeq_object_f <- prune_taxa(taxaToKeep, physeq_object)


#funtion to filter otu fraction < 0.01
filter_0.001 <- function(physeq, frac = 0.001){
  ## Estimate total abundance of OTUs
  total <- sum(phyloseq::taxa_sums(physeq))
  ## Remove OTUs
  res <- phyloseq::filter_taxa(physeq, function(x){ ( sum(x)/total ) > frac }, prune = TRUE)
  return(res)
}

physeq_object_f <- filter_0.001(physeq_object)

summarize_phyloseq(physeq_object_f)

merged <- collapse_replicates(physeq_object_f, method = "sample", replicate_fields = c("description", "surface"))
merged<- transform_sample_counts(merged, function(x)  x/sum(x))


pbar <- plot_bar(merged,x="surface", fill = "Order")
pbar+ geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack") + scale_fill_manual(values=CPCOLS) + 
  scale_color_manual(values= CPCOLS )+ theme_classic()+facet_wrap (.~timepoint)



