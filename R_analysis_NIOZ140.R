##%######################################################%##
#                                                          #
####              Exploratory Analysis of               ####
####            16S Amplicon Sequencing Data            ####
#                                                          #
##%######################################################%##

##%######################################################%##
#                                                          #
####           Project: foils statia NIOZ 140           ####
#                                                          #
##%######################################################%##




################Packages_init###################
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
library(ggpubr)
library("wesanderson")
library("plotrix")
library("FactoMineR")
library("factoextra")
library(usedist)
library("heatmaply")


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

##########getting rid of wonky taxonomy assignments ###########
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

#loops to redefine weird taxonomy to a single common character string "unassigned" 
taxo <- as.data.frame(physeq_object@tax_table)
for (i in 1:nrow(taxo)) {
  for (y in 1:ncol(taxo)) {
    if 
    (any(str_detect(taxo[i,y], c("uncultured","Uncultured","metagenome", "Metagenome","unknown", "Unknown","NA")))) {
      taxo[i,y] <- "unassigned" }
  }
}

taxo <- tax_table(as.matrix(taxo)) #re-define as tax table object

physeq_object <- merge_phyloseq(physeq_object@otu_table, taxo, map) # merge updated taxonomy

#euk <- subset_taxa(physeq_object, Domain=="Eukaryota")
#arch <- subset_taxa(physeq_object, Domain=="Archaea")
#bact <- subset_taxa(physeq_object, Domain=="Bacteria")


############### alpha div ###################
summarize_phyloseq(physeq_object)
alpha_tab <-microbiome::alpha(physeq_object, index = "all")
write.csv(alpha_tab, file = "./alpha_div/alpha_div_indexes_microbiome_2.csv")
metad <- data.frame(physeq_object@sam_data) 
metad$Shannon <- alpha_tab$diversity_shannon 
metad$evenness_simpson <- alpha_tab$evenness_simpson 
m <- subset_samples(physeq_object, timepoint %in% c("T1", "T6"))

p <- ggboxplot(metad, x = "Material", y = "evenness_simpson",
               color = "Material", palette =c("#5FB233FF" ,"#6A7F93FF" ,"#F57206FF" ,"#EB0F13FF", "#8F2F8BFF", "#1396DBFF"),
               add = "jitter", shape = "treatment", size = 1) + facet_wrap(~timepoint)
p


a <- phyloseq::estimate_richness(physeq_object)
write.csv(a, file = "./alpha_div/alpha_div_indexes_phyloseq_2.csv")
plot <- plot_richness(m, "Material", "treatment", measures="Chao1")+facet_grid(treatment~timepoint)
plot + geom_boxplot(data=plot$data, aes(Material,value,color=NULL), alpha=0.3)+ labs(title = "Alpha Diversity", subtitle ="Chao1", x =NULL , y = NULL )+theme_light()

microbiome::plot_taxa_prevalence(physeq_object, "Phylum")+ theme(legend.position = "none") #prevalence


#########merge samples per surface all replicates together##########

#getting the phyloseq object as tidy tibble 
tidy_psmelt <- function(physeq) {
  ### INSERT Initial variable and rank name checking and modding from `psmelt`
  # Get the OTU table with taxa as rows
  rankNames = rank_names(physeq, FALSE)
  sampleVars = sample_variables(physeq, FALSE) 
   otutab <- otu_table(physeq)
  if (!taxa_are_rows(otutab)) {
    otutab <- t(otutab)
  }
  # Convert the otu table to a tibble in tidy form
   tb <- otutab %>% 
     as("matrix") %>%
     tibble::as_tibble(rownames = "OTU") %>%
    tidyr::gather("Sample", "Abundance", -OTU)
  # Add the sample data if it exists
  if (!is.null(sampleVars)) {
    sam <- sample_data(physeq) %>%
      as("data.frame") %>% 
      tibble::as_tibble(rownames = "Sample")
    tb <- tb %>%
      dplyr::left_join(sam, by = "Sample")
  }
  # Add the tax table if it exists
  if (!is.null(rankNames)) {
    tax <- tax_table(physeq) %>%
      as("matrix") %>%
      tibble::as_tibble(rownames = "OTU")
    tb <- tb %>%
      dplyr::left_join(tax, by = "OTU")
  }
  tb %>%
    arrange(desc(Abundance))
  # Optional conversion to a data frame doesn't affect the speed/memory usage
  # %>% as.data.frame
}


tidy_physeq <- tidy_psmelt(physeq_object)
tidy_physeq$LinkerPrimerSequence <- NULL
tidy_physeq$ReversePrimerSequence <- NULL
tidy_physeq$InputFileName <- NULL
tidy_physeq$BarcodeSequence <- NULL
tidy_physeq$BarcodeSequence_1 <- NULL



t2 <- tidy_physeq  %>% group_by(Sample) %>% mutate(Sample_rel_abund = Abundance / sum(Abundance)) %>% #relative abundance of each otu per sample
  ungroup() %>%
  group_by(description) %>% mutate( rep_rel_abund = Sample_rel_abund / sum(Sample_rel_abund)) %>% #relative abundance of each otu per number of samples in replicates
  ungroup() %>% 
  #Domain_section
  group_by(Sample, Domain) %>% 
  mutate(Domain_rel_abund_Sample = sum(Sample_rel_abund)) %>%  #domain relative abundance per sample 
  ungroup() %>% 
  group_by(description, Domain) %>% 
  mutate(Domain_st_dev_abund_samples = sd(Domain_rel_abund_Sample)) %>% # standard dev of domain relative abundances between replicates of description (ployner_timepoint_treatment)
  mutate(Domain_rep_rel_abund = sum(rep_rel_abund)) %>% #domain relative abundance per samples of desc 
  ungroup() %>%
  #Phylum_section
  group_by(Sample, Phylum) %>% 
  mutate(Phylum_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(description, Phylum) %>% 
  mutate(st_dev_Phylum_abund = sd(Phylum_rel_abund_Sample)) %>%
  mutate(Phyla_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Class_section
  group_by(Sample, Class) %>% 
  mutate(Class_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(description, Class) %>% 
  mutate(st_dev_Class_abund = sd(Class_rel_abund_Sample)) %>%
  mutate(Class_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Order_section
  group_by(Sample, Order) %>% 
  mutate(Order_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(description, Order) %>% 
  mutate(st_dev_Order_abund = sd(Order_rel_abund_Sample)) %>%
  mutate(Order_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Family_section
  group_by(Sample, Family) %>% 
  mutate(Family_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(description, Family) %>% 
  mutate(st_dev_Family_abund = sd(Family_rel_abund_Sample)) %>%
  mutate(Family_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Genus_section
  group_by(Sample, Genus) %>% 
  mutate(Genus_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(description, Genus) %>% 
  mutate(st_dev_Genus_abund = sd(Genus_rel_abund_Sample)) %>%
  mutate(Genus_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Species_section
  group_by(Sample, Species) %>% 
  mutate(Species_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(description, Species) %>% 
  mutate(st_dev_Species_abund = sd(Species_rel_abund_Sample)) %>%
  mutate(Species_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup()  

  
write_csv(t2, "./data/tidy_data_NIOZ140.csv")  
  

##############tables for each rank#############

Domain <- t2  %>%  select(timepoint,treatment, Material, description, Domain, Domain_rep_rel_abund,Domain_st_dev_abund_samples)%>% 
  distinct() 
write_csv(Domain, "./data/domain.csv")

Domain <- Domain %>% filter ( timepoint %in% c("T1", "T6")) #%>% mutate(timepoint = ifelse(Material =="negative_c", "T1", timepoint))

#colors_vector_to_personalize#
CPCOLS <- c("#199442", "#ED1F1F", "#F5EE2C", "#B636D6", "#3D68E0", "#EBA53D", "#00688B", "#00EE76", "#CD9B9B", "#00BFFF", "#FFF68F", "#FF7F50", "#68228B", "#ADFF2F", "#CD0000", "#0000FF", "#CD9B1D", "#FF34B3", "#BBFFFF", "#191970") 

ggplot(Domain, aes(x=Material, y= Domain_rep_rel_abund, fill=Domain))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = CPCOLS)+
  theme_classic2()+  facet_grid (timepoint~treatment)


#creating table phyla
Phyla <- t2 %>% select(treatment, timepoint, Material,description, Phylum, Phyla_rep_rel_abund ,st_dev_Phylum_abund)%>% 
            distinct() 
head(Phyla)
write_csv(Phyla, "./data/phyla.csv")

Phyla <- Phyla %>% mutate(Phylum = ifelse(Phyla_rep_rel_abund<0.01, "others<0.01", Phylum)) %>%
        filter ( timepoint %in% c("T1", "T6")) #%>% mutate(timepoint = ifelse(Material =="negative_c", "T1", timepoint))


ggplot(Phyla, aes(x=Material, y= Phyla_rep_rel_abund, fill=Phylum))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = CPCOLS)+
  theme_classic2()+  facet_grid (timepoint~treatment)


Class <- t2 %>% select(treatment, timepoint, Material, description, Class, Class_rep_rel_abund, st_dev_Class_abund )%>% 
          distinct()  
head(Class)
write_csv(Class, "./data/class.csv")

Class <- Class %>% mutate(Class = ifelse(Class_rep_rel_abund<0.01, "others<0.01", Class)) %>%
  filter ( timepoint %in% c("T1", "T6"))

ggplot(Class, aes(x=Material, y= Class_rep_rel_abund, fill=Class))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = CPCOLS)+
  theme_classic()+  facet_grid (timepoint~treatment)

Order <- t2 %>% select(treatment, timepoint, Material, description, Order, Order_rep_rel_abund, st_dev_Order_abund )%>% 
          distinct() #
head(Order)
write_csv(Order, "./data/order.csv")

Order <- Order %>% mutate(Order = ifelse(Order_rep_rel_abund<0.01, "others<0.01", Order)) %>% 
  filter ( timepoint %in% c("T1", "T6"))


ggplot(Order, aes(x=Material, y= Order_rep_rel_abund, fill=Order))+
  geom_bar(stat="identity", position="stack")+ scale_color_brewer()+ theme_classic()+  facet_grid (timepoint~treatment)

Family <- t2 %>% select(treatment, timepoint, Material, description, Family, Family_rep_rel_abund, st_dev_Family_abund )%>% 
  distinct()
head(Family)
write_csv(Family, "./data/family.csv")

Family <- Family %>%  mutate(Family = ifelse(Family_rep_rel_abund<0.01, "others<0.01", Family)) %>% 
  filter ( timepoint %in% c("T1", "T6"))

ggplot(Family, aes(x=Material, y= Family_rep_rel_abund, fill=Family))+
  geom_bar(stat="identity", position="stack")+ scale_color_brewer()+
  theme_classic()+  facet_grid (timepoint~treatment)

Genus <- t2 %>% select(treatment, timepoint, Material, description, Genus, Genus_rep_rel_abund, st_dev_Genus_abund )%>%
  distinct()
head(Genus)
write_csv(Genus, "./data/genus.csv")

Genus <- Genus %>% mutate(Genus = ifelse(Genus_rep_rel_abund<0.01, "others<0.01", Genus))  %>% 
  filter ( timepoint %in% c("T1", "T6"))

ggplot(Genus, aes(x=Material, y= Genus_rep_rel_abund, fill=Genus))+
  geom_bar(stat="identity", position="stack")+ scale_color_brewer()+ 
  theme_classic()+  facet_grid (timepoint~treatment)

Species <- t2 %>% select(treatment, timepoint, Material, description, Species, Species_rep_rel_abund, st_dev_Species_abund)%>% 
  distinct()
  
head(Species)
write_csv(Species, "./data/species.csv") 
  
Species <- Species %>%  mutate(Species = ifelse(Species_rep_rel_abund<0.01, "others<0.01", Species))  %>% 
  filter ( timepoint %in% c("T1", "T6"))

ggplot(Species, aes(x=Material, y= Species_rep_rel_abund, fill=Species))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = CPCOLS)+
  theme_classic()+  facet_grid (timepoint~treatment)+
  theme (axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold") ) 


################## heatmap on genus & no NAs #############
t2_no_na_genus <- filter(t2, !Genus==" NA")%>% group_by(description, Genus) %>% mutate(Genus_rep_rel_abund = sum(rep_rel_abund)) %>% 
  select(description, Genus, Genus_rep_rel_abund, treatment, timepoint, Material)%>% distinct() %>% 
  arrange(desc(Genus_rep_rel_abund))

t2_no_na_genus <- filter(t2_no_na, Genus %in% (unique(t2_no_na$Genus)[1:20]))


pal <- wes_palette("Zissou1", 100, type = "continuous")
ggplot(t2_no_na_genus,aes(x=description,y=Genus,fill=Genus_rep_rel_abund))+
  geom_tile(colour="white",size=0.25)+ labs(x="",y="")+
  geom_text(aes(label = round(Genus_rep_rel_abund, 3)), colour = "Black" , size = 3)+ scale_fill_gradientn(colours =  pal)+
  theme (axis.text.x = element_text(face="bold", angle=90), axis.text.y = element_text(face="bold") ) 


###########beta_div####################

####ord_amp_vis####
library(ampvis2)


p_otu<-  as(otu_table(physeq_object), "matrix")
map2 <- as(sample_data(physeq_object), "data.frame")
taxa <- as(tax_table(physeq_object), "matrix")
colnames(taxa)[1] <- "Kingdom"
for (t in 1:length(map2[,"timepoint"])) {
if (is.na(map2[t,"timepoint"])==T) {
map2[t,"timepoint"] <- "pcr"  
}  
}

for (t in 1:length(map2[,"treatment"])) {
  if (is.na(map2[t,"treatment"])==T) {
    map2[t,"treatment"] <- "-"  
  }  
}
colnames(map2)
map2 <- rownames_to_column(map2)
colnames(map2)[1] <- "sample ID"
p_otu <- data.frame(p_otu)
p_otu <- rownames_to_column(p_otu)
colnames(p_otu)[1] <- "OTU"
taxa <- data.frame(taxa)
taxa <- rownames_to_column(taxa)
colnames(taxa)[1] <- "OTU"


data <-  amp_load(p_otu, metadata = map2, taxonomy = taxa, check.names=F)
data$abund 

bray_PCOA <- amp_ordinate(data , transform = "none", type = "PCOA", distmeasure= "bray", sample_label_size=4, sample_color_by = "treatment", sample_label_by= "Material",sample_shape_by ="timepoint", species_plot = F, detailed_output=T)

#venn
venn <- amp_venn(data = data, group_by = "treatment", cut_a = 0, cut_f = 0, detailed_output = T)
venn2 <- amp_venn(data = data, group_by = "timepoint", cut_f = 70, detailed_output = T)




bray_PCOA$plot+ labs(title = "PCoA", subtitle ="distance: Bray-Curtis")+theme_pubr()


############heatmap##############
transposed_otu <-  t(as(otu_table(physeq_object), "matrix"))
#bray_distance<-vegdist(as(otu_table(physeq_object), "matrix"), method="bray")
bray_distance2<-vegdist(transposed_otu, method="bray")

#devtools::install_github("kylebittinger/usedist")

sample_bray <- (as(bray_distance2, "matrix"))
row.names(sample_bray) <- map2$description
colnames(sample_bray) <- map2$surface
      
pal <- wes_palette("Zissou1",70, type = "continuous")
          
heat <- heatmap(sample_bray, scale = "none", col = pal)


heatmaply::heatmaply(sample_bray, row_text_angle = 0,
                     column_text_angle = 90)



###############################################ASV#############################################################################

asvTable<-read_tsv("../foils_statia_2019/runs/dada2/asv/taxonomy_dada2/asvTable_noSingletons_no_1st_line.txt")
rows <- asvTable$`#OTU ID`
tax <- asvTable$taxonomy
tax <- read_csv2(tax, col_names = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
tax <- as.data.frame(tax,  row.names = rows) %>% as.matrix() %>% tax_table()
asv
otu <- as.matrix(read.delim("./data/otuTable_noSingletons_no_first_line.txt", row.names = 1))
otu <- otu_table(otu, taxa_are_rows = T)
map <- sample_data(read.delim("./data/mapping_file_details.txt", row.names = 1, na.strings = c("NA", "")))
physeq_object = merge_phyloseq(otu, tax, map)                 



####basic_info_asv#####
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

#####subset T3 & merge asv####
sub1 <- subset_samples(physeq_object, timepoint %in% c("T1", "T6"))
sub2 <- subset_samples(physeq_object, surface=="negative_c")
physeq_object <- merge_phyloseq(sub1, sub2)

########## wonky taxonomy assignments asv###########
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

taxo <- as.data.frame(physeq_object@tax_table)
for (i in 1:nrow(taxo)) {
  for (y in 1:ncol(taxo)) {
    if 
    (any(str_detect(taxo[i,y], c("uncultured","Uncultured","metagenome", "Metagenome","unknown", "Unknown","NA")))) {
      taxo[i,y] <- "unassigned" }
  }
}

taxo <- tax_table(as.matrix(taxo))


physeq_object <- merge_phyloseq(physeq_object@otu_table, taxo, map)


#t1 <- subset_samples (physeq_object, timepoint=="T1")

#t6 <- subset_samples(physeq_object, timepoint=="T6")


#merge_samples(GlobalPatterns, group = factor(as.character(unlist(sample_data(GlobalPatterns)[,"SampleType"]))))
#euk <- subset_taxa(physeq_object, Domain=="Eukaryota")
#arch <- subset_taxa(physeq_object, Domain=="Archaea")
#bact <- subset_taxa(physeq_object, Domain=="Bacteria")




############### alpha div asv###################
summarize_phyloseq(physeq_object)
alpha_tab <-microbiome::alpha(physeq_object, index = "all")
write.csv(alpha_tab, file = "./alpha_div/alpha_div_indexes_microbiome_2.csv")
metad <- data.frame(physeq_object@sam_data) 
metad$Shannon <- alpha_tab$diversity_shannon 
metad$evenness_simpson <- alpha_tab$evenness_simpson 
m <- subset_samples(physeq_object, timepoint %in% c("T1", "T6"))

p <- ggboxplot(metad, x = "Material", y = "evenness_simpson",
               color = "Material", palette =c("#5FB233FF" ,"#6A7F93FF" ,"#F57206FF" ,"#EB0F13FF", "#8F2F8BFF", "#1396DBFF"),
               add = "jitter", shape = "treatment", size = 1) + facet_wrap(~timepoint)
p


a <- phyloseq::estimate_richness(physeq_object)
write.csv(a, file = "./alpha_div/alpha_div_indexes_phyloseq_2.csv")
plot <- plot_richness(m, "Material", "treatment", measures="Chao1")+facet_grid(treatment~timepoint)
plot + geom_boxplot(data=plot$data, aes(Material,value,color=NULL), alpha=0.3)+ labs(title = "Alpha Diversity", subtitle ="Chao1", x =NULL , y = NULL )+theme_light()

microbiome::plot_taxa_prevalence(physeq_object, "Phylum")+ theme(legend.position = "none") #prevalence
