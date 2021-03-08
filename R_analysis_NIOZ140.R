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


tidy_psmelt <- function(physeq) {
  ### INSERT Initial variable and rank name checking and modding from `psmelt` source HERE
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

#merged <- collapse_replicates(physeq_object, method = "sample", replicate_fields = c("description", "surface"))
tidy_physeq <- tidy_psmelt(physeq_object)

t2 <- tidy_physeq  %>% group_by(Sample) %>% mutate(Sample_relative_abundance = Abundance / sum(Abundance))
t2 <- t2  %>% group_by(description) %>% mutate( relative_abundance = Sample_relative_abundance / sum(Sample_relative_abundance))

t2$LinkerPrimerSequence <- NULL
t2$ReversePrimerSequence <- NULL
t2$InputFileName <- NULL
t2$BarcodeSequence <- NULL
t2$BarcodeSequence_1 <- NULL



#physeq_object = merge_phyloseq(otu, tax_table(physeq_object), sample_data(physeq_object))    
####### let's make bar plots #########
#t3 <- filter(t2, relative_abundance>=0.01) %>% ungroup()
#others <- filter(t2, relative_abundance<0.01) %>% ungroup()



#t4 <- full_join(t3,others)
#t5 <- filter (t4, timepoint %in% c("T1", "T6")) 

Phyla <- t2 %>% group_by(description, Phylum) %>% mutate(Phyla_relative_abundance = sum(relative_abundance)) %>% select(description, Phylum, Phyla_relative_abundance, treatment, timepoint, Material)%>% 
            distinct() %>%mutate(Phylum = ifelse(Phyla_relative_abundance<0.01, "others<0.01", Phylum)) %>%   filter ( timepoint %in% c("T1", "T6")) #%>% mutate(timepoint = ifelse(Material =="negative_c", "T1", timepoint))

ggplot(Phyla, aes(x=Material, y= Phyla_relative_abundance, fill=Phylum))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = CPCOLS)+
  theme_classic2()+  facet_grid (timepoint~treatment)


Class <- t2 %>% group_by(description, Class) %>% mutate(Class_relative_abundance = sum(relative_abundance)) %>% select(description, Class, Class_relative_abundance, treatment, timepoint, Material)%>% 
          distinct()  %>%mutate(Class = ifelse(Class_relative_abundance<0.01, "others<0.01", Class)) %>%  filter ( timepoint %in% c("T1", "T6"))


ggplot(Class, aes(x=Material, y= Class_relative_abundance, fill=Class))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = CPCOLS)+
  theme_classic()+  facet_grid (timepoint~treatment)

Order <- t2 %>% group_by(description, Order) %>% mutate(Order_relative_abundance = sum(relative_abundance)) %>% select(description, Order, Order_relative_abundance, treatment, timepoint, Material)%>% 
          distinct() %>% mutate(Order = ifelse(Order_relative_abundance<0.01, "others<0.01", Order))   %>% filter ( timepoint %in% c("T1", "T6"))


ggplot(Order, aes(x=Material, y= Order_relative_abundance, fill=Order))+
  geom_bar(stat="identity", position="stack")+ scale_color_brewer()+ theme_classic()+  facet_grid (timepoint~treatment)

Family <- t2 %>% group_by(description, Family) %>% mutate(Family_relative_abundance = sum(relative_abundance)) %>% select(description, Family, Family_relative_abundance, treatment, timepoint, Material)%>% 
  distinct() %>% mutate(Family = ifelse(Family_relative_abundance<0.01, "others<0.01", Family)) %>% filter ( timepoint %in% c("T1", "T6"))

ggplot(Family, aes(x=Material, y= Family_relative_abundance, fill=Family))+
  geom_bar(stat="identity", position="stack")+ scale_color_brewer()+
  theme_classic()+  facet_grid (timepoint~treatment)

Genus <- t2 %>% group_by(description, Genus) %>% mutate(Genus_relative_abundance = sum(relative_abundance)) %>% select(description, Genus, Genus_relative_abundance, treatment, timepoint, Material)%>%
  distinct()%>% mutate(Genus = ifelse(Genus_relative_abundance<0.01, "others<0.01", Genus))  %>% 
  filter ( timepoint %in% c("T1", "T6"))

ggplot(Genus, aes(x=Material, y= Genus_relative_abundance, fill=Genus))+
  geom_bar(stat="identity", position="stack")+ scale_color_brewer()+ 
  theme_classic()+  facet_grid (timepoint~treatment)

Species <- t2 %>% group_by(description, Species) %>% mutate(Species_relative_abundance = sum(relative_abundance)) %>% select(description, Species, Species_relative_abundance, treatment, timepoint, Material)%>% 
  distinct() %>% mutate(Species = ifelse(Species_relative_abundance<0.01, "others<0.01", Species))  %>% 
  filter ( timepoint %in% c("T1", "T6"))

ggplot(Species, aes(x=Material, y= Species_relative_abundance, fill=Species))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = CPCOLS)+
  theme_classic()+  facet_grid (timepoint~treatment)+
  theme (axis.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold") ) 

#colors_vector_to_personalize#
CPCOLS <- c("#199442", "#ED1F1F", "#F5EE2C", "#B636D6", "#3D68E0", "#EBA53D", "#00688B", "#00EE76", "#CD9B9B", "#00BFFF", "#FFF68F", "#FF7F50", "#68228B", "#ADFF2F", "#CD0000", "#0000FF", "#CD9B1D", "#FF34B3", "#BBFFFF", "#191970") 

################## heatmap on genus & no NAs #############
t2_no_na_genus <- filter(t2, !Genus==" NA")%>% group_by(description, Genus) %>% mutate(Genus_relative_abundance = sum(relative_abundance)) %>% 
  select(description, Genus, Genus_relative_abundance, treatment, timepoint, Material)%>% distinct() %>% 
  arrange(desc(Genus_relative_abundance))

t2_no_na_genus <- filter(t2_no_na, Genus %in% (unique(t2_no_na$Genus)[1:20]))


pal <- wes_palette("Zissou1", 100, type = "continuous")
ggplot(t2_no_na_genus,aes(x=description,y=Genus,fill=Genus_relative_abundance))+
  geom_tile(colour="white",size=0.25)+ labs(x="",y="")+
  geom_text(aes(label = round(Genus_relative_abundance, 3)), colour = "Black" , size = 3)+ scale_fill_gradientn(colours =  pal)+
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
bray_PCOA <- amp_ordinate(data , transform = "none", type = "PCOA", distmeasure= "bray", sample_label_size=4, sample_color_by = "treatment", sample_label_by= "Material",sample_shape_by ="timepoint", species_plot = F, detailed_output=T)

#venn
venn <- amp_venn(data = data, group_by = "treatment", cut_f = 50, detailed_output = T)
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

