#####################################################################################
#Diversity calculation CALCULATION Honours project 2022/23
#####################################################################################
#############################################################
# Clean-up the memory and start a new session
#############################################################
#https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################
library("ggplot2")
library("vegan")
library("phyloseq") 
library("tidyverse") 
library("RColorBrewer") 

#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()

#set the working directory
setwd("")
getwd()

#########################
##Phyloseq
########################
#We need design, the data_count and the data_taxonomy file) Optional: data_tree (phylogenentic tree)

#The phyloseq package recapitulates 
#a) The OTU Table counts. This is a n rows x m columns dataframe, where n represents the individual isolates, with their associated abundance counts based on morphotypes counts
#b) The taxonomy information. For each of the n OTUs, a representative sequence (based in morphotypes was Sanger seqeunced) 
#c) The mapping file (often reported as a design file or metadata): it says what the samples are, from where they come from and also it includes additional attributes that can be used in the analysis 


##Manually created OTU table based on ITS identification and count of morphotypes

OTU_table<- read.delim("abundance_otu_ce.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)


#design file
design <- read.delim("MAP3_BG.txt", sep = "\t", header=TRUE, row.names=1)
design

#generarte the count file
dat_count <- OTU_table[ ,1:10]
rownames(dat_count) <- rownames(OTU_table)
dat_count[ ,1:10]
OTU_table[ ,1:10]

#The taxonomy information
#it is a tab-delimited file with 8 columns, the column headers are: OTU id; "Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species",
#and a new taxa table


Fungal_dat_tax <- read.delim ("taxonomy_ce.txt", sep = "\t", row.names=1, header=T, blank.lines.skip = TRUE)
Fungal_taxa <- tax_table(as.matrix(Fungal_dat_tax))
dim(Fungal_taxa)

#save the above file and in excel we will create a new tax table where each column represents a taxonomic rank
#write.table(dat_taxa, file="dat_taxa_test.txt", sep="\t")

####################################################
##Species Accumulation Curves (do not feature in the Thesis)
####################################################
#https://search.r-project.org/CRAN/refmans/vegan/html/specaccum.html

##gives the expected number of observed species or distinct classes as a function of sampling effort
#Number of fungal species (y-axis) vs Number of plants (x-axis)


##Data has species as rows and genotypes as columns, but we need the opposite. Thus, transpose the original data

dat_count_t<-t(dat_count)

#Species accumulation method (partial match). Method "collector" adds sites in the order they happen to be in the data, 
#"random" adds sites in random order, "exact" finds the expected (mean) species richness, 
#"coleman" finds the expected richness following Coleman et al. 1982, and 
#"rarefaction" finds the mean when accumulating individuals instead of sites. 

rare_curv<- specaccum(dat_count_t,method = "rarefaction")

#plot the curve with some predefined settings
plot(rare_curv,ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

####################################################################################################

#The OTU Table counts
Fungal_OTU <- otu_table(dat_count, taxa_are_rows=TRUE)


#The taxonomy information
#it is a tab-delimited file with 8 columns, the column headers are: OTU id; "Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species",

dim(Fungal_taxa)

#The mapping file 
Fungal_map <- sample_data(design[1:10, ])

#merge the files and create the phyloseq object
Fungal_data_phyloseq <- merge_phyloseq(Fungal_OTU, Fungal_taxa, Fungal_map)
Fungal_data_phyloseq

#inspect the generated data
Fungal_data_phyloseq
sum(colSums(otu_table(Fungal_data_phyloseq)))
dim(dat_count)
sum(colSums(dat_count))


####################################
##Figure 3.2: Beta_diversity  
####################################

#ungal_data_phyloseq<- subset_samples(Fungal_data_phyloseq, genotype == "Barke", genotype == "124_52")
#Fungal_data_phyloseq

#constrained ordinatiton: consatrained for Description

Fungal_endosphere_CAP <- ordinate(Fungal_data_phyloseq, "CAP", "bray", ~ Genotype)
plot_ordination(Fungal_data_phyloseq, Fungal_endosphere_CAP, color = "Genotype")

#assign shapes to Soil type and color to Sample type
p=plot_ordination(Fungal_data_phyloseq, Fungal_endosphere_CAP , color = "Genotype")
p = p + geom_point(size = 6, alpha = 0.80, stroke =1)
p = p + scale_colour_manual(values = c("#E69F00", "#0072B2"))
p + ggtitle("CAP Fungal endhophytic data, Bray distance samples")

#ANOVA on the axis
anova(Fungal_endosphere_CAP, permutations=5000)

#BC distance
BC <- phyloseq::distance(Fungal_data_phyloseq, "bray")
adonis2(BC ~ Genotype, data= design , permutations = 5000)


  #################################
## Figure 3.3: Alpha diversity
##################################

FE_alpha <-  estimate_richness(Fungal_data_phyloseq, measures = c("Observed", "Shannon", "Chao1")) 
FE_alpha 

FE_otu_table<-otu_table(Fungal_data_phyloseq)
#Data frame Genotype_Description 
colData = design[colnames(FE_otu_table), ]
rownames(colData)
colData


#Description 
design_genotype  <- as.data.frame(colData[, 1]) 
rownames(design_genotype) <- rownames(colData) 
colnames(design_genotype) <- c("Genotype") 
design_genotype  

###############################
#### Alpha diversity OBSERVED 
##############################
#Observed ASVs 
FE_alpha_Observed <- as.data.frame(FE_alpha[ ,1]) 
rownames(FE_alpha_Observed) <- rownames(FE_alpha) 
colnames(FE_alpha_Observed) <- c("Observed") 

#Combine the dataset sample description and Observed OTUs 
FE_alpha_Observed_TD <- cbind(design_genotype , FE_alpha_Observed) 
FE_alpha_Observed_TD <- as.data.frame(FE_alpha_Observed_TD)
FE_alpha_Observed_TD$Genotype


#Order the levels according to a defined order 
FE_alpha_Observed_TD$Genotype <- ordered(FE_alpha_Observed_TD$Genotype, levels=c("Barke","124_52"))  

#Plotting
with(FE_alpha_Observed_TD, boxplot(Observed~Genotype, xlab = "Genotypes", ylab = "Number of taxa", main = "Taxa observed richness", col=c("#0072B2", "#E69F00")))
with(FE_alpha_Observed_TD, stripchart(Observed~Genotype, xlab = "Genotypes", ylab = "Number of taxa", vertical=T, add=T, method ="jitter", col=c('black'), pch=16))


#ANOVA 
Observed_OTUs_stats <- t.test(FE_alpha_Observed, data =FE_alpha_Observed_TD) 
Observed_OTUs_stats  

###########################
#### Figure 3.4: Alpha diversity CHAO1 
##########################
#Chao1 ASVs 
FE_alpha_Chao1 <- as.data.frame(FE_alpha[ ,2]) 
rownames(FE_alpha_Chao1) <- rownames(FE_alpha) 
colnames(FE_alpha_Chao1) <- c("Chao1") 


#Combine the dataset sample description and Chao1 OTUs 
FE_alpha_Chao1_TD <- cbind(design_genotype , FE_alpha_Chao1) 
FE_alpha_Chao1_TD <- as.data.frame(FE_alpha_Chao1_TD) 
FE_alpha_Chao1_TD$Genotype


#Order the levels according to a defined order 
FE_alpha_Chao1_TD$Genotype <- ordered(FE_alpha_Chao1_TD$Genotype, levels=c("Barke","124_52"))  

#Plotting
with(FE_alpha_Chao1_TD, boxplot(Chao1~Genotype, xlab = "Genotypes", ylab = "Number of Taxa", main = "Taxa Chao richness", col=c("#0072B2","#E69F00", "#D55E00")))
with(FE_alpha_Chao1_TD, stripchart(Chao1~Genotype, xlab = "Genotypes", ylab = "Number of ASVs", vertical=T, add=T, method ="jitter", col=c('black'), pch=16))

#ANOVA 
Chao1_OTUs_stats <- t.test(FE_alpha_Chao1, data = FE_alpha_Chao1_TD) 
Chao1_OTUs_stats

#############################
#### Alpha diversity SHANNON (Does not feature in the Thesis)
#############################
#Shannon ASvs 
FE_alpha_Shannon <- as.data.frame(FE_alpha[ ,4]) 
rownames(FE_alpha_Shannon) <- rownames(FE_alpha) 
colnames(FE_alpha_Shannon) <- c("Shannon") 


#Combine the dataset sample description and Shannon OTUs 
FE_alpha_Shannon_TD <- cbind(design_genotype , FE_alpha_Shannon) 
FE_alpha_Shannon_TD <- as.data.frame(FE_alpha_Shannon_TD) 
FE_alpha_Shannon_TD$Genotype

#Order the levels according to a defined order 
FE_alpha_Shannon_TD$Genotype <- ordered(FE_alpha_Shannon_TD$Genotype, levels=c("Barke","52"))  

#Plotting
with(FE_alpha_Shannon_TD, boxplot(Shannon~Genotype, xlab = "Genotypes", ylab = "Number of Taxa", main = "ASVs Shannon richness",col=c("#0072B2","#E69F00")))
with(FE_alpha_Shannon_TD, stripchart(Shannon~Genotype, xlab = "Genotypes", ylab = "Number of ASVs", vertical=T, add=T, method ="jitter", col=c('black'), pch=16))


#ANOVA 
Shannon_OTUs_stats <- t.test(FE_alpha_Shannon, data = FE_alpha_Shannon_TD) 
Shannon_OTUs_stats

##############################################################################
# Figure 3.5: Top 5 Genus
#############################################################################

#agglomerate at family level
FE_Genus <-  tax_glom(Fungal_data_phyloseq, taxrank="Genera")
FE_Genus

sample_data(FE_Genus)
#transform in to relative abundance
ps_genus_0 <- transform_sample_counts(FE_Genus, function(x) x / sum(x))
#abundance of all samples plot
plot_bar(ps_genus_0, fill="Genera")
#merge samples by treatment
ps_genus_1 <- merge_samples(ps_genus_0, "Genotype")
#transform to relative abudance
ps_genus_2 <- transform_sample_counts(ps_genus_1, function(x) x / sum(x))
#plot ra of all phyla
plot_bar(ps_genus_2, fill="Genera")

sample_data(ps_genus_2)
tax_table(ps_genus_2)

#plotting based on sample type
df_genus <- psmelt(ps_genus_2)
#write.table(df_phylum, file="df_phylum.txt ", sep="\t")


top_genus <- df_genus
group_by("Genotype", Genera) %>% 
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
top_genus

top_genus <- df_genus %>%
  group_by("Genotype", Genera) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean)
top_genus



##End