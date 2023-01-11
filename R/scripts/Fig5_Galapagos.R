#################################################################################
#             Microclimate is a strong predictor of the native and
#           invasive plant-associated soil microbiota on San Cristóbal Island, 
#                         Galápagos archipelago
#
# By: Caroline McVicar and Vanja Klepac-Ceraj
# code adapted from C. Rojas
################################################################################

## CODE FOR: 
##  - generating PCoA plots based on microbial community dissimilarity matrices
##  - Bray-Curtis & Jaccard 
##  - testing significance of patterns using PERMANOVAs
##  - creating boxplots of average Bray-Curtis dissimilarity

################################################################################
#             Configure the workspace for subsequent R project scripts                 
################################################################################

#set conditions for R session
rm(list=ls());
options(scipen=999);
options(stringsAsFactors = FALSE);
options(show.error.locations = TRUE); # show line numbers on error

#load necessary packages
library(pacman);
pacman::p_load("car","MASS","dplyr","tidyr","reshape2","vegan","ggplot2",
               "lme4","lmtest","multcomp","FSA","tidyverse");


################################################################################
#             Communicate the R version and package versions to reader                 
################################################################################
print("This code was developed with R version 4.0.3 (2020-10-10)");

print("The packages used and their versions were: multcomp_1.4-15| lmtest_0.9-38| 
lme4_1.1-26| ggplot2_3.3.3| vegan_2.5-7| reshape2_1.4.4| tidyr_1.1.2| dplyr_1.0.3| 
MASS_7.3-53| car_3.0-10| pacman_0.5.1| FSA_0.8.32");


################################################################################
#                Load OTU table and sample metadata
#           Chloroplast, Mitochondria, and Unknown were removed
################################################################################

# 16S data:
taxa16S=read.table("data/Galapagos_L6_16S.txt", sep="\t", header=T);
meta_data16S=read.table("data/GalapagosMeta.txt", sep="\t", header=T);
meta_data16S=meta_data16S[order(meta_data16S$sample),];

meta_data16S$Site<-factor(meta_data16S$Site, levels=c("Mirador", "Cerro Alto","el junco"));
meta_data16S$Native <- factor(meta_data16S$Native, labels = c("no", "yes"));
meta_data16S$Sample_Depth <- factor(meta_data16S$Sample_Depth, labels = c("rhizo", "5cm", "15cm"));

# ITS data:
taxaITS=read.table("data/Galapagos_L6_ITS.txt", sep="\t", header=T);
meta_dataITS=read.table("data/GalapagosMeta.txt", sep="\t", header=T);
meta_dataITS=meta_dataITS[order(meta_dataITS$sample),];
#need to remove rows with samples C1_B, C3_C and M6_A
#meta_dataITS %>%
meta_dataITS = subset(meta_dataITS, sample != "C1_B")
meta_dataITS = subset(meta_dataITS, sample != "C3_C")
meta_dataITS = subset(meta_dataITS, sample != "M6_A")
  
meta_dataITS$Site<-factor(meta_dataITS$Site, levels=c("Mirador", "Cerro Alto","el junco"))
meta_dataITS$Native <- factor(meta_dataITS$Native, labels = c("no", "yes"))
meta_dataITS$Sample_Depth <- factor(meta_dataITS$Sample_Depth, labels = c("rhizo", "5cm", "15cm"))


################################################################################
#              Generate Bray-Curtis distance matrix
################################################################################

#transpose OTU table so OTUs are columns -- 16S data
taxadf16S=taxa16S[,7:ncol(taxa16S)];
taxadf16S=t(taxadf16S);
taxadf16S=taxadf16S[order(row.names(taxadf16S)),];
 
###BRAY-CURTIS distance
bray16S<-apply(taxadf16S, 1, function(i) (i/sum(i)));
bray16S=as.data.frame(t(bray16S));
#print(rowSums(bray));
bray.dist16S=vegdist(bray16S, method="bray");

#transpose OTU table so OTUs are columns -- ITS data
taxadfITS=taxaITS[,7:ncol(taxaITS)];
taxadfITS=t(taxadfITS);
taxadfITS=taxadfITS[order(row.names(taxadfITS)),];

###BRAY-CURTIS distance
brayITS<-apply(taxadfITS, 1, function(i) (i/sum(i)));
brayITS=as.data.frame(t(brayITS));
#print(rowSums(bray));
bray.distITS=vegdist(brayITS, method="bray");


################################################################################
#             3. Conduct PERMANOVAs by Site*soil_depth*species_invasiveness
################################################################################
#does microbial community structure vary with the site (microclimate), soil depth and 
#species invasiveness for samples collected at St. Cristobal, Galapagos? 
#Test for the effects as well as their interaction.
set.seed(19760403)
#Bray-Curtis 16S
print(adonis(bray.dist16S~Site * Native * Sample_Depth,             
             data=meta_data16S, method = "bray",         
             permutations = 9999)); 

#Bray-Curtis ITS
print(adonis(bray.distITS~Site * Native * Sample_Depth,             
             data=meta_dataITS, method = "bray",         
             permutations = 9999));

################################################################################
#             5. Generate NMDS plots
################################################################################

NMDS16S <- metaMDS(bray.dist16S)
NMDS16S$points
scores(NMDS16S) %>%
  as_tibble(rownames = "sample") %>%
  inner_join(., meta_data16S, by="sample" ) %>%
  
ggplot(aes(x=NMDS1, y=NMDS2, color = Site))+
  geom_point(size = 3, aes(shape = Sample_Depth,alpha = Native, fill = Site, show.legend = FALSE)) +
  stat_ellipse() +
  theme_bw() +
  scale_shape_manual(values = c(23, 24, 22)) +
  scale_alpha_manual(values = c("no" = 0.5, "yes" = 1)) +
  scale_color_manual(values = c("Mirador" = "orange2", "Cerro Alto" = "#008E51", "el junco" = "#043FFF")) +
  labs(title = "BACT NMDS Site Plot", x = "NMDS1", y="NMDS2")
  
ggsave(filename="Figure_5A_16SNMDS.pdf",
       device="pdf",path="./images",
       width=6,
       height=5,
       units="in",
       dpi=500);

#NMDS - ITS data
set.seed(1)
NMDSITS <- metaMDS(bray.distITS)
NMDSITS$points
scores(NMDSITS) %>%
  as_tibble(rownames = "sample") %>%
  inner_join(., meta_dataITS, by="sample" ) %>%
  
  ggplot(aes(x=NMDS1, y=NMDS2, color = Site))+
  geom_point(size = 3, aes(shape = Sample_Depth,alpha = Native, fill = Site, show.legend =)) +
  stat_ellipse() +
  theme_bw() +
  scale_shape_manual(values = c(23, 24, 22)) +
  scale_alpha_manual(values = c("no" = 0.5, "yes" = 1)) +
  scale_color_manual(values = c("Mirador" = "orange2", "Cerro Alto" = "#008E51", "el junco" = "#043FFF")) +
  labs(title = "ITS NMDS Site Plot", x = "NMDS1", y="NMDS2")

ggsave(filename="Figure_5B_ITSNMDS.pdf",
       device="pdf",path="./images",
       width=6,
       height=5,
       units="in",
       dpi=500);


################################################################################
#             6. differential relative abundance
################################################################################


