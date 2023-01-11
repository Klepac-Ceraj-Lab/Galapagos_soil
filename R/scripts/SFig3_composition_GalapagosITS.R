#################################################################################
#
#                     GALAPAGOS
#                      
#    
#
# By: Vanja Klepac-Ceraj and Caroline MacVicar
# Code adapted and modified from Connie Rojas
################################################################################

## CODE FOR: generating stacked bar plots of microbial community composition for 
# BOKASHI data at the fungal phylum, order, and genus taxonomic levels


source(file="scripts/00_background.R"); #load necessary packages and specifications

################################################################################
#             1.  Load OTU table and sample metadata
#          
################################################################################
taxa=read.table("data/R_Oct_Galapagos_L6_ITS.txt", sep="\t", header=T);
meta_data=read.table("data/GalapagosMeta.txt", sep="\t", header=T);

################################################################################
#             2. Create Phylum level composition barplots                 
################################################################################
#select your taxonomic level
#1-Kingdom 2-Phylum 3-Class 4-Order 5-Family 6-Genus
taxa2<-taxa[, c(2, 7:ncol(taxa))];
colnames(taxa2)[1]="Taxa";

#order [experiments]
#meta_data$Sample_Depth<-factor(meta_data$Sample_Depth, levels=c("rhizo", "5cm","15cm"))
#meta_data$Elevation_m <- factor(meta_data$Elevation_m, levels = c(320, 520, 690), labels = c("320", "520", "690"))
#meta_data$Plant<-factor(meta_data$Plant, levels=c("guava", "scalesia", "vervain", "miconia"))
meta_data$Site_Plant<-factor(meta_data$Site_Plant, levels = c("Mirador_Guava", "Mirador_Scalesia", "Cerro Alto_Guava",
                                                              "Cerro Alto_Vervain", "el junco_Guava", "el junco_Miconia"))

#build stacked bar plots of taxa abundances


#calculate taxa relative abundances 
taxa_rel=aggregate(.~Taxa, taxa2, sum);
taxa_rel[,-1] <- lapply(taxa_rel[,-1], function(x) (x/sum(x))*100); 
print(colSums(taxa_rel[-1]));  


#keep taxa >1% relative abundance across samples
taxa_rel$AVG=rowMeans(taxa_rel[,-1]);
taxa_rel=taxa_rel[taxa_rel$AVG>1,];
taxa_rel$AVG=NULL;

#denote the rest of taxa as "Other"
newrow=c(NA, 100-colSums(taxa_rel[2:ncol(taxa_rel)])); 
taxa_rel=rbind(taxa_rel, newrow); 
taxa_rel$Taxa=as.character(taxa_rel$Taxa);
taxa_rel[nrow(taxa_rel),1]="1_Other";

#melt data frame for ggplot
pbar<-reshape2::melt(taxa_rel, id.vars="Taxa",value.name = "abun");
colnames(pbar)[2]="sample";
pbar=merge(pbar, meta_data, by="sample");

#set color palette
phy_col=c("grey", "#66c2a5","#fc8d62","#7570b3","#e78ac3",'#a6d854','#ffd92f',
          "#6baed6","#fc9272","#238b45","#525252","#e7298a","#08589e",
          "#bf812d");

#create plot (group samples by Elevation_m)
barphy=ggplot(data=pbar, 
              mapping=aes(x=plant_nr,y=abun, fill=Taxa))+
  geom_bar(stat = "identity")+
  facet_grid(Sample_Depth ~ Site_Plant, scales='free')+
  theme_bw() +
  labs(x = "Plant_Nr",
       y = "Relative Abundance (%)",
       fill="Fungal Phylum")+
  scale_fill_manual(values=phy_col)+
  theme(legend.position="right", 
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text( angle = 90,  hjust = 1 ),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        strip.text = element_text(size =8, face="bold"),
        panel.background = element_rect(size=2))+
  guides(fill=guide_legend(ncol=1));

plot(barphy);

##save image 
ggsave(filename="SFig3_galapagosITS_stacked_phyla.pdf",
       device="pdf",path="./images",
       plot=barphy,
       width=15,
       height=8,
       units="in",
       dpi=500);


################################################################################
#             3. Create Order level composition barplots                 
################################################################################
#select your taxonomic level
#1-Kingdom 2-Phylum 3-Class 4-Order 5-Family 6-Genus
taxa2<-taxa[, c(4, 7:ncol(taxa))];
colnames(taxa2)[1]="Taxa";

#calculate taxa relative abundances 
taxa_rel=aggregate(.~Taxa, taxa2, sum);
taxa_rel[,-1] <- lapply(taxa_rel[,-1], function(x) (x/sum(x))*100);
#print(colSums(taxa_rel[-1]));

#keep taxa >1% relative abundance across samples
taxa_rel$AVG=rowMeans(taxa_rel[,-1]);
taxa_rel=taxa_rel[taxa_rel$AVG>1,];
taxa_rel$AVG=NULL;

#denote the rest of taxa as "Other"
newrow=c(NA, 100-colSums(taxa_rel[2:ncol(taxa_rel)])); 
taxa_rel=rbind(taxa_rel, newrow); 
taxa_rel$Taxa=as.character(taxa_rel$Taxa);
taxa_rel[nrow(taxa_rel),1]="1_Other";

#melt data frame for ggplot
obar<-reshape2::melt(taxa_rel, id.vars="Taxa",value.name = "abun");
colnames(obar)[2]="sample";
obar=merge(obar, meta_data, by="sample");

#set color palette
ord_col=c("grey", "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", 
          "#737373","#117744","#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411","#AA7744", 
          "#DDAA77", "#771122", "#AA4455", "#DD7788","black","coral", "cornsilk", "yellow", "ivory", "lightskyblue1", "slateblue1", "orange2", "blue4", "lightgreen" )


#create plot (group samples by Site_Plant)
barord=ggplot(data=obar, 
              mapping=aes(x=plant_nr,y=abun, fill=Taxa))+
  geom_bar(stat = "identity")+
  facet_grid(Sample_Depth ~ Site_Plant, scales='free')+
  theme_bw() +
  labs(x = "Plant_Nr",
       y = "Relative Abundance (%)",
       fill="Fungal Order")+
  scale_fill_manual(values=ord_col)+
  theme(legend.position="right", 
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text( angle = 90,  hjust = 1 ),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        strip.text = element_text(size =8, face="bold"),
        panel.background = element_rect(size=2))+
  guides(fill=guide_legend(ncol=1));

plot(barord);

##save image 
ggsave(filename="SFig3_GalapagosITS_stacked_order.pdf",
       device="pdf",path="./images",
       plot=barord,
       width=15,
       height=8,
       units="in",
       dpi=500);


################################################################################
#             4. Create Genus level composition barplots                 
################################################################################
#select your taxonomic level
#1-Kingdom 2-Phylum 3-Class 4-Order 5-Family 6-Genus
taxa2<-taxa[, c(6, 7:ncol(taxa))];
colnames(taxa2)[1]="Taxa";

#calculate taxa relative abundances 
taxa_rel=aggregate(.~Taxa, taxa2, sum);
taxa_rel[,-1] <- lapply(taxa_rel[,-1], function(x) (x/sum(x))*100);
#print(colSums(taxa_rel[-1]));

#keep taxa >1% relative abundance across samples
taxa_rel$AVG=rowMeans(taxa_rel[,-1]);
taxa_rel=taxa_rel[taxa_rel$AVG>1,];
taxa_rel$AVG=NULL;

#denote the rest of taxa as "Other"
newrow=c(NA, 100-colSums(taxa_rel[2:ncol(taxa_rel)])); 
taxa_rel=rbind(taxa_rel, newrow); 
taxa_rel$Taxa=as.character(taxa_rel$Taxa);
taxa_rel[nrow(taxa_rel),1]="1_Other";

#melt data frame for ggplot
gbar<-reshape2::melt(taxa_rel, id.vars="Taxa",value.name = "abun");
colnames(gbar)[2]="sample";
gbar=merge(gbar, meta_data, by="sample");

#set color palette
gen_col=c("grey", "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", 
          "#737373","#117744","#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", 
          "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","brown","coral", "orange4", "purple3", "lightblue", "slateblue1" )


#create plot (group samples by Site_Plant)
bargen=ggplot(data=gbar, 
              mapping=aes(x=plant_nr,y=abun, fill=Taxa))+
  geom_bar(stat = "identity")+
  facet_grid(Sample_Depth ~ Site_Plant, scales='free')+
  theme_bw() +
  labs(x = "Plant_Nr",
       y = "Relative Abundance (%)",
       fill="Fungal Genus")+
  scale_fill_manual(values=gen_col)+
  theme(legend.position="right", 
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text( angle = 90,  hjust = 1 ),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank(),
        strip.text = element_text(size =8, face="bold"),
        panel.background = element_rect(size=2))+
  guides(fill=guide_legend(ncol=1));

plot(bargen);

##save image 
ggsave(filename="SFig3_Galapagos_ITS_stacked_gen.pdf",
       device="pdf",path="./images",
       plot=bargen,
       width=15,
       height=8,
       units="in",
       dpi=500);


################################################################################
#             5. Download table of taxa abundances                 
################################################################################
#calculate abundances
taxa2<-taxa[, c(6, 7:ncol(taxa))];
colnames(taxa2)[1]="Taxa";
taxa_rel=aggregate(.~Taxa, taxa2, sum);
taxa_rel[,-1] <- lapply(taxa_rel[,-1], function(x) (x/sum(x))*100);

#add sample Site_Plant labels
colnames(taxa_rel)=meta_data$Site_Plant[match(colnames(taxa_rel),
                                              meta_data$sample)];
colnames(taxa_rel)[1]=1;
taxa_rel=taxa_rel[, order(names(taxa_rel))];
colnames(taxa_rel)[1]="Taxa";

#write file
write.csv(taxa_rel, file="data/01_GalapagosITS_genus.csv",row.names=F);


#calculate abundances
taxa2<-taxa[, c(4, 7:ncol(taxa))];
colnames(taxa2)[1]="Taxa";
taxa_rel=aggregate(.~Taxa, taxa2, sum);
taxa_rel[,-1] <- lapply(taxa_rel[,-1], function(x) (x/sum(x))*100);

#add sample Site_Plant labels
#colnames(taxa_rel)=meta_data$Site_Plant[match(colnames(taxa_rel),
#                                              meta_data$sample)];
colnames(taxa_rel)[1]=1;
taxa_rel=taxa_rel[, order(names(taxa_rel))];
colnames(taxa_rel)[1]="Taxa_Order";

#write file
write.csv(taxa_rel, file="data/01_GalapagosITS_order.csv",row.names=F);
