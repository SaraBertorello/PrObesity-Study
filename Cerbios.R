############### CERBIOS ############### 
############### LOAD PACKAGES ############### 
{ 
  library("phyloseq")
  library("ggplot2")
  library("ggvenn")
  library("egg")
  library("ggh4x")
  library("ggbreak") 
  library("vegan")
  library("stringr")
  library("xlsx")  
  library("qiime2R")
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(readxl)
  library(dplyr)
  library(DESeq2)
  library(ggpubr)
  library(devtools)
  library(tidyr)
  library(qiime2R)
  library(dunn.test)
  library(pairwiseAdonis)
  library(cluster)
  library(ggh4x) 
}

options(scipen=100)

setwd('../CERBIOS')

dir.create("Results_CERBIOS")
dir.create("Results_CERBIOS/Alpha_Div")
dir.create("Results_CERBIOS/Beta_Div")
dir.create("Results_CERBIOS/Control")
dir.create("Results_CERBIOS/Deseq2")
dir.create("Results_CERBIOS/Abundances")
dir.create("Results_CERBIOS/PCoA")
dir.create("Results_CERBIOS/TOP_5_Phyla")
dir.create("Results_CERBIOS/TOP_5_Genus")


############### LOAD OBJECTS ############### 
#taxonomy
data<-qza_to_phyloseq(taxonomy="taxonomy.qza")
View(data)
write.csv(x = data, file='taxonomy.csv')
data<-read.csv('taxonomy.csv')
r<- data$X
data<-data[,-1]
row.names(data)<-r

#otu table
data1<-qza_to_phyloseq(features="Featuredada.qza")
View(data1)
colnames(data1)
write.csv(x = data1, file='otu_table.csv')
colnames(data1) <- sub(".*(BC-\\d+).*", "\\1", colnames(data1))
colnames(data1) <- gsub("-", "", colnames(data1))
print(colnames(data1))

#metadata
metadata<-read.xlsx('metadata.xlsx', sheetIndex = 1)
colnames(metadata)
metadata<-metadata[1:72,]
rows<-metadata$Sample
row.names(metadata)<-rows
identical(row.names(metadata),colnames(data1))

data<- phyloseq(otu_table(as.matrix(data1), taxa_are_rows = TRUE), tax_table(as.matrix(data)))
head(otu_table(data))
tax_table(data)
sample_names()
sample_data(data)<-metadata
sample_data(data)

############### preliminary PCoA ############### 
data.prop_first <- transform_sample_counts(data, function(x) x/sum(x)*100)
otu_table(data.prop_first)[is.na(otu_table(data.prop_first))] <- 0
which(is.na(otu_table(data.prop_first)))

#### on ASV - T0
sample_data(data)
stool_pcoa<-subset_samples(data.prop_first, Time=='T0')
data.sqrt_prop<-transform_sample_counts(stool_pcoa, sqrt)
DistEU <- phyloseq::distance(data.sqrt_prop, method = "euclidean") 
ordEU <- ordinate(data.sqrt_prop, method = "PCoA", distance = DistEU)

eigval<-ordEU$values$Eigenvalues   
eigval<- round((eigval/sum(eigval))*100,1) #to extract value of 2 PC axes

color=c("INTERV" = "cyan", "CONTR" = "green3")
plot_ordination(data.sqrt_prop, ordEU, color = "Condition") +
  geom_point(size=2.5) + 
  geom_text(aes(label = Sample),color="black",show.legend=F, size=2) +
  theme_classic() +
  stat_ellipse() + 
  labs(title="PCoA with Euclidean distance pre filter on ASV percent abundance T0", 
       color="Condition", 
       x=paste("PC1: ",eigval[1],"% variation"), 
       y=paste("PC2: ",eigval[2],"% variation")) +
  scale_color_manual(values = color)

ggsave(file="Results_CERBIOS/PCoA/PCoA_Euclidean_ASV_pre_filter_T0.png", width = 7, height = 5, dpi=300)

#### on ASV - T10
stool_pcoa<-subset_samples(data.prop_first, Time=='T1')
data.sqrt_prop<-transform_sample_counts(stool_pcoa, sqrt)
sample_data(data.sqrt_prop)
DistEU <- phyloseq::distance(data.sqrt_prop, method = "euclidean") 
ordEU <- ordinate(data.sqrt_prop, method = "PCoA", distance = DistEU)

eigval<-ordEU$values$Eigenvalues   
eigval<- round((eigval/sum(eigval))*100,1) #to extract value of 2 PC axes

plot_ordination(data.sqrt_prop, ordEU, color = "Condition") +
  geom_point(size=2.5) + 
  geom_text( aes(label = Sample),color="black",show.legend=F, size=2) +
  theme_classic() +
  stat_ellipse() + 
  labs(title="PCoA with Euclidean distance pre filter on ASV percent abundance T10", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation")) +
  scale_color_manual(values = color)

ggsave(file="Results_CERBIOS/PCoA/PCoA_Euclidean_ASV_pre_filter_T10.png", width = 7, height = 5, dpi=300)

############### FILTERING AND CLEANING ##############
sample_data (data)
data<-subset_samples(data, Sample!='BC3')
data<-subset_samples(data, Sample!='BC19')
data<-subset_samples(data, Sample!='BC39')
data<-subset_samples(data, Sample!='BC55')
sample_data (data)

# control Domain
tax_table(data)
data.temp <- tax_glom(data, taxrank = "Kingdom", NArm = F)
cbind(otu_table(data.temp),tax_table(data.temp)[,1])
write.csv2(cbind(tax_table(data.temp)[,1], otu_table(data.temp)), file="Results_CERBIOS/Control/Kingdom_proportions.csv", quote=F, row.names = F)

tax_table(data.temp)
otu_table(data.temp)

data<-subset_taxa(data, Kingdom=="d__Bacteria")
length(which(tax_table(data)=="d__Bacteria")) 
data <- subset_taxa(data, Order!=" Chloroplast")

data.temp1 <- tax_glom(data, taxrank = "Genus", NArm = F)

data.temp1 <- transform_sample_counts(data.temp1, function(x) (x/sum(x)) *100) # relative abundance
otu_table(data.temp1)[is.na(otu_table(data.temp1))] <- 0
head(otu_table(data.temp1))
sample_data(data.temp1)

original_taxa_number<-length(taxa_names(data.temp1))
original_taxa_number 

data.bad<-filter_taxa(data.temp1, function(x) mean(x)<=0.005, TRUE )  # manteining only those with rel abundance lower than 0.005
unique(as.character(tax_table(data.bad)[,"Genus"]))
otu_table(data.bad)
write.csv2(cbind(tax_table(data.bad), otu_table(data.bad)), file="Results_CERBIOS/Control/Contaminant_genera_under_0005_average.csv", quote=F, row.names = F)

data.filter<-filter_taxa(data.temp1, function(x) mean(x)>=0.005, TRUE )
head(otu_table(data.filter))
head(tax_table(data.filter))

if(! "filter_marker" %in% ls() ){
  data.original<-data
  data<-subset_taxa(data, Genus %in% as.character(tax_table(data.filter)[ ,"Genus"])) 
  filter_marker<- "Filtering done"
  rm(data.temp, data.temp1)
} else { 
  cat("\n The dataset is already filtered \n\n")
}

save(data, file = "data_cerbios_filtered.RData") # after all filter
otu_table(data)
setwd('SARA/CERBIOS')
load(file = "data_cerbios_filtered.RData")
sample_data(data)

############################## RARECURVE ##############################################

evalslopes<-function(x,names,lim=0.5,t=10,cex=0.5) {
  #x: the rarefaction curve as generated by rarecurve (with label=F)
  #lim: the threshold of the slope value to accept saturation
  #b: how long the rarefaction tail should be evaluated (e.g. the last 10 points)
  #names: the labels (the same used of the original samples (and in the same order!!!)
  sat<-0
  for (i in 1:length(x)) {
    v<-as.vector(x[[i]])
    dx<-as.numeric(gsub("N","",names(x[[1]])[2]))-1
    check<-"red"
    if(length(x) < 10) {
      points(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),v[1],col="cyan",pch=16,cex=1)
    } else {
      #the slope is estimated (average) at the last b rarefaction points
      slope<-mean(diff(v[(length(v)-t):length(v)])/dx)
      if(slope<lim) {
        check<-"blue"
        sat = sat+1
      }
      cat(i,slope,check,"\n")
      # text(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),v[1],labels=slope,col=check,cex=0.5)
      # points(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),v[1],col=check,pch=16,cex=1)
      text(as.numeric(rev(gsub("N","",names(x[[i]])))[1]),rev(v)[1],col=check,pch=16,cex=cex,labels=names[i])
    }
  }
  legend("bottomright",paste(sat,"saturated samples"),bty="n")
}

png(file="Results_CERBIOS/Control/Rarefaction_curve.png",width=3000,height=2100, res=300)
r<-rarecurve(t(as(otu_table(data),"matrix")), step=100,label=F)
evalslopes(r,sample_names(data),lim=0.001,cex=1)
dev.off() 

sample_data(data) #all samples saturated greatly
########################## ABUNDANCES #################################
{
  data.phy <- tax_glom(data, taxrank = "Phylum", NArm = F)
  data.genus <- tax_glom(data, taxrank = "Genus", NArm = F)
  
  data.prop <- transform_sample_counts(data, function(x) x/sum(x)*100)
  otu_table(data.prop)[is.na(otu_table(data.prop))] <- 0
  data.phy.prop <- transform_sample_counts(data.phy, function(x) x/sum(x)*100)
  otu_table(data.phy.prop)[is.na(otu_table(data.phy.prop))] <- 0
  data.genus.prop <- transform_sample_counts(data.genus, function(x) x/sum(x)*100)
  otu_table(data.genus.prop)[is.na(otu_table(data.genus.prop))] <- 0
  
}
which(is.na(otu_table(data.prop)))
head( as.data.frame(tax_table(data.genus)) ) 
head( cbind(as(otu_table(data.phy.prop),"matrix"), as(tax_table(data.phy),"matrix")) )

write.csv2( cbind(as(otu_table(data),"matrix"), as(tax_table(data),"matrix")), file="Results_CERBIOS/Abundances/ASV_abundances.csv", quote=F)
write.csv2( cbind(as(otu_table(data.phy.prop),"matrix"), as(tax_table(data.phy),"matrix")), file="Results_CERBIOS/Abundances/Phyla_abundances.csv", quote=F)
write.csv2( cbind(as(otu_table(data.genus.prop),"matrix"), as(tax_table(data.genus),"matrix")), file="Results_CERBIOS/Abundances/Genera_abundances.csv", quote=F)

########################### PCoA after filter ########################
#### on ASV - T0
data.sqrt<-transform_sample_counts(data.prop, sqrt) # sqrt of proportion
stool_pcos_post<-subset_samples(data.sqrt, Time=='T0')

DistBC = phyloseq::distance(stool_pcos_post, method = "euclidean")
ordBC = ordinate(stool_pcos_post, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues   
eigval<- round((eigval/sum(eigval))*100,1) #to extract value of 2 PC axes

plot_ordination(stool_pcos_post, ordBC, color = "Condition") +
  geom_point(size=2.5) + 
  #geom_text(aes(label = Sample), color="black",show.legend=F, size=2) +
  theme_classic() +
  stat_ellipse() + 
  labs(title="PCoA with Euclidean distance after filter on ASV percent abundance T0 ", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation")) +
  scale_color_manual(values = color)

ggsave(file="Results_CERBIOS/PCoA/PCoA_Euclidean_ASV_after_filter_T0.png", width = 7, height = 5, dpi=300)

#### on ASV - T1
data.sqrt<-transform_sample_counts(data.prop, sqrt) # sqrt of proportion
stool_pcos_post<-subset_samples(data.sqrt, Time=='T1')

DistBC = phyloseq::distance(stool_pcos_post, method = "euclidean")
ordBC = ordinate(stool_pcos_post, method = "PCoA", distance = DistBC)
eigval<-ordBC$values$Eigenvalues   
eigval<- round((eigval/sum(eigval))*100,1) #to extract value of 2 PC axes

plot_ordination(stool_pcos_post, ordBC, color = "Condition") +
  geom_point(size=2.5) + 
  #geom_text(aes(label = Sample_number), color="black",show.legend=F, size=2) +
  theme_classic() +
  stat_ellipse() + 
  labs(title="PCoA with Euclidean distance after filter on ASV percent abundance T10", 
       color="Condition", x=paste("PC1: ",eigval[1],"% variation"), y=paste("PC2: ",eigval[2],"% variation")) +
  scale_color_manual(values = color)

ggsave(file="Results_CERBIOS/PCoA/PCoA_Euclidean_ASV_after_filter_T10.png", width = 7, height = 5, dpi=300)

######################  TOP 5 PHYLA e GENUS  ################################################################
###################### PHYLA #############################
### INTERV T0
sample_data (data.phy.prop)
stool.phyla.prop<-subset_samples(data.phy.prop, Time=='T0'& Condition=='INTERV')
sample_data(stool.phyla.prop)

top5_phyla <- names(sort(taxa_sums(stool.phyla.prop), decreasing=TRUE))[1:5] 
prune.dat_top5_phyla <- prune_taxa(top5_phyla, stool.phyla.prop)
otu_table(prune.dat_top5_phyla)

tax_table(stool.phyla.prop)
others_phyla<-taxa_names(stool.phyla.prop)
others_phyla
others_phyla<-others_phyla[!(others_phyla %in% top5_phyla)]
prune.others_phyla<- prune_taxa(others_phyla, stool.phyla.prop)

tabella_top_phyla<-psmelt(prune.dat_top5_phyla)
tabella_others_phyla<-psmelt(prune.others_phyla)
tabella_others_phyla$Phylum<-"Others"
tabella_phyla<-rbind.data.frame(tabella_top_phyla,tabella_others_phyla)
phyla_hc<-tabella_top_phyla

write.csv(file = paste0("Results_CERBIOS/TOP_5_Phyla/TOP5_phyla_INTERV_T0.csv"), row.names = F,
          cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top5_phyla),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top5_phyla))[["Phylum"]]))

fill_color_5<-c('Proteobacteria'="coral",'Verrucomicrobiota'="springgreen3",'Bacteroidota'="gold3",'Firmicutes'="firebrick3",'Actinobacteriota'="deepskyblue3",'Others' ="darkslategray3") 
ggplot(data=tabella_phyla, aes(x=Sample, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=14) + 
  scale_fill_manual(values= fill_color_5) +
  #scale_fill_manual(values=c("#DF536B", "#61D04F" ,"#2297E6" ,"#28E2E5" ,"#CD0BBC", "#F5C710")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Subject", y="Relative abundance", title = paste0("Five most abundant Phyla Intervention T0"))

ggsave(filename = paste0("Results_CERBIOS/TOP_5_Phyla/Bar_plot_5_top_phyla_INTERV_T0.png"), width = 9, height = 5, dpi = 300)

###  INTERV T1
sample_data (data.phy.prop)
stool.phyla.prop<-subset_samples(data.phy.prop, Time=='T1'& Condition=='INTERV')
sample_data(stool.phyla.prop)

top5_phyla <- names(sort(taxa_sums(stool.phyla.prop), decreasing=TRUE))[1:5] 
prune.dat_top5_phyla <- prune_taxa(top5_phyla, stool.phyla.prop)
otu_table(prune.dat_top5_phyla)

tax_table(stool.phyla.prop)
others_phyla<-taxa_names(stool.phyla.prop)
others_phyla
others_phyla<-others_phyla[!(others_phyla %in% top5_phyla)]
prune.others_phyla<- prune_taxa(others_phyla, stool.phyla.prop)

tabella_top_phyla<-psmelt(prune.dat_top5_phyla)
tabella_others_phyla<-psmelt(prune.others_phyla)
tabella_others_phyla$Phylum<-"Others"
tabella_phyla<-rbind.data.frame(tabella_top_phyla,tabella_others_phyla)
phyla_hc<-tabella_top_phyla

write.csv(file = paste0("Results_CERBIOS/TOP_5_Phyla/TOP5_phyla_INTERV_T1.csv"), row.names = F,
          cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top5_phyla),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top5_phyla))[["Phylum"]]))

fill_color_5<-c('Proteobacteria'="coral",'Verrucomicrobiota'="springgreen3",'Bacteroidota'="gold3",'Firmicutes'="firebrick3",'Actinobacteriota'="deepskyblue3",'Others' ="darkslategray3")
ggplot(data=tabella_phyla, aes(x=Sample, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=14) + 
  scale_fill_manual(values= fill_color_5) +
  #scale_fill_manual(values=c("#DF536B", "#61D04F" ,"#2297E6" ,"#28E2E5" ,"#CD0BBC", "#F5C710")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Subject", y="Relative abundance", title = paste0("Five most abundant Phyla Intervention T1"))

ggsave(filename = paste0("Results_CERBIOS/TOP_5_Phyla/Bar_plot_5_top_phyla_INTERV_T1.png"), width = 9, height = 5, dpi = 300)

### CONTR T0
sample_data (data.phy.prop)
stool.phyla.prop<-subset_samples(data.phy.prop, Time=='T0'& Condition=='CONTR')
sample_data(stool.phyla.prop)

top5_phyla <- names(sort(taxa_sums(stool.phyla.prop), decreasing=TRUE))[1:5] 
prune.dat_top5_phyla <- prune_taxa(top5_phyla, stool.phyla.prop)
otu_table(prune.dat_top5_phyla)

tax_table(stool.phyla.prop)
others_phyla<-taxa_names(stool.phyla.prop)
others_phyla
others_phyla<-others_phyla[!(others_phyla %in% top5_phyla)]
prune.others_phyla<- prune_taxa(others_phyla, stool.phyla.prop)

tabella_top_phyla<-psmelt(prune.dat_top5_phyla)
tabella_others_phyla<-psmelt(prune.others_phyla)
tabella_others_phyla$Phylum<-"Others"
tabella_phyla<-rbind.data.frame(tabella_top_phyla,tabella_others_phyla)
phyla_hc<-tabella_top_phyla

write.csv(file = paste0("Results_CERBIOS/TOP_5_Phyla/TOP5_phyla_CONTR_T0.csv"), row.names = F,
          cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top5_phyla),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top5_phyla))[["Phylum"]]))

fill_color_5<-c('Proteobacteria'="coral",'Verrucomicrobiota'="springgreen3",'Bacteroidota'="gold3",'Firmicutes'="firebrick3",'Actinobacteriota'="deepskyblue3",'Others' ="darkslategray3") 
ggplot(data=tabella_phyla, aes(x=Sample, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=14) + 
  scale_fill_manual(values= fill_color_5) +
  #scale_fill_manual(values=c("#DF536B", "#61D04F" ,"#2297E6" ,"#28E2E5" ,"#CD0BBC", "#F5C710")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Subject", y="Relative abundance", title = paste0("Five most abundant Phyla Control T0"))

ggsave(filename = paste0("Results_CERBIOS/TOP_5_Phyla/Bar_plot_5_top_phyla_CONTR_T0.png"), width = 9, height = 5, dpi = 300)

###  CONTR T1
sample_data (data.phy.prop)
stool.phyla.prop<-subset_samples(data.phy.prop, Time=='T1'& Condition=='CONTR')
sample_data(stool.phyla.prop)

top5_phyla <- names(sort(taxa_sums(stool.phyla.prop), decreasing=TRUE))[1:5] 
prune.dat_top5_phyla <- prune_taxa(top5_phyla, stool.phyla.prop)
otu_table(prune.dat_top5_phyla)

tax_table(stool.phyla.prop)
others_phyla<-taxa_names(stool.phyla.prop)
others_phyla
others_phyla<-others_phyla[!(others_phyla %in% top5_phyla)]
prune.others_phyla<- prune_taxa(others_phyla, stool.phyla.prop)

tabella_top_phyla<-psmelt(prune.dat_top5_phyla)
tabella_others_phyla<-psmelt(prune.others_phyla)
tabella_others_phyla$Phylum<-"Others"
tabella_phyla<-rbind.data.frame(tabella_top_phyla,tabella_others_phyla)
phyla_hc<-tabella_top_phyla

write.csv(file = paste0("Results_CERBIOS/TOP_5_Phyla/TOP5_phyla_CONTR_T1.csv"), row.names = F,
          cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top5_phyla),1,mean)), "phylum"= as.data.frame(tax_table(prune.dat_top5_phyla))[["Phylum"]]))

fill_color_5<-c('Proteobacteria'="coral",'Verrucomicrobiota'="springgreen3",'Bacteroidota'="gold3",'Firmicutes'="firebrick3",'Actinobacteriota'="deepskyblue3",'Others' ="darkslategray3")
ggplot(data=tabella_phyla, aes(x=Sample, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=14) + 
  scale_fill_manual(values= fill_color_5) +
  #scale_fill_manual(values=c("#DF536B", "#61D04F" ,"#2297E6" ,"#28E2E5" ,"#CD0BBC", "#F5C710")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Subject", y="Relative abundance", title = paste0("Five most abundant Phyla Control T1"))

ggsave(filename = paste0("Results_CERBIOS/TOP_5_Phyla/Bar_plot_5_top_phyla_CONTR_T1.png"), width = 9, height = 5, dpi = 300)

######## GENERA ########
###  INTERV T0
sample_data (data.genus.prop)
stool.Genus.prop<-subset_samples(data.genus.prop, Time=='T0'& Condition=='INTERV')
sample_data(stool.Genus.prop)

top5_Genus <- names(sort(taxa_sums(stool.Genus.prop), decreasing=TRUE))[1:5] 
prune.dat_top5_Genus <- prune_taxa(top5_Genus, stool.Genus.prop)
otu_table(prune.dat_top5_Genus)

tax_table(stool.Genus.prop)
others_Genus<-taxa_names(stool.Genus.prop)
others_Genus
others_Genus<-others_Genus[!(others_Genus %in% top5_Genus)]
prune.others_Genus<- prune_taxa(others_Genus, stool.Genus.prop)

tabella_top_Genus<-psmelt(prune.dat_top5_Genus)
tabella_others_Genus<-psmelt(prune.others_Genus)
tabella_others_Genus$Genus<-"Others"
tabella_Genus<-rbind.data.frame(tabella_top_Genus,tabella_others_Genus)
Genus_hc<-tabella_top_Genus

write.csv(file = paste0("Results_CERBIOS/TOP_5_Genus/TOP5_Genus_INTERV_T0.csv"), row.names = F,
          cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top5_Genus),1,mean)), "Genus"= as.data.frame(tax_table(prune.dat_top5_Genus))[["Genus"]]))

fill_color_5<-c('Bifidobacterium'="coral",'Collinsella'="springgreen3",'Faecalibacterium'="gold3",'Blautia'="firebrick3",'Bacteroides'="deepskyblue3",'Others' ="darkslategray3") 
ggplot(data=tabella_Genus, aes(x=Sample, y=Abundance, fill=Genus)) + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=14) + 
  scale_fill_manual(values= fill_color_5) +
  #scale_fill_manual(values=c("#DF536B", "#61D04F" ,"#2297E6" ,"#28E2E5" ,"#CD0BBC", "#F5C710")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Subject", y="Relative abundance", title = paste0("Five most abundant Genus Intervention T0"))

ggsave(filename = paste0("Results_CERBIOS/TOP_5_Genus/Bar_plot_5_top_Genus_INTERV_T0.png"), width = 9, height = 5, dpi = 300)

###  INTERV T1
sample_data (data.genus.prop)
stool.Genus.prop<-subset_samples(data.genus.prop, Time=='T1'& Condition=='INTERV')
sample_data(stool.Genus.prop)

top5_Genus <- names(sort(taxa_sums(stool.Genus.prop), decreasing=TRUE))[1:5] 
prune.dat_top5_Genus <- prune_taxa(top5_Genus, stool.Genus.prop)
otu_table(prune.dat_top5_Genus)

tax_table(stool.Genus.prop)
others_Genus<-taxa_names(stool.Genus.prop)
others_Genus
others_Genus<-others_Genus[!(others_Genus %in% top5_Genus)]
prune.others_Genus<- prune_taxa(others_Genus, stool.Genus.prop)

tabella_top_Genus<-psmelt(prune.dat_top5_Genus)
tabella_others_Genus<-psmelt(prune.others_Genus)
tabella_others_Genus$Genus<-"Others"
tabella_Genus<-rbind.data.frame(tabella_top_Genus,tabella_others_Genus)
Genus_hc<-tabella_top_Genus

write.csv(file = paste0("Results_CERBIOS/TOP_5_Genus/TOP5_Genus_INTERV_T1.csv"), row.names = F,
          cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top5_Genus),1,mean)), "Genus"= as.data.frame(tax_table(prune.dat_top5_Genus))[["Genus"]]))

fill_color_5<-c('Bifidobacterium'="coral",'Collinsella'="springgreen3",'Faecalibacterium'="gold3",'Blautia'="firebrick3",'Bacteroides'="deepskyblue3",'Others' ="darkslategray3")  
ggplot(data=tabella_Genus, aes(x=Sample, y=Abundance, fill=Genus)) + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=14) + 
  scale_fill_manual(values= fill_color_5) +
  #scale_fill_manual(values=c("#DF536B", "#61D04F" ,"#2297E6" ,"#28E2E5" ,"#CD0BBC", "#F5C710")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Subject", y="Relative abundance", title = paste0("Five most abundant Genus Intervention T1"))

ggsave(filename = paste0("Results_CERBIOS/TOP_5_Genus/Bar_plot_5_top_Genus_INTERV_T1.png"), width = 9, height = 5, dpi = 300)

###  CONTR T0
sample_data (data.genus.prop)
stool.Genus.prop<-subset_samples(data.genus.prop, Time=='T0'& Condition=='CONTR')
sample_data(stool.Genus.prop)

top5_Genus <- names(sort(taxa_sums(stool.Genus.prop), decreasing=TRUE))[1:5] 
prune.dat_top5_Genus <- prune_taxa(top5_Genus, stool.Genus.prop)
otu_table(prune.dat_top5_Genus)

tax_table(stool.Genus.prop)
others_Genus<-taxa_names(stool.Genus.prop)
others_Genus
others_Genus<-others_Genus[!(others_Genus %in% top5_Genus)]
prune.others_Genus<- prune_taxa(others_Genus, stool.Genus.prop)

tabella_top_Genus<-psmelt(prune.dat_top5_Genus)
tabella_others_Genus<-psmelt(prune.others_Genus)
tabella_others_Genus$Genus<-"Others"
tabella_Genus<-rbind.data.frame(tabella_top_Genus,tabella_others_Genus)
Genus_hc<-tabella_top_Genus

write.csv(file = paste0("Results_CERBIOS/TOP_5_Genus/TOP5_Genus_CONTR_T0.csv"), row.names = F,
          cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top5_Genus),1,mean)), "Genus"= as.data.frame(tax_table(prune.dat_top5_Genus))[["Genus"]]))

fill_color_5<-c('Bifidobacterium'="coral",'Subdoligranulum'="pink3",'Faecalibacterium'="gold3",'Blautia'="firebrick3",'Bacteroides'="deepskyblue3",'Others' ="darkslategray3") 
ggplot(data=tabella_Genus, aes(x=Sample, y=Abundance, fill=Genus)) + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=14) + 
  scale_fill_manual(values= fill_color_5) +
  #scale_fill_manual(values=c("#DF536B", "#61D04F" ,"#2297E6" ,"#28E2E5" ,"#CD0BBC", "#F5C710")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Subject", y="Relative abundance", title = paste0("Five most abundant Genus Control T0"))

ggsave(filename = paste0("Results_CERBIOS/TOP_5_Genus/Bar_plot_5_top_Genus_CONTR_T0.png"), width = 9, height = 5, dpi = 300)

###  CONTR T1
sample_data (data.genus.prop)
stool.Genus.prop<-subset_samples(data.genus.prop, Time=='T1'& Condition=='CONTR')
sample_data(stool.Genus.prop)

top5_Genus <- names(sort(taxa_sums(stool.Genus.prop), decreasing=TRUE))[1:5] 
prune.dat_top5_Genus <- prune_taxa(top5_Genus, stool.Genus.prop)
otu_table(prune.dat_top5_Genus)

tax_table(stool.Genus.prop)
others_Genus<-taxa_names(stool.Genus.prop)
others_Genus
others_Genus<-others_Genus[!(others_Genus %in% top5_Genus)]
prune.others_Genus<- prune_taxa(others_Genus, stool.Genus.prop)

tabella_top_Genus<-psmelt(prune.dat_top5_Genus)
tabella_others_Genus<-psmelt(prune.others_Genus)
tabella_others_Genus$Genus<-"Others"
tabella_Genus<-rbind.data.frame(tabella_top_Genus,tabella_others_Genus)
Genus_hc<-tabella_top_Genus

write.csv(file = paste0("Results_CERBIOS/TOP_5_Genus/TOP5_Genus_CONTR_T1.csv"), row.names = F,
          cbind.data.frame("mean"=as.numeric(apply(otu_table(prune.dat_top5_Genus),1,mean)), "Genus"= as.data.frame(tax_table(prune.dat_top5_Genus))[["Genus"]]))

fill_color_5<-c('Bifidobacterium'="coral",'Subdoligranulum'="pink3",'Faecalibacterium'="gold3",'Blautia'="firebrick3",'Bacteroides'="deepskyblue3",'Others' ="darkslategray3")  
ggplot(data=tabella_Genus, aes(x=Sample, y=Abundance, fill=Genus)) + 
  geom_bar(stat="identity", position="stack") + theme_classic(base_size=14) + 
  scale_fill_manual(values= fill_color_5) +
  #scale_fill_manual(values=c("#DF536B", "#61D04F" ,"#2297E6" ,"#28E2E5" ,"#CD0BBC", "#F5C710")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5), 
        legend.key.size = unit(0.4, "cm"),legend.text = element_text ( size = 12 )) + 
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) +
  labs(x="Subject", y="Relative abundance", title = paste0("Five most abundant Genus Control T1"))

ggsave(filename = paste0("Results_CERBIOS/TOP_5_Genus/Bar_plot_5_top_Genus_CNT_T1.png"), width = 9, height = 5, dpi = 300)

############### ALPHA DIVERSITY ############### 
sample_data(data)

#totale tabella indici alpha diversity
alpha<-estimate_richness(data)
vec<-row.names(alpha)
alpha$Sample<-vec
H<-alpha[,c(6,10)]
obs<-alpha[,c(1,10)]
# adding evenness
identical(H$Sample, obs$Sample) # TRUE
ev<-H
ev$Evenness<-(H$Shannon)/log((obs$Observed))
ev$Observed_richness<-obs$Observed
ev
vec1<-sample_data(data)$Condition
vec2<-sample_data(data)$Time
ev$Condition<-vec1
ev$Time<-vec2
ev
write.xlsx (ev,'Results_CERBIOS/Alpha_Div/Alpha_di_value.xlsx')

#T0
T0<-subset_samples(data, Time=='T0') 

{pAlpha<-plot_richness(T0, measures=c("Shannon", "Observed"), x="Condition", color="Condition")
  pAlpha
  # plot_richness( ) compute diversity like estimate_diversity( )
  H<-dplyr::filter(pAlpha$data, variable=="Shannon")
  obs<-dplyr::filter(pAlpha$data, variable=="Observed")
  # adding evenness
  identical(H$Sample, obs$Sample) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha +
  geom_boxplot(data=pAlpha$data, aes(x=Condition, y=value, color=NULL), alpha=0.1) +
  #geom_jitter(data=pAlpha$data, aes(x=Condition, y=value, shape=Sex), 
              #position=position_jitter(width=0.2), size=2) + # Add jittered points with different shapes for 'sex'
  theme_classic2(base_size = 11) + 
  scale_color_manual(values = c('CONTR'='lightgreen', 'INTERV'='cyan')) +
  labs(x="Condition",
       title="Alpha diversity Contr T0 - Interv T0") +
  guides(fill="none", color="none", shape=guide_legend(title="Sex")) + # Add legend for shape
  stat_compare_means(aes(group = Condition), method = "wilcox.test", paired=TRUE, 
                     label="p.format", label.x=1.5, size=3.5, label.y.npc="top", 
                     vjust=-0.5, hjust=0.45) +
  theme(panel.grid.major.y = element_line(size=0.4),
        panel.grid.minor.y = element_line(size=0.2), 
        axis.text.x = element_text(angle=28, vjust=1, hjust=1, size=10),
        axis.title.x = element_text(vjust = -1)
  ) +
  scale_color_manual(values = c('CONTR'='lightgreen', 'INTERV'='cyan')) 

ggsave(file=paste("Results_CERBIOS/Alpha_Div/Alpha diversity T0.png", sep=""),  width = 8, height = 6, dpi=300)  

#T1
T1<-subset_samples(data, Time=='T1') 
sample_data(T1)

{pAlpha<-plot_richness(T1, measures=c("Shannon", "Observed"), x="Condition", color="Condition")
  pAlpha
  # plot_richness( ) compute diversity like estimate_diversity( )
  H<-dplyr::filter(pAlpha$data, variable=="Shannon")
  obs<-dplyr::filter(pAlpha$data, variable=="Observed")
  # adding evenness
  identical(H$Sample, obs$Sample) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha +
  geom_boxplot(data=pAlpha$data, aes(x=Condition, y=value, color=NULL), alpha=0.1) +
  #geom_jitter(data=pAlpha$data, aes(x=Condition, y=value, shape=Sex), 
  #position=position_jitter(width=0.2), size=2) + # Add jittered points with different shapes for 'sex'
  theme_classic2(base_size = 11) + 
  scale_color_manual(values = c('CONTR'='green3', 'INTERV'='cyan')) +
  labs(x="Condition",
       title="Alpha diversity Contr T1 - Interv T1") +
  guides(fill="none", color="none", shape=guide_legend(title="Sex")) + # Add legend for shape
  stat_compare_means(aes(group = Condition), method = "wilcox.test", paired=TRUE, 
                     label="p.format", label.x=1.5, size=3.5, label.y.npc="top", 
                     vjust=-0.5, hjust=0.45) +
  theme(panel.grid.major.y = element_line(size=0.4),
        panel.grid.minor.y = element_line(size=0.2), 
        axis.text.x = element_text(angle=28, vjust=1, hjust=1, size=10),
        axis.title.x = element_text(vjust = -1)
  ) +
  scale_color_manual(values = c('CONTR'='green3', 'INTERV'='cyan')) 
ggsave(file=paste("Results_CERBIOS/Alpha_Div/Alpha diversity T1.png", sep=""),  width = 8, height = 6, dpi=300)

#CONTR
CONTR<-subset_samples(data, Condition=='CONTR') 
sample_data(CONTR)

{pAlpha<-plot_richness(CONTR, measures=c("Shannon", "Observed"), x="Time", color="Time")
  pAlpha
  # plot_richness( ) compute diversity like estimate_diversity( )
  H<-dplyr::filter(pAlpha$data, variable=="Shannon")
  obs<-dplyr::filter(pAlpha$data, variable=="Observed")
  # adding evenness
  identical(H$Sample, obs$Sample) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha +
  geom_boxplot(data=pAlpha$data, aes(x=Time, y=value, color=NULL), alpha=0.1) +
  #geom_jitter(data=pAlpha$data, aes(x=Time, y=value, shape=Sex), 
  #position=position_jitter(width=0.2), size=2) + # Add jittered points with different shapes for 'sex'
  theme_classic2(base_size = 11) + 
  scale_color_manual(values = c('T0'='lightgreen', 'T1'='green3')) +
  labs(x="Time",
       title="Alpha diversity Control T0 vs T1") +
  guides(fill="none", color="none", shape=guide_legend(title="Sex")) + # Add legend for shape
  stat_compare_means(aes(group = Time), method = "wilcox.test", paired=TRUE, 
                     label="p.format", label.x=1.5, size=3.5, label.y.npc="top", 
                     vjust=-0.5, hjust=0.45) +
  theme(panel.grid.major.y = element_line(size=0.4),
        panel.grid.minor.y = element_line(size=0.2), 
        axis.text.x = element_text(angle=28, vjust=1, hjust=1, size=10),
        axis.title.x = element_text(vjust = -1)
  ) +
  scale_color_manual(values = c('T0'='lightgreen', 'T1'='green3')) 

ggsave(file=paste("Results_CERBIOS/Alpha_Div/Alpha diversity Control.png", sep=""),  width = 8, height = 6, dpi=300)

######## *** MODIFICA PER REVISORI *** ######## *** 
# Filtra per una metrica, es: Shannon
df_shannon <- pAlpha$data %>%
  filter(variable == "Shannon")
df_shannon$ID

# Dividi la colonna ID in due: Subject e Time
df_shannon <- df_shannon %>%
  separate(ID, into = c("Subject", "Time"), sep = " ")

# Plot
plot_shannon<-ggplot(df_shannon, aes(x = Time, y = value)) +
  geom_line(aes(group = Subject), color = "gray70", alpha = 0.6) +
  geom_boxplot(aes(group = Time), width = 0.4, outlier.shape = NA, alpha = 0.2) +
  geom_point(aes(color = Time), size = 2) +
  theme_classic() +
  labs(title = "Shannon",
       x = "Time", y = "") +
  stat_compare_means(paired = TRUE, method = "wilcox.test",
                     label = "p.format", label.y = max(df_shannon$value) + 0.2, size = 4) +
  scale_color_manual(values = c('T0' = 'lightgreen', 'T1' = 'green3'))

## --- Evenness ---
df_evenness <- pAlpha$data %>%
  filter(variable == "Evenness") %>%
  separate(ID, into = c("Subject", "Time"), sep = " ")

plot_evenness <- ggplot(df_evenness, aes(x = Time, y = value)) +
  geom_line(aes(group = Subject), color = "gray70", alpha = 0.6) +
  geom_boxplot(aes(group = Time), width = 0.4, outlier.shape = NA, alpha = 0.2) +
  geom_point(aes(color = Time), size = 2) +
  theme_classic() +
  labs(title = "Evenness", x = "Time", y = "") +
  stat_compare_means(paired = TRUE, method = "wilcox.test",
                     label = "p.format", label.y = max(df_evenness$value) + 0.05, size = 4) +
  scale_color_manual(values = c('T0' = 'lightgreen', 'T1' = 'green3'))

## --- Richness ---
df_richness <- pAlpha$data %>%
  filter(variable == "Observed richness") %>%  # oppure "Richness", dipende come si chiama
  separate(ID, into = c("Subject", "Time"), sep = " ")

plot_richness <- ggplot(df_richness, aes(x = Time, y = value)) +
  geom_line(aes(group = Subject), color = "gray70", alpha = 0.6) +
  geom_boxplot(aes(group = Time), width = 0.4, outlier.shape = NA, alpha = 0.2) +
  geom_point(aes(color = Time), size = 2) +
  theme_classic() +
  labs(title = "Observed Richness", x = "Time", y = "Alpha diversity measure") +
  stat_compare_means(paired = TRUE, method = "wilcox.test",
                     label = "p.format", label.y = max(df_richness$value) + 2, size = 4) +
  scale_color_manual(values = c('T0' = 'lightgreen', 'T1' = 'green3'))

# Carica patchwork
library(patchwork)

combined_plot <- (plot_richness | plot_shannon | plot_evenness) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

combined_plot

ggsave("Alpha_div_matched_control.png", plot = combined_plot, width = 12, height = 10, dpi = 300)


#INTERV
INTERV<-subset_samples(data, Condition=='INTERV') 
sample_data(INTERV)

{pAlpha<-plot_richness(INTERV, measures=c("Shannon", "Observed"), x="Time", color="Time")
  pAlpha
  # plot_richness( ) compute diversity like estimate_diversity( )
  H<-dplyr::filter(pAlpha$data, variable=="Shannon")
  obs<-dplyr::filter(pAlpha$data, variable=="Observed")
  # adding evenness
  identical(H$Sample, obs$Sample) # TRUE
  ev<-H
  ev$value<-(H$value)/log((obs$value))
  ev$variable<-rep("Evenness")
  New_data<-rbind.data.frame(obs,H,ev)
  pAlpha$data<-New_data
  pAlpha$data$variable<-gsub("Observed","Observed richness",pAlpha$data$variable)
  pAlpha$data$variable<-factor(pAlpha$data$variable,levels = c("Observed richness","Shannon","Evenness"))
}
pAlpha +
  geom_boxplot(data=pAlpha$data, aes(x=Time, y=value, color=NULL), alpha=0.1) +
  #geom_jitter(data=pAlpha$data, aes(x=Time, y=value, shape=Sex), 
  #position=position_jitter(width=0.2), size=2) + # Add jittered points with different shapes for 'sex'
  theme_classic2(base_size = 11) + 
  scale_color_manual(values = c('T0'='cyan', 'T1'='lightblue')) +
  labs(x="Time",
       title="Alpha diversity Intervention T0 vs T1") +
  guides(fill="none", color="none", shape=guide_legend(title="Sex")) + # Add legend for shape
  stat_compare_means(aes(group = Time), method = "wilcox.test", paired=TRUE, 
                     label="p.format", label.x=1.5, size=3.5, label.y.npc="top", 
                     vjust=-0.5, hjust=0.45) +
  theme(panel.grid.major.y = element_line(size=0.4),
        panel.grid.minor.y = element_line(size=0.2), 
        axis.text.x = element_text(angle=28, vjust=1, hjust=1, size=10),
        axis.title.x = element_text(vjust = -1)
  ) +
  scale_color_manual(values = c('T0'='cyan', 'T1'='lightblue')) 

ggsave(file=paste("Results_CERBIOS/Alpha_Div/Alpha diversity Intervention.png", sep=""),  width = 8, height = 6, dpi=300)

######## *** MODIFICA PER REVISORI *** ######## *** 
# Filtra per una metrica, es: Shannon
df_shannon <- pAlpha$data %>%
  filter(variable == "Shannon")
df_shannon$ID

# Dividi la colonna ID in due: Subject e Time
df_shannon <- df_shannon %>%
  separate(ID, into = c("Subject", "Time"), sep = " ")

# Plot
plot_shannon<-ggplot(df_shannon, aes(x = Time, y = value)) +
  geom_line(aes(group = Subject), color = "gray70", alpha = 0.6) +
  geom_boxplot(aes(group = Time), width = 0.4, outlier.shape = NA, alpha = 0.2) +
  geom_point(aes(color = Time), size = 2) +
  theme_classic() +
  labs(title = "Shannon",
       x = "Time", y = "") +
  stat_compare_means(paired = TRUE, method = "wilcox.test",
                     label = "p.format", label.y = max(df_shannon$value) + 0.2, size = 4) +
  scale_color_manual(values = c('T0' = 'cyan', 'T1' = 'lightblue'))

## --- Evenness ---
df_evenness <- pAlpha$data %>%
  filter(variable == "Evenness") %>%
  separate(ID, into = c("Subject", "Time"), sep = " ")

plot_evenness <- ggplot(df_evenness, aes(x = Time, y = value)) +
  geom_line(aes(group = Subject), color = "gray70", alpha = 0.6) +
  geom_boxplot(aes(group = Time), width = 0.4, outlier.shape = NA, alpha = 0.2) +
  geom_point(aes(color = Time), size = 2) +
  theme_classic() +
  labs(title = "Evenness", x = "Time", y = "") +
  stat_compare_means(paired = TRUE, method = "wilcox.test",
                     label = "p.format", label.y = max(df_evenness$value) + 0.05, size = 4) +
  scale_color_manual(values =c('T0' = 'cyan', 'T1' = 'lightblue'))

## --- Richness ---
df_richness <- pAlpha$data %>%
  filter(variable == "Observed richness") %>%  # oppure "Richness", dipende come si chiama
  separate(ID, into = c("Subject", "Time"), sep = " ")

plot_richness <- ggplot(df_richness, aes(x = Time, y = value)) +
  geom_line(aes(group = Subject), color = "gray70", alpha = 0.6) +
  geom_boxplot(aes(group = Time), width = 0.4, outlier.shape = NA, alpha = 0.2) +
  geom_point(aes(color = Time), size = 2) +
  theme_classic() +
  labs(title = "Observed Richness", x = "Time", y = "Alpha diversity measure") +
  stat_compare_means(paired = TRUE, method = "wilcox.test",
                     label = "p.format", label.y = max(df_richness$value) + 2, size = 4) +
  scale_color_manual(values = c('T0' = 'cyan', 'T1' = 'lightblue'))

# Carica patchwork
library(patchwork)

combined_plot <- (plot_richness | plot_shannon | plot_evenness) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

combined_plot

ggsave("Alpha_div_matched_intervention.png", plot = combined_plot, width = 12, height = 10, dpi = 300)


########### BETA DIVERSITY & PERMANOVA ###############
# T0
b.stool<-subset_samples(data.genus.prop , Time=='T0')
sample_data(b.stool)
data.sqrt_prop<-transform_sample_counts(b.stool, sqrt)
DistBC <- phyloseq::distance(data.sqrt_prop, method = "euclidean") 
ordBC <- ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)

PERM<-adonis2(DistBC ~ sample_data(b.stool)[['Condition']], data=as(sample_data(data.sqrt_prop),"data.frame"), permutations = 9999)
DISP<-permutest(betadisper(DistBC, sample_data(data.sqrt_prop)[['Condition']]))
temp<-cbind(PERM$`Pr(>F)`[1],DISP$tab$`Pr(>F)`[1])

color3= c('CONTR'='green3','INTERV'='blue2')
plot_ordination(data.sqrt_prop, ordBC, color = 'Condition') +
  geom_point(size=3) + theme_bw() + stat_ellipse()+
  labs(title=paste("PCoA on genera with Hellinger distance \n(euclidean of sqrt of proportions) T0"), 
       subtitle = paste0("Pr(>F): ",PERM$`Pr(>F)`[1]),color3= 'Condition') +
  scale_color_manual(values = color3)

ggsave(file=paste("Results_CERBIOS/Beta_Div/PCoA Hellinger Beta diversity T0.png", sep=""),  width = 8, height = 6, dpi=300)

# T1
b.stool<-subset_samples(data.genus.prop , Time=='T1')
sample_data(b.stool)
data.sqrt_prop<-transform_sample_counts(b.stool, sqrt)
DistBC <- phyloseq::distance(data.sqrt_prop, method = "euclidean") 
ordBC <- ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)

PERM<-adonis2(DistBC ~ sample_data(b.stool)[['Condition']], data=as(sample_data(data.sqrt_prop),"data.frame"), permutations = 9999)
DISP<-permutest(betadisper(DistBC, sample_data(data.sqrt_prop)[['Condition']]))
temp<-cbind(PERM$`Pr(>F)`[1],DISP$tab$`Pr(>F)`[1])

color3= c('CONTR'='green3','INTERV'='blue2')
plot_ordination(data.sqrt_prop, ordBC, color = 'Condition') +
  geom_point(size=3) + theme_bw() + stat_ellipse()+
  labs(title=paste("PCoA on genera with Hellinger distance \n(euclidean of sqrt of proportions) T1"), 
       subtitle = paste0("Pr(>F): ",PERM$`Pr(>F)`[1]),color3= 'Condition') +
  scale_color_manual(values = color3)

ggsave(file=paste("Results_CERBIOS/Beta_Div/PCoA Hellinger Beta diversity T1.png", sep=""),  width = 8, height = 6, dpi=300)

# CONTROL
b.stool<-subset_samples(data.genus.prop , Condition=='CONTR')
sample_data(b.stool)
data.sqrt_prop<-transform_sample_counts(b.stool, sqrt)
DistBC <- phyloseq::distance(data.sqrt_prop, method = "euclidean") 
ordBC <- ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)

PERM<-adonis2(DistBC ~ sample_data(b.stool)[['Time']], data=as(sample_data(data.sqrt_prop),"data.frame"), permutations = 9999)
DISP<-permutest(betadisper(DistBC, sample_data(data.sqrt_prop)[['Time']]))
temp<-cbind(PERM$`Pr(>F)`[1],DISP$tab$`Pr(>F)`[1])

color3= c('T0'='lightgreen','T1'='green3')
plot_ordination(data.sqrt_prop, ordBC, color = 'Time') +
  geom_point(size=3) + theme_bw() + stat_ellipse()+
  labs(title=paste("PCoA on genera with Hellinger distance \n(euclidean of sqrt of proportions) Control"), 
       subtitle = paste0("Pr(>F): ",PERM$`Pr(>F)`[1]),color3= 'Time') +
  scale_color_manual(values = color3)

ggsave(file=paste("Results_CERBIOS/Beta_Div/PCoA Hellinger Beta diversity CNT.png", sep=""),  width = 8, height = 6, dpi=300)

######### *** BETA DIVERSITY REVISORI *** #########
# Ordination
DistBC <- phyloseq::distance(data.sqrt_prop, method = "euclidean") 
ordBC <- ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)

# Estrai coordinate
scores_df <- as.data.frame(ordBC$vectors[,1:2])

scores_df$SampleID <- rownames(scores_df)

metadata <- data.frame(sample_data(data.sqrt_prop))
metadata$SampleID <- rownames(metadata)

plot_data <- merge(scores_df, metadata, by = "SampleID")

plot_data <- plot_data %>%
  separate(ID, into = c("Subject", "Time"), sep = " ", remove = FALSE)

ggplot(plot_data, aes(x = Axis.1, y = Axis.2, color = Time)) +
  geom_point(size=3) +
  geom_path(aes(group = Subject), color = "gray50", alpha = 0.7) + 
  stat_ellipse() +
  theme_bw() +
  labs(title="PCoA on genera with Hellinger distance",
       subtitle = paste0("PERMANOVA p = ", round(PERM$`Pr(>F)`[1], 4))) +
  scale_color_manual(values = c('T0'='lightgreen','T1'='green3'))
ggsave(file="PCoA Hellinger Beta diversity CNT matched.png",  width = 8, height = 6, dpi=300)


# INTERV
b.stool<-subset_samples(data.genus.prop , Condition=='INTERV')
sample_data(b.stool)
data.sqrt_prop<-transform_sample_counts(b.stool, sqrt)
DistBC <- phyloseq::distance(data.sqrt_prop, method = "euclidean") 
ordBC <- ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)

PERM<-adonis2(DistBC ~ sample_data(b.stool)[['Time']], data=as(sample_data(data.sqrt_prop),"data.frame"), permutations = 9999)
DISP<-permutest(betadisper(DistBC, sample_data(data.sqrt_prop)[['Time']]))
temp<-cbind(PERM$`Pr(>F)`[1],DISP$tab$`Pr(>F)`[1])

color3= c('T0'='cyan','T1'='blue2')
plot_ordination(data.sqrt_prop, ordBC, color = 'Time') +
  geom_point(size=3) + theme_bw() + stat_ellipse()+
  labs(title=paste("PCoA on genera with Hellinger distance \n(euclidean of sqrt of proportions) Intervention"), 
       subtitle = paste0("Pr(>F): ",PERM$`Pr(>F)`[1]),color3= 'Time') +
  scale_color_manual(values = color3)

ggsave(file=paste("Results_CERBIOS/Beta_Div/PCoA Hellinger Beta diversity INTERV.png", sep=""),  width = 8, height = 6, dpi=300)

######### *** BETA DIVERSITY REVISORI *** #########
# Ordination
DistBC <- phyloseq::distance(data.sqrt_prop, method = "euclidean") 
ordBC <- ordinate(data.sqrt_prop, method = "PCoA", distance = DistBC)

# Estrai coordinate
scores_df <- as.data.frame(ordBC$vectors[,1:2])

scores_df$SampleID <- rownames(scores_df)

metadata <- data.frame(sample_data(data.sqrt_prop))
metadata$SampleID <- rownames(metadata)

plot_data <- merge(scores_df, metadata, by = "SampleID")

plot_data <- plot_data %>%
  separate(ID, into = c("Subject", "Time"), sep = " ", remove = FALSE)

ggplot(plot_data, aes(x = Axis.1, y = Axis.2, color = Time)) +
  geom_point(size=3) +
  geom_path(aes(group = Subject), color = "gray50", alpha = 0.7) + 
  stat_ellipse() +
  theme_bw() +
  labs(title="PCoA on genera with Hellinger distance",
       subtitle = paste0("PERMANOVA p = ", round(PERM$`Pr(>F)`[1], 4))) +
  scale_color_manual(values = c('T0'='cyan','T1'='blue2'))
ggsave(file="PCoA Hellinger Beta diversity INT matched.png",  width = 8, height = 6, dpi=300)

############# DESEQ2 ############# 
# T1
T1 <- subset_samples(data, Time=='T1')
sample_data(T1)

tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus")
all_results <- list()
tax_table_data <- tax_table(T1)
meta <- as.data.frame(sample_data(T1))
tax_table_data <- as.data.frame(tax_table(T1))
all_results <- list()

for (tax_level in tax_levels) {
  tax_data <- tax_glom(T1, taxrank=tax_level, NArm=FALSE)
  tax_data_df <- as.data.frame(otu_table(tax_data))
  dds <- DESeqDataSetFromMatrix(countData = tax_data_df,
                                colData = meta,
                                design = ~ Condition)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  res <- res[order(res$padj, na.last=NA), ] 
  res <- res[(res$padj < 0.05) & abs(res$log2FoldChange) > 1, ]
  res <- as.data.frame(res)
  res_filtered <- subset(res, baseMean > 100)
  if (nrow(res_filtered) == 0) {
    cat("Nessun risultato significativo per il livello", tax_level, "\n")
    next
  }
  res_filtered$OTU <- rownames(res_filtered)
  res_filtered$OTU <- gsub(".*\\.", "", res_filtered$OTU)
  res_filtered$TaxLevel <- tax_level
  res_filtered$Taxonomy <- tax_table_data[res_filtered$OTU, tax_level]
  rownames(res_filtered) <- NULL
  all_results[[tax_level]] <- res_filtered
}

final_results <- do.call(rbind, all_results)
final_results
write.xlsx(final_results, 'Results_CERBIOS/Deseq2/DE_res_T1.xlsx')

#Genus
final_results
Dati_box<-final_results
Dati_box_g<-Dati_box[Dati_box$TaxLevel=="Genus","OTU"] 
Dati_box_g 
data.order<- tax_glom(T1, taxrank = "Genus", NArm = F)
data.order <- transform_sample_counts(data.order, function(x) x/sum(x)*100)
Dati_box_g<-prune_taxa(Dati_box_g, data.order)
Dati_box_g
Dati_box_g<-psmelt(Dati_box_g)
Dati_box_g

#
tabella_g <- Dati_box_g
tabella_g$Taxa<-"Genus"
tabella_g[,c("Phylum","Order","Class","Family")]<-NULL
colnames(tabella_g)[colnames(tabella_g)=="Genus"]<-"Bacteria"

tabella_g
ggplot(tabella_g, aes(x= Bacteria, y=Abundance, fill=Condition)) + 
  facet_grid(~factor(Taxa,levels = "Genus"), scales = "free_x", space="free") +
  scale_fill_manual(values=c("INTERV"="cyan","CONTR"="green3")) +
  scale_color_manual(values=c("INTERV"="cyan","CONTR"="green3")) +
  geom_boxplot(width=0.9, size= 0.4, alpha= 0.2, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.9, jitter.width = 0.4),
             aes(color=Condition), size= 0.8, alpha= 0.65) +  
  theme_classic2(base_size = 11) + 
  theme(strip.text.x=element_text(size=11.7,colour="black"),
        legend.margin=margin(-25, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=10), 
        axis.text.y = element_text(size=8.8),
        plot.title= element_text(size=15) ,
        legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=14),
        panel.grid.major.y = element_line(size=0.12, color="gray")
  ) +
  guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa T1", y="Percentual Abundance", fill="Condition", x="") + 
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(0.1, 0.5 ,1, seq(2,20,2),seq(10,max(tabella_g$Abundance)+4,3)))
ggsave(filename = "Results_CERBIOS/Deseq2/DE_T1.png", width = 7.5, height = 6, dpi=300)

# CONTR
CONTR <- subset_samples(data, Condition=='CONTR')
sample_data(CONTR)

tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus")
all_results <- list()
tax_table_data <- tax_table(CONTR)
meta <- as.data.frame(sample_data(CONTR))
tax_table_data <- as.data.frame(tax_table(CONTR))
all_results <- list()

for (tax_level in tax_levels) {
  tax_data <- tax_glom(CONTR, taxrank=tax_level, NArm=FALSE)
  tax_data_df <- as.data.frame(otu_table(tax_data))
  dds <- DESeqDataSetFromMatrix(countData = tax_data_df,
                                colData = meta,
                                design = ~ Time)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  res <- res[order(res$padj, na.last=NA), ] 
  res <- res[(res$padj < 0.05) & abs(res$log2FoldChange) > 1, ]
  res <- as.data.frame(res)
  res_filtered <- subset(res, baseMean > 100)
  if (nrow(res_filtered) == 0) {
    cat("Nessun risultato significativo per il livello", tax_level, "\n")
    next
  }
  res_filtered$OTU <- rownames(res_filtered)
  res_filtered$OTU <- gsub(".*\\.", "", res_filtered$OTU)
  res_filtered$TaxLevel <- tax_level
  res_filtered$Taxonomy <- tax_table_data[res_filtered$OTU, tax_level]
  rownames(res_filtered) <- NULL
  all_results[[tax_level]] <- res_filtered
}

final_results <- do.call(rbind, all_results)
final_results #NULL

# INTERV
INTERV <- subset_samples(data, Condition=='INTERV')
sample_data(INTERV)

tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus")
all_results <- list()
tax_table_data <- tax_table(INTERV)
meta <- as.data.frame(sample_data(INTERV))
tax_table_data <- as.data.frame(tax_table(INTERV))
all_results <- list()

for (tax_level in tax_levels) {
  tax_data <- tax_glom(INTERV, taxrank=tax_level, NArm=FALSE)
  tax_data_df <- as.data.frame(otu_table(tax_data))
  dds <- DESeqDataSetFromMatrix(countData = tax_data_df,
                                colData = meta,
                                design = ~ Time)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  res <- res[order(res$padj, na.last=NA), ] 
  res <- res[(res$padj < 0.05) & abs(res$log2FoldChange) > 1, ]
  res <- as.data.frame(res)
  res_filtered <- subset(res, baseMean > 100)
  if (nrow(res_filtered) == 0) {
    cat("Nessun risultato significativo per il livello", tax_level, "\n")
    next
  }
  res_filtered$OTU <- rownames(res_filtered)
  res_filtered$OTU <- gsub(".*\\.", "", res_filtered$OTU)
  res_filtered$TaxLevel <- tax_level
  res_filtered$Taxonomy <- tax_table_data[res_filtered$OTU, tax_level]
  rownames(res_filtered) <- NULL
  all_results[[tax_level]] <- res_filtered
}

final_results <- do.call(rbind, all_results)
final_results 

write.xlsx(final_results, 'Results_CERBIOS/Deseq2/DE_res_T1.xlsx')

#Order
final_results
Dati_box<-final_results
Dati_box_o<-Dati_box[Dati_box$TaxLevel=="Order","OTU"] 
Dati_box_o 
data.order<- tax_glom(INTERV, taxrank = "Order", NArm = F)
data.order <- transform_sample_counts(data.order, function(x) x/sum(x)*100)
Dati_box_o<-prune_taxa(Dati_box_o, data.order)
Dati_box_o
Dati_box_o<-psmelt(Dati_box_o)
Dati_box_o

#
head(colnames(tax_table(data)))
tabella_o <- Dati_box_o
tabella_o$Taxa<-"Order"
tabella_o[,c("Phylum","Class")]<-NULL
colnames(tabella_o)[colnames(tabella_o)=="Order"]<-"Bacteria"

tabella_o
ggplot(tabella_o, aes(x= Bacteria, y=Abundance, fill=Time)) + 
  facet_grid(~factor(Taxa,levels = "Order"), scales = "free_x", space="free") +
  scale_fill_manual(values=c("T0"="cyan","T1"="blue2")) +
  scale_color_manual(values=c("T0"="cyan","T1"="blue2")) +
  geom_boxplot(width=0.9, size= 0.4, alpha= 0.2, outlier.alpha = 0) +
  geom_point(position = position_jitterdodge(seed = 1994, dodge.width = 0.9, jitter.width = 0.4),
             aes(color=Time), size= 0.8, alpha= 0.65) +  
  theme_classic2(base_size = 11) + 
  theme(strip.text.x=element_text(size=11.7,colour="black"),
        legend.margin=margin(-25, 0, 0, 0), legend.position="bottom" ,
        axis.text.x = element_text(angle = 22, vjust=1, hjust=1, size=10), 
        axis.text.y = element_text(size=8.8),
        plot.title= element_text(size=15) ,
        legend.key.size=unit(0.8,"cm"), 
        legend.text=element_text(size=14),
        panel.grid.major.y = element_line(size=0.12, color="gray")
  ) +
  guides( fill=guide_legend(nrow=1) ) +
  labs(title= "Differently abundant Taxa Intervention T0 - T1", y="Percentual Abundance", fill="Time", x="") + 
  scale_x_discrete(expand=c(-0.2, 1)) + 
  scale_y_sqrt(breaks=c(10,20,30, seq(20,100,20),seq(10,max(tabella_o$Abundance)+15,3)))
ggsave(filename = "Results_CERBIOS/Deseq2/DE_INTERV.png", width = 7.5, height = 6, dpi=300)


########## Acidi Grassi ########## 
###################### PREPARING DATA #########################
{library(Hmisc)
  library(ggpubr)
  library(ggplot2)
  library(vegan)
  library(pairwiseAdonis)
  library(mixOmics)
  library(ecodist)
  library(reshape2)
  library(readxl)
}


options(scipen=100)
setwd('SARA/CERBIOS')
FFA <- read_excel("FFA.xlsx", sheet=1) 
colnames(FFA)
#rimuovo gli stessi sample rimossi nel microbiota
FFA<-subset.data.frame(FFA, Sample!='BC3')
FFA<-subset.data.frame(FFA, Sample!='BC19')
FFA<-subset.data.frame(FFA, Sample!='BC39')
FFA<-subset.data.frame(FFA, Sample!='BC55')

FFA[! colnames(FFA) %in% c("Sample","Condition",'Time')]<-apply(FFA[!colnames(FFA)%in% c("Sample","Condition","Time")], MARGIN = 2, as.numeric)
FFA$Condition<-factor(FFA$Condition, levels = c("INTERV","CONTR"))
FFA$Time<-factor(FFA$Time, levels = c("T0","T1"))

SCFA<-FFA[,1:10] # from acetic to valeric 
MCFA<-FFA[,c(1:3,11:20)] # to dodecanoic
LCFA<-FFA[,c(1:3,21:23)] # to Octadecanoic, before the raw total columns
colnames(LCFA)

######################### NORMALIZATIONS (stool, già in % quindi non serve) ##############################
SCFA_raw<-as.data.frame(SCFA)
MCFA_raw<-as.data.frame(MCFA)
LCFA_raw<- as.data.frame(LCFA)

####################### MANN WITHNEY TESTs ####################
setwd('Results_CERBIOS')
dir.create("FFA")

##### FFA raw abundances CONTR
head(FFA)
FFA_cnt<-subset.data.frame(FFA, Condition=='CONTR')
SCFA<-FFA_cnt[,1:10] # from acetic to valeric 
MCFA<-FFA_cnt[,c(1:3,11:20)] # to dodecanoic
LCFA<-FFA_cnt[,c(1:3,21:23)] # to Octadecanoic, before the raw total columns
colnames(LCFA)

SCFA_raw<-as.data.frame(SCFA)
MCFA_raw<-as.data.frame(MCFA)
LCFA_raw<- as.data.frame(LCFA)

list<-list("SCFA"=SCFA_raw, "MCFA"=MCFA_raw, "LCFA"=LCFA_raw)
names(list)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("Sample","Condition","Time")]) {
    test<-wilcox.test(temp[[x]]~temp$Time)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  write.csv2(results, file=paste0("FFA/",names(list)[y],"_raw_abund_Mann_Whit_CNT.csv"), row.names = F)
  
  table <- melt(data = temp, id.vars = c("Sample", "Condition","Time"))
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Time)) +
    scale_fill_manual(values = c("T0"="#32CD32","T1"="#FF4500")) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  ggsave(filename = paste0("FFA/",names(list)[y],"_raw_abundance_CNT.png"), height=3, width = 6, dpi=300)
  
  if(y==2){ # if MCFA
    table_MCFA<-table }
  
  if(y==3){ # if LCFA
    table_LCFA<-table }
}


# specific plot for SCFA
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Acetic", "Propionic","Butyric","Valeric"), "", "")
signif_y<-tapply(table_SCFA$value, table_SCFA$variable, max )
table_SCFA$signif_y<-as.numeric(signif_y[table_SCFA$variable]) + 10 # plus 5 to put the sign higher

ggplot(data = table_SCFA, mapping = aes(x = variable, y = value, fill = Time)) +
  scale_fill_manual(values = c("T0" = alpha("lightgreen", 0.5), "T1" = alpha("green3", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label = signif, y = signif_y, group = variable), 
            color = "darkgray", size = 8) +
  theme(axis.text.x = element_text(angle = -20, size = 9.5, vjust = 1, hjust = 0.05),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 11),
        panel.grid.minor.y = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.35),
        title = element_text(size = 10.5)) +
  scale_y_log10(breaks = c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50, seq(100, 500, 100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width = 0.65, size = 0.5, outlier.size = 1) + 
  labs(title = "SCFA concentration", x = "", fill = "", y = "Percentage (%)")
ggsave(filename = "FFA/SCFA_zoom_plot_CNT.png", height=4.5, width = 7, dpi=300)


# specific plot for MCFA
head(table_MCFA, n=2)
unique(table_MCFA$variable)
table_MCFA$signif<-ifelse(table_MCFA$variable %in% c("Octanoic","Decanoic", "Dodecanoic"), "", "")   # DA CAMBIARE!
signif_y<-tapply(table_MCFA$value, table_MCFA$variable, max )
table_MCFA$signif_y<-as.numeric(signif_y[table_MCFA$variable]) + 2 # plus 5 to put the sign higher

ggplot(data=table_MCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("lightgreen", 0.5), "T1" = alpha("green3", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="MCFA concentration", x="", fill="", y="Concentration (µM/L)")
ggsave(filename = "FFA/MCFA_zoom_plot_CNT.png", height=4.5, width = 7, dpi=300)


# specific plot for LCFA
head(table_LCFA, n=2)
table_LCFA$signif<-ifelse(table_LCFA$variable %in% c("Hexadecanoic","Octadecanoic"), "", "")   # DA CAMBIARE!
signif_y<-tapply(table_LCFA$value, table_LCFA$variable, max )
table_LCFA$signif_y<-as.numeric(signif_y[table_LCFA$variable]) + 10 # plus 5 to put the sign higher

# View(table_LCFA)
ggplot(data=table_LCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("lightgreen", 0.5), "T1" = alpha("green3", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 700,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="LCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA/LCFA_zoom_plot_CNT.png", height=4.5, width = 7, dpi=300)

########### PCoA ########### 
list_PCoA <- list("SCFA" = SCFA_raw, "MCFA" = MCFA_raw, "LCFA" = LCFA_raw)
LCFA_raw$Time

p_values <- NULL
for (x in 1:length(list_PCoA)) {
  rm(data, coord)
  data <- as.data.frame(list_PCoA[[x]])
  row.names(data) <- data$Sample
  
  if (x == 3) {
    data <- data[!is.na(data$Octadecanoic),] # rimuovere NA in LCFA
  }
  
  dist <- vegdist(data[, !colnames(data) %in% c("Sample", "Condition","Time")], method = "bray") # calcolare dissimilarità Bray-Curtis
  coord <- cmdscale(dist, k = 2, eig = TRUE) # PCoA
  
  # Convertire le coordinate in un dataframe e aggiungere Subject e Time
  coord_df <- as.data.frame(coord$points)
  colnames(coord_df) <- c("PC1", "PC2")
  coord_df$Subject <- row.names(coord_df)
  coord_df$Time <- data$Time
  
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Time)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("T0" = alpha("lightgreen", 0.5), "T1" = alpha("green3", 0.5))) +
    stat_ellipse(size = 0.2) +
    geom_text(aes(label = Subject), size = 2, color = "black") +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of Control subjects\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("PCoA_BRAY_", names(list_PCoA)[x], "CNT.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )     
  
  # again but with no names
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Time)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("T0" = alpha("lightgreen", 0.5), "T1" = alpha("green3", 0.5))) +
    stat_ellipse(size = 0.2) +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of Control subjects\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("PCoA_BRAY_", names(list_PCoA)[x], "CNT_noname.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )                        
  
  # checking for significant dispersion
  suppressWarnings(rm(p_value_diver, p_value_disp))
  p_value_diver<-adonis(dist ~ data$Time, data = data, permutations = 9999)
  p_value_disp<-permutest(betadisper(dist, data$Time), permutations = 9999)
  p_value_disp<-as.data.frame(p_value_disp$tab[1,])
  colnames(p_value_disp)<-colnames(p_value_diver$aov.tab[1,])
  
  new_values<-rbind(p_value_diver$aov.tab[1,],p_value_disp)
  row.names(new_values)<-c(paste("Beta_diversity",names(list_PCoA)[x], sep = "_"),
                           paste("Beta_dispersion",names(list_PCoA)[x], sep = "_"))
  p_values<-rbind(p_values,new_values)
}  

p_values$adjusted_BH_p_values<-p.adjust(p_values$`Pr(>F)`, method = "BH")
p_values
write.csv2(p_values, file="p_values_BRAY_CNT.csv", quote = F,row.names = T)

############### Total FFAs ##############
setwd('../')
FFA_tot<-read_excel("FFA.xlsx", sheet=2) 
FFA_tot<-subset.data.frame(FFA_tot, Sample!='BC3')
FFA_tot<-subset.data.frame(FFA_tot, Sample!='BC19')
FFA_tot<-subset.data.frame(FFA_tot, Sample!='BC39')
FFA_tot<-subset.data.frame(FFA_tot, Sample!='BC55')
FFA_tot[! colnames(FFA_tot) %in% c("Sample","Condition","Time")]<-apply(FFA_tot[!colnames(FFA_tot)%in% c("Sample","Condition","Time")], MARGIN = 2, as.numeric)
FFA_tot$Condition<-factor(FFA_tot$Condition, levels = c("INTERV","CONTR"))
FFA_tot$Time<-factor(FFA_tot$Time, levels = c("T0","T1"))

head(FFA_tot)
FFA_cnt<-subset.data.frame(FFA_tot, Condition=='CONTR')
SCFA<-FFA_cnt[,1:4] # from acetic to valeric 
MCFA<-FFA_cnt[,c(1:3,5)] # to dodecanoic
LCFA<-FFA_cnt[,c(1:3,6)] # to Octadecanoic, before the raw total columns
colnames(LCFA)

SCFA_raw<-as.data.frame(SCFA)
MCFA_raw<-as.data.frame(MCFA)
LCFA_raw<- as.data.frame(LCFA)

####################### MANN WITHNEY TESTs ####################
setwd('Results_CERBIOS')
dir.create("FFA_tot")

##### FFA_tot raw abundances
list<-list("SCFA"=SCFA_raw, "MCFA"=MCFA_raw, "LCFA"=LCFA_raw)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("Sample","Condition","Time")]) {
    test<-wilcox.test(temp[[x]]~temp$Time)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  results
  write.csv2(results, file=paste0("FFA_tot/",names(list)[y],"_raw_abund_Mann_Whit_TOT_CNT.csv"), row.names = F)
  
  table <- melt(data = temp)
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Time)) +
    scale_fill_manual(values = c("T0" = alpha("lightgreen", 0.5), "T1" = alpha("green3", 0.5))) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  ggsave(filename = paste0("FFA_tot/",names(list)[y],"_raw_abundance_TOT_CNT.png"), height=3, width = 6, dpi=300)
  
  if(y==2){ # if MCFA
    table_MCFA<-table }
  
  if(y==3){ # if MCFA
    table_LCFA<-table }
}

# specific plot for SCFA
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Total SCFAs"), "", "")
signif_y<-tapply(table_SCFA$value, table_SCFA$variable, max )
table_SCFA$signif_y<-as.numeric(signif_y[table_SCFA$variable]) + 5 # plus 5 to put the sign higher

ggplot(data=table_SCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("lightgreen", 0.5), "T1" = alpha("green3", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="SCFAs concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA_tot/SCFA_zoom_plot_TOT_CNT.png", height=4.5, width = 7, dpi=300)

# specific plot for MCFA
head(table_MCFA, n=2)
unique(table_MCFA$variable)
table_MCFA$signif<-ifelse(table_MCFA$variable %in% c("Total MCFAs"), "", "")   # DA CAMBIARE!
signif_y<-tapply(table_MCFA$value, table_MCFA$variable, max )
table_MCFA$signif_y<-as.numeric(signif_y[table_MCFA$variable]) + 5 # plus 5 to put the sign higher

ggplot(data=table_MCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("lightgreen", 0.5), "T1" = alpha("green3", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="MCFAs concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA_tot/MCFA_zoom_plot_tot_CNT.png", height=4.5, width = 7, dpi=300)

# specific plot for LCFA
head(table_LCFA, n=2)
table_LCFA$signif<-ifelse(table_LCFA$variable %in% c("Total LCFAs"), "", "")   # DA CAMBIARE!
signif_y<-tapply(table_LCFA$value, table_LCFA$variable, max )
table_LCFA$signif_y<-as.numeric(signif_y[table_LCFA$variable]) + 10 # plus 5 to put the sign higher
# View(table_LCFA)

ggplot(data=table_LCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("lightgreen", 0.5), "T1" = alpha("green3", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 700,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="LCFAs concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA_tot/LCFA_zoom_plot_tot_CNT.png", height=4.5, width = 7, dpi=300)

##### FFA raw abundances INTERV
head(FFA)
FFA_INT<-subset.data.frame(FFA, Condition=='INTERV')
FFA_INT
SCFA<-FFA_INT[,1:10] # from acetic to valeric 
MCFA<-FFA_INT[,c(1:3,11:20)] # to dodecanoic
LCFA<-FFA_INT[,c(1:3,21:23)] # to Octadecanoic, before the raw total columns
colnames(LCFA)

SCFA_raw<-as.data.frame(SCFA)
MCFA_raw<-as.data.frame(MCFA)
LCFA_raw<- as.data.frame(LCFA)

list<-list("SCFA"=SCFA_raw, "MCFA"=MCFA_raw, "LCFA"=LCFA_raw)
names(list)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("Sample","Condition","Time")]) {
    test<-wilcox.test(temp[[x]]~temp$Time)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  write.csv2(results, file=paste0("FFA/",names(list)[y],"_raw_abund_Mann_Whit_INT.csv"), row.names = F)
  
  table <- melt(data = temp, id.vars = c("Sample", "Condition","Time"))
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Time)) +
    scale_fill_manual(values = c("T0"="#32CD32","T1"="#FF4500")) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  ggsave(filename = paste0("FFA/",names(list)[y],"_raw_abundance_INT.png"), height=3, width = 6, dpi=300)
  
  if(y==2){ # if MCFA
    table_MCFA<-table }
  
  if(y==3){ # if LCFA
    table_LCFA<-table }
}


# specific plot for SCFA
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Acetic", "Propionic","Butyric","Valeric"), "", "")
signif_y<-tapply(table_SCFA$value, table_SCFA$variable, max )
table_SCFA$signif_y<-as.numeric(signif_y[table_SCFA$variable]) + 10 # plus 5 to put the sign higher

ggplot(data = table_SCFA, mapping = aes(x = variable, y = value, fill = Time)) +
  scale_fill_manual(values = c("T0" = alpha("cyan", 0.5), "T1" = alpha("blue2", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label = signif, y = signif_y, group = variable), 
            color = "darkgray", size = 8) +
  theme(axis.text.x = element_text(angle = -20, size = 9.5, vjust = 1, hjust = 0.05),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 11),
        panel.grid.minor.y = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.35),
        title = element_text(size = 10.5)) +
  scale_y_log10(breaks = c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50, seq(100, 500, 100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width = 0.65, size = 0.5, outlier.size = 1) + 
  labs(title = "SCFA concentration", x = "", fill = "", y = "Percentage (%)")
ggsave(filename = "FFA/SCFA_zoom_plot_INT.png", height=4.5, width = 7, dpi=300)

# specific plot for MCFA
head(table_MCFA, n=2)
unique(table_MCFA$variable)
table_MCFA$signif<-ifelse(table_MCFA$variable %in% c("Octanoic","Decanoic", "Dodecanoic"), "", "")   # DA CAMBIARE!
signif_y<-tapply(table_MCFA$value, table_MCFA$variable, max )
table_MCFA$signif_y<-as.numeric(signif_y[table_MCFA$variable]) + 2 # plus 5 to put the sign higher

ggplot(data=table_MCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values =c("T0" = alpha("cyan", 0.5), "T1" = alpha("blue2", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="MCFA concentration", x="", fill="", y="Concentration (µM/L)")
ggsave(filename = "FFA/MCFA_zoom_plot_INT.png", height=4.5, width = 7, dpi=300)


# specific plot for LCFA
head(table_LCFA, n=2)
table_LCFA$signif<-ifelse(table_LCFA$variable %in% c("Hexadecanoic","Octadecanoic"), "", "")   # DA CAMBIARE!
signif_y<-tapply(table_LCFA$value, table_LCFA$variable, max )
table_LCFA$signif_y<-as.numeric(signif_y[table_LCFA$variable]) + 10 # plus 5 to put the sign higher

# View(table_LCFA)
ggplot(data=table_LCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("cyan", 0.5), "T1" = alpha("blue2", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 700,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="LCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA/LCFA_zoom_plot.png", height=4.5, width = 7, dpi=300)

########### PCoA ########### 
list_PCoA <- list("SCFA" = SCFA_raw, "MCFA" = MCFA_raw, "LCFA" = LCFA_raw)
LCFA_raw$Time

p_values <- NULL
for (x in 1:length(list_PCoA)) {
  rm(data, coord)
  data <- as.data.frame(list_PCoA[[x]])
  row.names(data) <- data$Sample
  
  if (x == 3) {
    data <- data[!is.na(data$Octadecanoic),] # rimuovere NA in LCFA
  }
  
  dist <- vegdist(data[, !colnames(data) %in% c("Sample", "Condition","Time")], method = "bray") # calcolare dissimilarità Bray-Curtis
  coord <- cmdscale(dist, k = 2, eig = TRUE) # PCoA
  
  # Convertire le coordinate in un dataframe e aggiungere Subject e Time
  coord_df <- as.data.frame(coord$points)
  colnames(coord_df) <- c("PC1", "PC2")
  coord_df$Subject <- row.names(coord_df)
  coord_df$Time <- data$Time
  
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Time)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("T0" = alpha("cyan", 0.5), "T1" = alpha("blue2", 0.5))) +
    stat_ellipse(size = 0.2) +
    geom_text(aes(label = Subject), size = 2, color = "black") +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of Intervention subjects\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("PCoA_BRAY_", names(list_PCoA)[x], "INT.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )     
  
  # again but with no names
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Time)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("T0" = alpha("cyan", 0.5), "T1" = alpha("blue2", 0.5))) +
    stat_ellipse(size = 0.2) +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of Intervention subjects\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("PCoA_BRAY_", names(list_PCoA)[x], "INT_noname.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )                        
  
  # checking for significant dispersion
  suppressWarnings(rm(p_value_diver, p_value_disp))
  p_value_diver<-adonis(dist ~ data$Time, data = data, permutations = 9999)
  p_value_disp<-permutest(betadisper(dist, data$Time), permutations = 9999)
  p_value_disp<-as.data.frame(p_value_disp$tab[1,])
  colnames(p_value_disp)<-colnames(p_value_diver$aov.tab[1,])
  
  new_values<-rbind(p_value_diver$aov.tab[1,],p_value_disp)
  row.names(new_values)<-c(paste("Beta_diversity",names(list_PCoA)[x], sep = "_"),
                           paste("Beta_dispersion",names(list_PCoA)[x], sep = "_"))
  p_values<-rbind(p_values,new_values)
}  

p_values$adjusted_BH_p_values<-p.adjust(p_values$`Pr(>F)`, method = "BH")
p_values
write.csv2(p_values, file="p_values_BRAY_INT.csv", quote = F,row.names = T)

############### Total FFAs ##############
head(FFA_tot)
FFA_INT<-subset.data.frame(FFA_tot, Condition=='INTERV')
SCFA<-FFA_INT[,1:4] # from acetic to valeric 
MCFA<-FFA_INT[,c(1:3,5)] # to dodecanoic
LCFA<-FFA_INT[,c(1:3,6)] # to Octadecanoic, before the raw total columns
colnames(LCFA)

SCFA_raw<-as.data.frame(SCFA)
MCFA_raw<-as.data.frame(MCFA)
LCFA_raw<- as.data.frame(LCFA)

####################### MANN WITHNEY TESTs ####################
##### FFA_tot raw abundances
list<-list("SCFA"=SCFA_raw, "MCFA"=MCFA_raw, "LCFA"=LCFA_raw)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("Sample","Condition","Time")]) {
    test<-wilcox.test(temp[[x]]~temp$Time)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  results
  write.csv2(results, file=paste0("FFA_tot/",names(list)[y],"_raw_abund_Mann_Whit_TOT_INT.csv"), row.names = F)
  
  table <- melt(data = temp)
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Time)) +
    scale_fill_manual(values = c("T0" = alpha("cyan", 0.5), "T1" = alpha("blue2", 0.5))) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  ggsave(filename = paste0("FFA_tot/",names(list)[y],"_raw_abundance_TOT_INT.png"), height=3, width = 6, dpi=300)
  
  if(y==2){ # if MCFA
    table_MCFA<-table }
  
  if(y==3){ # if MCFA
    table_LCFA<-table }
}

# specific plot for SCFA
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Total SCFAs"), "", "")
signif_y<-tapply(table_SCFA$value, table_SCFA$variable, max )
table_SCFA$signif_y<-as.numeric(signif_y[table_SCFA$variable]) + 5 # plus 5 to put the sign higher

ggplot(data=table_SCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("cyan", 0.5), "T1" = alpha("blue2", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="SCFAs concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA_tot/SCFA_zoom_plot_TOT_INT.png", height=4.5, width = 7, dpi=300)


# specific plot for MCFA
head(table_MCFA, n=2)
unique(table_MCFA$variable)
table_MCFA$signif<-ifelse(table_MCFA$variable %in% c("Total MCFAs"), "", "")   # DA CAMBIARE!
signif_y<-tapply(table_MCFA$value, table_MCFA$variable, max )
table_MCFA$signif_y<-as.numeric(signif_y[table_MCFA$variable]) + 5 # plus 5 to put the sign higher

ggplot(data=table_MCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("cyan", 0.5), "T1" = alpha("blue2", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="MCFAs concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA_tot/MCFA_zoom_plot_tot_INT.png", height=4.5, width = 7, dpi=300)


# specific plot for LCFA
head(table_LCFA, n=2)
table_LCFA$signif<-ifelse(table_LCFA$variable %in% c("Total LCFAs"), "", "")   # DA CAMBIARE!
signif_y<-tapply(table_LCFA$value, table_LCFA$variable, max )
table_LCFA$signif_y<-as.numeric(signif_y[table_LCFA$variable]) + 10 # plus 5 to put the sign higher
# View(table_LCFA)

ggplot(data=table_LCFA, mapping=aes(x=variable, y=value, fill=Time)) +
  scale_fill_manual(values = c("T0" = alpha("cyan", 0.5), "T1" = alpha("blue2", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 700,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="LCFAs concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA_tot/LCFA_zoom_plot_tot_INT.png", height=4.5, width = 7, dpi=300)

##### FFA raw abundances T1
head(FFA)
FFA_T1<-subset.data.frame(FFA, Time=='T1')
FFA_T1
SCFA<-FFA_T1[,1:10] # from acetic to valeric 
MCFA<-FFA_T1[,c(1:3,11:20)] # to dodecanoic
LCFA<-FFA_T1[,c(1:3,21:23)] # to Octadecanoic, before the raw total columns
colnames(LCFA)

SCFA_raw<-as.data.frame(SCFA)
MCFA_raw<-as.data.frame(MCFA)
LCFA_raw<- as.data.frame(LCFA)

list<-list("SCFA"=SCFA_raw, "MCFA"=MCFA_raw, "LCFA"=LCFA_raw)
names(list)

for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("Sample","Condition","Time")]) {
    test<-wilcox.test(temp[[x]]~temp$Condition)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  write.csv2(results, file=paste0("FFA/",names(list)[y],"_raw_abund_Mann_Whit_T1.csv"), row.names = F)
  
  table <- melt(data = temp, id.vars = c("Sample", "Condition","Time"))
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Condition)) +
    scale_fill_manual(values = c("CONTR"="#32CD32","INTERV"="#FF4500")) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  ggsave(filename = paste0("FFA/",names(list)[y],"_raw_abundance_T1.png"), height=3, width = 6, dpi=300)
  
  if(y==2){ # if MCFA
    table_MCFA<-table }
  
  if(y==3){ # if LCFA
    table_LCFA<-table }
}


# specific plot for SCFA
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Acetic", "Propionic","Butyric","Valeric"), "", "")
signif_y<-tapply(table_SCFA$value, table_SCFA$variable, max )
table_SCFA$signif_y<-as.numeric(signif_y[table_SCFA$variable]) + 10 # plus 5 to put the sign higher

ggplot(data = table_SCFA, mapping = aes(x = variable, y = value, fill = Condition)) +
  scale_fill_manual(values = c("CONTR" = alpha("green3", 0.5), "INTERV" = alpha("blue2", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label = signif, y = signif_y, group = variable), 
            color = "darkgray", size = 8) +
  theme(axis.text.x = element_text(angle = -20, size = 9.5, vjust = 1, hjust = 0.05),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 11),
        panel.grid.minor.y = element_line(linewidth = 0.25),
        panel.grid.major.y = element_line(linewidth = 0.35),
        title = element_text(size = 10.5)) +
  scale_y_log10(breaks = c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50, seq(100, 500, 100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width = 0.65, size = 0.5, outlier.size = 1) + 
  labs(title = "SCFA concentration", x = "", fill = "", y = "Percentage (%)")
ggsave(filename = "FFA/SCFA_zoom_plot_T1.png", height=4.5, width = 7, dpi=300)


# specific plot for MCFA
head(table_MCFA, n=2)
unique(table_MCFA$variable)
table_MCFA$signif<-ifelse(table_MCFA$variable %in% c("Octanoic","Decanoic", "Dodecanoic"), "", "")   # DA CAMBIARE!
signif_y<-tapply(table_MCFA$value, table_MCFA$variable, max )
table_MCFA$signif_y<-as.numeric(signif_y[table_MCFA$variable]) + 2 # plus 5 to put the sign higher

ggplot(data=table_MCFA, mapping=aes(x=variable, y=value, fill=Condition)) +
  scale_fill_manual(values = c("CONTR" = alpha("green3", 0.5), "INTERV" = alpha("blue2", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="MCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA/MCFA_zoom_plot_T1.png", height=4.5, width = 7, dpi=300)


# specific plot for LCFA
head(table_LCFA, n=2)
table_LCFA$signif<-ifelse(table_LCFA$variable %in% c("Hexadecanoic","Octadecanoic"), "", "")   # DA CAMBIARE!
signif_y<-tapply(table_LCFA$value, table_LCFA$variable, max )
table_LCFA$signif_y<-as.numeric(signif_y[table_LCFA$variable]) + 10 # plus 5 to put the sign higher

# View(table_LCFA)
ggplot(data=table_LCFA, mapping=aes(x=variable, y=value, fill=Condition)) +
  scale_fill_manual(values = c("CONTR" = alpha("green3", 0.5), "INTERV" = alpha("blue2", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 700,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="LCFA concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA/LCFA_zoom_plot_T1.png", height=4.5, width = 7, dpi=300)

########### PCoA ########### 
list_PCoA <- list("SCFA" = SCFA_raw, "MCFA" = MCFA_raw, "LCFA" = LCFA_raw)
LCFA_raw$Condition

p_values <- NULL
for (x in 1:length(list_PCoA)) {
  rm(data, coord)
  data <- as.data.frame(list_PCoA[[x]])
  row.names(data) <- data$Sample
  
  if (x == 3) {
    data <- data[!is.na(data$Octadecanoic),] # rimuovere NA in LCFA
  }
  
  dist <- vegdist(data[, !colnames(data) %in% c("Sample", "Condition","Time")], method = "bray") # calcolare dissimilarità Bray-Curtis
  coord <- cmdscale(dist, k = 2, eig = TRUE) # PCoA
  
  # Convertire le coordinate in un dataframe e aggiungere Subject e Time
  coord_df <- as.data.frame(coord$points)
  colnames(coord_df) <- c("PC1", "PC2")
  coord_df$Subject <- row.names(coord_df)
  coord_df$Condition <- data$Condition
  
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Condition)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values =c("CONTR" = alpha("green3", 0.5), "INTERV" = alpha("blue2", 0.5))) +
    stat_ellipse(size = 0.2) +
    geom_text(aes(label = Subject), size = 2, color = "black") +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of T1 subjects\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("PCoA_BRAY_", names(list_PCoA)[x], "T1.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )     
  
  # again but with no names
  ggplot(data = coord_df, aes(x = PC1, y = PC2, color = Condition)) + 
    geom_point(size = 3, alpha = 1) +
    theme_classic() + 
    scale_color_manual(values = c("CONTR" = alpha("green3", 0.5), "INTERV" = alpha("blue2", 0.5))) +
    stat_ellipse(size = 0.2) +
    labs(
      x = paste("PC1:", round(coord$eig[1] / sum(coord$eig) * 100, digits = 2), "%"),
      y = paste("PC2:", round(coord$eig[2] / sum(coord$eig) * 100, digits = 2), "%"),
      color = "", 
      title = paste("PCoA computed on", names(list_PCoA)[x], "quantities\n of T1 subjects\n with Bray-Curtis dissimilarity")
    )
  
  ggsave(
    file = paste("PCoA_BRAY_", names(list_PCoA)[x], "T1_noname.png", sep = "_"), 
    height = 5, 
    width = 6.5, 
    dpi = 300
  )                        
  
  # checking for significant dispersion
  suppressWarnings(rm(p_value_diver, p_value_disp))
  p_value_diver<-adonis(dist ~ data$Condition, data = data, permutations = 9999)
  p_value_disp<-permutest(betadisper(dist, data$Condition), permutations = 9999)
  p_value_disp<-as.data.frame(p_value_disp$tab[1,])
  colnames(p_value_disp)<-colnames(p_value_diver$aov.tab[1,])
  
  new_values<-rbind(p_value_diver$aov.tab[1,],p_value_disp)
  row.names(new_values)<-c(paste("Beta_diversity",names(list_PCoA)[x], sep = "_"),
                           paste("Beta_dispersion",names(list_PCoA)[x], sep = "_"))
  p_values<-rbind(p_values,new_values)
}  

p_values$adjusted_BH_p_values<-p.adjust(p_values$`Pr(>F)`, method = "BH")
p_values
write.csv2(p_values, file="p_values_BRAY_T1.csv", quote = F,row.names = T)

############### Total FFAs ##############
head(FFA_tot)
FFA_T1<-subset.data.frame(FFA_tot, Time=='T1')
SCFA<-FFA_T1[,1:4] # from acetic to valeric 
MCFA<-FFA_T1[,c(1:3,5)] # to dodecanoic
LCFA<-FFA_T1[,c(1:3,6)] # to Octadecanoic, before the raw total columns
colnames(LCFA)

SCFA_raw<-as.data.frame(SCFA)
MCFA_raw<-as.data.frame(MCFA)
LCFA_raw<- as.data.frame(LCFA)

####################### MANN WITHNEY TESTs ####################
##### FFA_tot raw abundances
list<-list("SCFA"=SCFA_raw, "MCFA"=MCFA_raw, "LCFA"=LCFA_raw)
SCFA_raw$Time
for(y in 1:length(list)){
  temp<-list[[y]]
  results<-NULL
  
  for(x in colnames(temp)[!colnames(temp)%in% c("Sample","Condition","Time")]) {
    test<-wilcox.test(temp[[x]]~temp$Condition)
    wilc<-cbind(paste(x),test$statistic,test$p.value)
    results<-rbind.data.frame(results,wilc)
  }
  colnames(results)<-c(paste(names(list)[y]),"W-statistic","p.value")
  results$'p.adj(BH)'<-p.adjust(results$p.value, method = "BH")
  row.names(results)<-NULL
  results
  write.csv2(results, file=paste0("FFA_tot/",names(list)[y],"_raw_abund_Mann_Whit_TOT_T1.csv"), row.names = F)
  
  table <- melt(data = temp)
  if(y==1){ # if SCFA
    table_SCFA<-table } # to build a specific plot
  ggplot(data=table, mapping=aes(x=variable, y=value, fill=Condition)) +
    scale_fill_manual(values =c("CONTR" = alpha("green3", 0.5), "INTERV" = alpha("blue2", 0.5))) +
    theme( axis.text.x = element_text(angle=-25, size = 8, vjust=1, hjust = 0.05))+
    geom_boxplot() + labs(title=paste("Raw",names(list)[y],"quantities",sep=" "), x="", y="Quantity")
  ggsave(filename = paste0("FFA_tot/",names(list)[y],"_raw_abundance_TOT_T1.png"), height=3, width = 6, dpi=300)
  
  if(y==2){ # if MCFA
    table_MCFA<-table }
  
  if(y==3){ # if MCFA
    table_LCFA<-table }
}

# specific plot for SCFA
head(table_SCFA, n=2)
table_SCFA$signif<-ifelse(table_SCFA$variable %in% c("Total SCFAs"), "", "")
signif_y<-tapply(table_SCFA$value, table_SCFA$variable, max )
table_SCFA$signif_y<-as.numeric(signif_y[table_SCFA$variable]) + 5 # plus 5 to put the sign higher

ggplot(data=table_SCFA, mapping=aes(x=variable, y=value, fill=Condition)) +
  scale_fill_manual(values = c("CONTR" = alpha("green3", 0.5), "INTERV" = alpha("blue2", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="SCFAs concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA_tot/SCFA_zoom_plot_TOT_T1.png", height=4.5, width = 7, dpi=300)


# specific plot for MCFA
head(table_MCFA, n=2)
unique(table_MCFA$variable)
table_MCFA$signif<-ifelse(table_MCFA$variable %in% c("Total MCFAs"), "", "")   # DA CAMBIARE!
signif_y<-tapply(table_MCFA$value, table_MCFA$variable, max )
table_MCFA$signif_y<-as.numeric(signif_y[table_MCFA$variable]) + 5 # plus 5 to put the sign higher

ggplot(data=table_MCFA, mapping=aes(x=variable, y=value, fill=Condition)) +
  scale_fill_manual(values = c("CONTR" = alpha("green3", 0.5), "INTERV" = alpha("blue2", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 500,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="MCFAs concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA_tot/MCFA_zoom_plot_tot_T1.png", height=4.5, width = 7, dpi=300)


# specific plot for LCFA
head(table_LCFA, n=2)
table_LCFA$signif<-ifelse(table_LCFA$variable %in% c("Total LCFAs"), "", "")   # DA CAMBIARE!
signif_y<-tapply(table_LCFA$value, table_LCFA$variable, max )
table_LCFA$signif_y<-as.numeric(signif_y[table_LCFA$variable]) + 10 # plus 5 to put the sign higher
# View(table_LCFA)

ggplot(data=table_LCFA, mapping=aes(x=variable, y=value, fill=Condition)) +
  scale_fill_manual(values = c("CONTR" = alpha("green3", 0.5), "INTERV" = alpha("blue2", 0.5))) +
  theme_bw(base_size = 15) +
  geom_text(aes(label=signif, y = signif_y,
                group= variable), 
            color="darkgray", size=8) +
  theme( axis.text.x = element_text(angle=-20, size = 9.5, vjust=1, hjust = 0.05),
         axis.text.y = element_text(size = 7),
         axis.title.y = element_text(size = 11),
         panel.grid.minor.y = element_line(linewidth = 0.25),
         panel.grid.major.y = element_line(linewidth = 0.35),
         title = element_text(size = 10.5)) +
  scale_y_log10(breaks= c(0.2, 0.5, 1, 2.5, 5, 10, 25, 50,
                          seq(100, 700,100))) +
  scale_x_discrete(expand = c(0, 0.4)) +
  geom_boxplot(width=0.65, size=0.5, outlier.size = 1) + 
  labs(title="LCFAs concentration", x="", fill="", y="Percentage (%)")
ggsave(filename = "FFA_tot/LCFA_zoom_plot_tot_T1.png", height=4.5, width = 7, dpi=300)

################## Tabella RA per ogni pz ################## 
#NB per farla di tutti devi lanciare lo script prima di rimuovere gli altri ma facendo il filter
physeq <- data
tax_table_physeq <- tax_table(physeq)
tax_table_physeq[ tax_table_physeq == "d__Bacteria" ] <- "Bacteria"
tax_table(physeq) <- tax_table_physeq
samples <- as.character(sample_data(physeq)$ID)  
tax_table(physeq)
sample_data(physeq)
# Directory to save output files
#setwd('Results_CERBIOS')
output_dir <- "RA_patients"
dir.create(output_dir, showWarnings = FALSE)
sample_data(physeq)

for (i in samples) {
  physeq_patient <- subset_samples(physeq, ID == i)
  
  physeq_patient_rel <- tax_glom(physeq_patient, taxrank = "Genus")
  physeq_patient_rel <- transform_sample_counts(physeq_patient_rel, function(x) x / sum(x))
 
  melted_data <- psmelt(physeq_patient_rel)
  otu_table_patient <- melted_data [,c("ID", "Condition","Kingdom","Phylum","Class","Order","Family","Genus","Abundance")]
  otu_table_patient <- otu_table_patient[order(-otu_table_patient$Abundance), ]
  
  write.xlsx(otu_table_patient,
             file = file.path(output_dir, paste0("Patient_", i, "_OTU_Table.xlsx")),
             row.names = FALSE)
}

cat("Files saved in directory:", output_dir)

#################### CORRELAZIONI #################### 
# estrazione RA Lachnospirales per correlazioni
sample_data(data)
sub<-tax_glom(data.prop, taxrank = 'Order')
sub<-subset_taxa(sub, Order=='Lachnospirales')
tax_table(sub)
o<-as.data.frame(otu_table(sub))
o<-t(o)

write.csv(o, 'RA_Lachnospirales.csv')

sample_data(data)

#FFA (fecali)
siero<-read_excel('FFA.xlsx', sheet = 4)

head(siero)
siero<-as.data.frame(siero)
row.names(siero)<-siero$Sample
siero$Sample

siero$Condition
siero_int<-subset.data.frame(siero, Condition=='INTERV')
siero_cnt<-subset.data.frame(siero, Condition=='CONTR')

#CLINICI body composition
#setwd('../../CERBIOS')
clinici<-read_excel('CERBIOS_clinical.xlsx', sheet = 2)
head(clinici)
clinici<-as.data.frame(clinici)
row.names(clinici)<-clinici$Sample
clinici$Sample
setdiff(siero$Sample,clinici$Sample)

clinici$Condition
clinici_int<-subset.data.frame(clinici, Condition=='INTERV')
clinici_cnt<-subset.data.frame(clinici, Condition=='CONTR')

length(siero_cnt$Sample) #15 cnt
length(clinici_cnt$Sample)
setdiff(siero_cnt$Sample,clinici_cnt$Sample)
length(siero_int$Sample) #17 int
length(clinici_int$Sample)

# FFA vs body composition 
nuovo<-cbind(siero_cnt,clinici_cnt)
nuovo<-nuovo[,c(3:23,27:30)]
results <- data.frame(matrix(ncol=4, nrow=0))
colnames(results)<-c("Variabile_1","Variabile_2","Correlazione","P_value")

for (i in 1:(ncol(nuovo) - 1)) {
  for (j in i:ncol(nuovo)) {
    cor.test_result <- cor.test(nuovo[, i], nuovo[, j], method = "spearman", complete.obs = TRUE)
    cor_value <- as.numeric(cor.test_result$estimate)
    if (is.numeric(cor_value) && cor_value >= -1 && cor_value <= 1) {
      results[nrow(results) + 1, ] <- c(
        names(nuovo)[i],
        names(nuovo)[j],
        cor_value,
        cor.test_result$p.value
      )
    } else {
      cat("Errore con i valori: ", names(nuovo)[i], names(nuovo)[j], "\n")
    }
  }
}


vec<-colnames(siero_cnt)
prova<- subset.data.frame(results, Variabile_1 %in% vec & !Variabile_2 %in%vec)

prova$padj<-p.adjust(prova$P_value, method = "BH")
prova$Significative<- ifelse(prova$padj<0.05, "*", "")
prova #no Results_Ceren/Correlazioni significative
#dir.create('Results_CERBIOS/Correlazioni')
write.xlsx (prova, 'Results_CERBIOS/Correlazioni/Corr_FFA_body_comp.xlsx')

# FFA vs body composition INTERV
nuovo<-cbind(siero_int,clinici_int)
length(siero_int$Sample)
length(clinici_int$Sample)
nuovo<-nuovo[,c(3:23,27:30)]
setdiff(siero_int$Sample,clinici_int$Sample)
results <- data.frame(matrix(ncol=4, nrow=0))
colnames(results)<-c("Variabile_1","Variabile_2","Correlazione","P_value")

for (i in 1:(ncol(nuovo) - 1)) {
  for (j in i:ncol(nuovo)) {
    cor.test_result <- cor.test(nuovo[, i], nuovo[, j], method = "spearman", complete.obs = TRUE)
    cor_value <- as.numeric(cor.test_result$estimate)
    if (is.numeric(cor_value) && cor_value >= -1 && cor_value <= 1) {
      results[nrow(results) + 1, ] <- c(
        names(nuovo)[i],
        names(nuovo)[j],
        cor_value,
        cor.test_result$p.value
      )
    } else {
      cat("Errore con i valori: ", names(nuovo)[i], names(nuovo)[j], "\n")
    }
  }
}

vec<-colnames(siero_cnt)
prova<- subset.data.frame(results, Variabile_1 %in% vec & !Variabile_2 %in%vec)

prova$padj<-p.adjust(prova$P_value, method = "BH")
prova$Significative<- ifelse(prova$padj<0.05, "*", "")
prova #no Results_Ceren/Correlazioni significative
write.xlsx (prova, 'Results_CERBIOS/Correlazioni/Corr_FFA_body_comp_int.xlsx')

#HEATMAP FFA AND BODY COMPOSITION
library(ggplot2)
library(reshape2)
library(forcats)  

dati_long <- read_xlsx('Results_CERBIOS/Correlazioni/Corr_FFA_body_comp_tot.xlsx',sheet=1)
dati_long$Correlation<-as.numeric(dati_long$Correlation)
dati_long$padj<-as.numeric(dati_long$padj)

#plot stool FFAs
ggplot(dati_long, aes(x = fct_inorder(Variabile_2), y = fct_inorder(Variabile_1), fill = Correlation)) +
  geom_tile() +  # Rappresenta i dati come quadrati
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Imposta la scala dei colori
  geom_text(aes(label = ifelse(padj < 0.05, ifelse(padj < 0.01, ifelse(padj < 0.001, "***", "**"), "*"), "")), size = 3)+  # Aggiungi l'asterisco per i p-value significativi
  labs(x = "", y = "Faecal fatty acids and Bacteria", title = "Heatmap Correlation - Body composition") +  # Etichette degli assi e del titolo
  theme_minimal() +  # Stile del grafico
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  # Rotazione etichette sull'asse x
  facet_wrap(~ Condition)  # Dividi le variabili in base alla colonna di trattamento

#ggsave(filename = paste0("Heatmap_HC_FSHD_corr_feci.png"), width = 9, height = 5, dpi = 300)

#FFA
siero<-read_excel('FFA.xlsx', sheet = 4)

head(siero)
siero<-as.data.frame(siero)
row.names(siero)<-siero$Sample
siero$Sample

siero$Condition
siero_int<-subset.data.frame(siero, Condition=='INTERV')
siero_cnt<-subset.data.frame(siero, Condition=='CONTR')

#CLINICI biochemical parameters
#setwd('../../CERBIOS')
clinici<-read_excel('CERBIOS_clinical.xlsx', sheet = 3)
head(clinici)
clinici<-as.data.frame(clinici)
row.names(clinici)<-clinici$Sample
clinici$Sample

clinici$Condition
clinici_int<-subset.data.frame(clinici, Condition=='INTERV')
clinici_cnt<-subset.data.frame(clinici, Condition=='CONTR')

length(siero_cnt$Sample)
length(clinici_cnt$Sample)
setdiff(siero_cnt$Sample,clinici_cnt$Sample)
length(siero_int$Sample)
length(clinici_int$Sample)

# FFA vs biochemical parameters
nuovo<-cbind(siero_cnt,clinici_cnt)
nuovo<-nuovo[,c(3:23,27:37)]
results <- data.frame(matrix(ncol=4, nrow=0))
colnames(results)<-c("Variabile_1","Variabile_2","Correlazione","P_value")

cor.test(nuovo$Plt_109L, nuovo$Lachnospirales, method = "spearman", complete.obs = TRUE)

for (i in 1:(ncol(nuovo) - 1)) {
  for (j in i:ncol(nuovo)) {
    cor.test_result <- cor.test(nuovo[, i], nuovo[, j], method = "spearman", complete.obs = TRUE)
    cor_value <- as.numeric(cor.test_result$estimate)
    if (is.numeric(cor_value) && cor_value >= -1 && cor_value <= 1) {
      results[nrow(results) + 1, ] <- c(
        names(nuovo)[i],
        names(nuovo)[j],
        cor_value,
        cor.test_result$p.value
      )
    } else {
      cat("Errore con i valori: ", names(nuovo)[i], names(nuovo)[j], "\n")
    }
  }
}

vec<-colnames(siero_cnt)
prova<- subset.data.frame(results, Variabile_1 %in% vec & !Variabile_2 %in%vec)

prova$padj<-p.adjust(prova$P_value, method = "BH")
prova$Significative<- ifelse(prova$padj<0.05, "*", "")
prova #no Results_Ceren/Correlazioni significative
#dir.create('Results_CERBIOS/Correlazioni')
write.xlsx (prova, 'Results_CERBIOS/Correlazioni/Corr_FFA_bioch_para.xlsx')

# FFA vs biochemical parameters INTERV
nuovo<-cbind(siero_int,clinici_int)
length(siero_int$Sample)
length(clinici_int$Sample)
nuovo<-nuovo[,c(3:23,27:37)]
setdiff(siero_int$Sample,clinici_int$Sample)
results <- data.frame(matrix(ncol=4, nrow=0))
colnames(results)<-c("Variabile_1","Variabile_2","Correlazione","P_value")

for (i in 1:(ncol(nuovo) - 1)) {
  for (j in i:ncol(nuovo)) {
    cor.test_result <- cor.test(nuovo[, i], nuovo[, j], method = "spearman", complete.obs = TRUE)
    cor_value <- as.numeric(cor.test_result$estimate)
    if (is.numeric(cor_value) && cor_value >= -1 && cor_value <= 1) {
      results[nrow(results) + 1, ] <- c(
        names(nuovo)[i],
        names(nuovo)[j],
        cor_value,
        cor.test_result$p.value
      )
    } else {
      cat("Errore con i valori: ", names(nuovo)[i], names(nuovo)[j], "\n")
    }
  }
}
vec<-colnames(siero_cnt)
prova<- subset.data.frame(results, Variabile_1 %in% vec & !Variabile_2 %in%vec)

prova$padj<-p.adjust(prova$P_value, method = "BH")
prova$Significative<- ifelse(prova$padj<0.05, "*", "")
prova #no Results_Ceren/Correlazioni significative
write.xlsx (prova, 'Results_CERBIOS/Correlazioni/Corr_FFA_body_bioch_par_INT.xlsx')

#HEATMAP FFA AND BODY COMPOSITION
library(ggplot2)
library(reshape2)
library(forcats)  

dati_long <- read_xlsx('Results_CERBIOS/Correlazioni/Corr_FFA_body_bioch_para_tot.xlsx',sheet=1)
dati_long$Correlation<-as.numeric(dati_long$Correlation)
dati_long$padj<-as.numeric(dati_long$padj)

pro$Correlation <- sign(pro$Correlation) * log10(abs(pro$Correlation) + 1)
#plot stool FFAs
ggplot(dati_long, aes(x = fct_inorder(Variabile_2), y = fct_inorder(Variabile_1), fill = Correlation)) +
  geom_tile() +  # Rappresenta i dati come quadrati
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Imposta la scala dei colori
  geom_text(aes(label = ifelse(padj < 0.05, ifelse(padj < 0.01, ifelse(padj < 0.001, "***", "**"), "*"), "")), size = 3)+  # Aggiungi l'asterisco per i p-value significativi
  labs(x = "", y = "Faecal fatty acids and Bacteria", title = "Heatmap Correlation - Biochemical parameters") +  # Etichette degli assi e del titolo
  theme_minimal() +  # Stile del grafico
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  # Rotazione etichette sull'asse x
  facet_wrap(~ Condition)  # Dividi le variabili in base alla colonna di trattamento

#ggsave(filename = paste0("Heatmap_HC_FSHD_corr_feci.png"), width = 9, height = 5, dpi = 300)

#FFA
siero<-read_excel('FFA.xlsx', sheet = 4)

head(siero)
siero<-as.data.frame(siero)
row.names(siero)<-siero$Sample
siero$Sample

siero$Condition
siero_int<-subset.data.frame(siero, Condition=='INTERV')
siero_cnt<-subset.data.frame(siero, Condition=='CONTR')

#CLINICI body composition
#setwd('../../CERBIOS')
clinici<-read_excel('CERBIOS_clinical.xlsx', sheet = 4)
head(clinici)
clinici<-as.data.frame(clinici)
row.names(clinici)<-clinici$Sample
clinici$Sample

clinici$Condition
clinici_int<-subset.data.frame(clinici, Condition=='INTERV')
clinici_cnt<-subset.data.frame(clinici, Condition=='CONTR')

length(siero_cnt$Sample)
length(clinici_cnt$Sample)
setdiff(siero_cnt$Sample,clinici_cnt$Sample)
length(siero_int$Sample)
length(clinici_int$Sample)

# FFA vs body composition 
nuovo<-cbind(siero_cnt,clinici_cnt)
nuovo<-nuovo[,c(3:23,27:31)]
results <- data.frame(matrix(ncol=4, nrow=0))
colnames(results)<-c("Variabile_1","Variabile_2","Correlazione","P_value")
summary(nuovo)
cor.test(nuovo$Acetic, nuovo$IL18, method = "spearman", complete.obs = TRUE)

for (i in 1:(ncol(nuovo) - 1)) {
  for (j in i:ncol(nuovo)) {
    cor.test_result <- cor.test(nuovo[, i], nuovo[, j], method = "spearman", complete.obs = TRUE)
    cor_value <- as.numeric(cor.test_result$estimate)
    if (is.numeric(cor_value) && cor_value >= -1 && cor_value <= 1) {
      results[nrow(results) + 1, ] <- c(
        names(nuovo)[i],
        names(nuovo)[j],
        cor_value,
        cor.test_result$p.value
      )
    } else {
      cat("Errore con i valori: ", names(nuovo)[i], names(nuovo)[j], "\n")
    }
  }
}

vec<-colnames(siero_cnt)
prova<- subset.data.frame(results, Variabile_1 %in% vec & !Variabile_2 %in%vec)

prova$padj<-p.adjust(prova$P_value, method = "BH")
prova$Significative<- ifelse(prova$padj<0.05, "*", "")
prova #no Results_Ceren/Correlazioni significative
#dir.create('Results_CERBIOS/Correlazioni')
write.xlsx (prova, 'Results_CERBIOS/Correlazioni/Corr_FFA_cito.xlsx')

# FFA vs body composition INTERV
nuovo<-cbind(siero_int,clinici_int)
length(siero_int$Sample)
length(clinici_int$Sample)
nuovo<-nuovo[,c(3:23,27:31)]
setdiff(siero_int$Sample,clinici_int$Sample)
results <- data.frame(matrix(ncol=4, nrow=0))
colnames(results)<-c("Variabile_1","Variabile_2","Correlazione","P_value")

for (i in 1:(ncol(nuovo) - 1)) {
  for (j in i:ncol(nuovo)) {
    cor.test_result <- cor.test(nuovo[, i], nuovo[, j], method = "spearman", complete.obs = TRUE)
    cor_value <- as.numeric(cor.test_result$estimate)
    if (is.numeric(cor_value) && cor_value >= -1 && cor_value <= 1) {
      results[nrow(results) + 1, ] <- c(
        names(nuovo)[i],
        names(nuovo)[j],
        cor_value,
        cor.test_result$p.value
      )
    } else {
      cat("Errore con i valori: ", names(nuovo)[i], names(nuovo)[j], "\n")
    }
  }
}

vec<-colnames(siero_cnt)
prova<- subset.data.frame(results, Variabile_1 %in% vec & !Variabile_2 %in%vec)

prova$padj<-p.adjust(prova$P_value, method = "BH")
prova$Significative<- ifelse(prova$padj<0.05, "*", "")
prova #no Results_Ceren/Correlazioni significative
write.xlsx (prova, 'Results_CERBIOS/Correlazioni/Corr_FFA_body_cito_INT.xlsx')

#HEATMAP FFA AND BODY COMPOSITION
library(ggplot2)
library(reshape2)
library(forcats)  

dati_long <- read_xlsx('Results_CERBIOS/Correlazioni/Corr_FFA_body_cito_tot.xlsx',sheet=1)
dati_long$Correlation<-as.numeric(dati_long$Correlation)
dati_long$padj<-as.numeric(dati_long$padj)

ggplot(dati_long, aes(x = fct_inorder(Variabile_2), y = fct_inorder(Variabile_1), fill = Correlation)) +
  geom_tile() +  # Rappresenta i dati come quadrati
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Imposta la scala dei colori
  geom_text(aes(label = ifelse(padj < 0.05, ifelse(padj < 0.01, ifelse(padj < 0.001, "***", "**"), "*"), "")), size = 3)+  # Aggiungi l'asterisco per i p-value significativi
  labs(x = "", y = "Faecal fatty acids and Bacteria", title = "Heatmap Correlation - Cytokines") +  # Etichette degli assi e del titolo
  theme_minimal() +  # Stile del grafico
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  # Rotazione etichette sull'asse x
  facet_wrap(~ Condition)  # Dividi le variabili in base alla colonna di trattamento

#ggsave(filename = paste0("Heatmap_HC_FSHD_corr_feci.png"), width = 9, height = 5, dpi = 300)