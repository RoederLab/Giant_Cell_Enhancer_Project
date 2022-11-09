#manually load in dataset "GCE_Image_processing_data_RATIOs_GCthresh50_T6.csv"

#Rename dataset to something smaller
Cell_counts = `GCE_Image_processing_data_RATIOs_GCthresh50_T6`

#Load in required libraries
library(dplyr)
library(ggplot2)
library(agricolae)
library(EnvStats)

#Subset the first set of data of interest
Bashing = subset(Cell_counts, Genotype %in% c("pAR111", "pAR254", "pAR257","pAR260","pAR261","pAR262"))

#Perform ANOVA and Tukey HSD Test on subsetted data
Signif_bash = aov(Small.Giant_ratio ~ Genotype, data = Bashing)
summary(Signif_bash)
hsd = HSD.test(Signif_bash, "Genotype")

#group data by genotype and extract the significance letters from HSD test
Signif_summary = Bashing %>% group_by(Genotype) %>% summarize(Mean.ratio=mean(Small.Giant_ratio))
signif_letters = data.frame(row.names(hsd$groups),hsd$groups$groups,hsd$groups$Small.Giant_ratio)

#reorder the the significance letters by alphabetical order of the genotypes so the order is the same as the plot
signif_letters = signif_letters %>% arrange(signif_letters$hsd.groups.Small.Giant_ratio)

#Make the Boxplot with Significance letters
######   FIGURE 2  #######
ggplot(Bashing, aes(x=reorder(Genotype, Small.Giant_ratio), y=Small.Giant_ratio, color = Genotype)) + 
  geom_boxplot(lwd = 1.5) +
  theme_classic() +
  #Color in bar graphs manually 
  scale_color_manual(values = c("blue","purple", "#01af50","purple","purple","blue")) +
  geom_text(data=Signif_summary, aes(x=signif_letters$row.names.hsd.groups.,y=75,label =signif_letters$hsd.groups.groups), color = "black", vjust=0, size = 8) +
  #labs(x ="Genotype", y = "Small Cell : Giant Cell Ratio") +
  theme(text = element_text(size = 20),axis.title.x=element_blank(),axis.title.y=element_blank(), legend.position = "none") +
  #Add Sample Sizes
  stat_n_text(size = 6)
#Exported as 900 x 900 Pixel Tif


#Subset the second set of data of interest
Dofs = subset(Cell_counts, Genotype %in% c("pAR111", "pAR254", "pAR257","pAR307","pAR308","pLH166"))

#Perform one way ANOVA and Tukey HSD Test on subsetted data
Signif_dofs = aov(Small.Giant_ratio ~ Genotype, data = Dofs)
summary(Signif_dofs)
hsd2 = HSD.test(Signif_dofs, "Genotype")

#group data by genotype and extract the significance letters from HSD test
Signif_summary2 = Dofs %>% group_by(Genotype) %>% summarize(Mean.ratio=mean(Small.Giant_ratio))
signif_letters2 = data.frame(row.names(hsd2$groups),hsd2$groups$groups,hsd2$groups$`Small.Giant_ratio`)

#reorder the the significance letters by alphabetical order of the genotypes so the order is the same as the plot
signif_letters2 = signif_letters2 %>% arrange(signif_letters2$hsd2.groups.Small.Giant_ratio)

######   FIGURE 4  #######
#Make the Boxplot with Significance letters
ggplot(Dofs, aes(x=reorder(Genotype, Small.Giant_ratio), y=Small.Giant_ratio, color = Genotype)) + 
  geom_boxplot(lwd = 1.5) +
  theme_classic() +
  scale_color_manual(values = c("blue","purple", "#01af50","#01af50","purple","purple")) +
  geom_text(data=Signif_summary2, aes(x=signif_letters2$row.names.hsd2.groups.,y=45,label =signif_letters2$hsd2.groups.groups),vjust=0, color = "black", size = 8) +
  labs(x ="Genotype", y = "Small Cell : Giant Cell Ratio") +
  theme(text = element_text(size = 20), axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
  stat_n_text(size = 6)
#Exported as 1350 x 900 Pixel Tif



##### FIGURE S5 TRIHELIX MUTATIONS ######

#Subset the first set of data of interest
Trihelix = subset(Cell_counts, Genotype %in% c("pAR111", "pAR307", "pAR388","pAR389")) #,"pLH166","pLH167"))

#Perform ANOVA and Tukey HSD Test on subsetted data
Signif_helix = aov(Small.Giant_ratio ~ Genotype, data = Trihelix)
summary(Signif_helix)
hsd = HSD.test(Signif_helix, "Genotype")

#group data by genotype and extract the significance letters from HSD test
Signif_summary = Trihelix %>% group_by(Genotype) %>% summarize(Mean.ratio=mean(Small.Giant_ratio))
signif_letters = data.frame(row.names(hsd$groups),hsd$groups$groups,hsd$groups$Small.Giant_ratio)

#reorder the the significance letters by alphabetical order of the genotypes so the order is the same as the plot
signif_letters = signif_letters %>% arrange(signif_letters$hsd.groups.Small.Giant_ratio)

#Compare groups with T-Tests
t.test(Small.Giant_ratio ~ Genotype, data = subset(Trihelix, Genotype %in% c("pAR111","pAR389"))) #.33 NS
t.test(Small.Giant_ratio ~ Genotype, data = subset(Trihelix, Genotype %in% c("pAR307","pAR388"))) #.58 NS

Trihelix$Genotype <- factor(Trihelix$Genotype,     # Reorder factor levels
                         c("pAR111", 'pAR389', "pAR307","pAR388"))#,"pLH166","pLH167"))

#Make the Boxplot with Significance letters
ggplot(Trihelix, aes(x=Genotype, y=Small.Giant_ratio, color = Genotype)) + 
  geom_boxplot() +
  theme_classic() +
  #Color in bar graphs manually 
  scale_color_manual(values = c("blue","blue","purple","purple")) +
  geom_text(data=Signif_summary, aes(x=signif_letters$row.names.hsd.groups.,y=45,label =signif_letters$hsd.groups.groups), color = "black", vjust=0, size = 8) +
  #labs(x ="Genotype", y = "Small Cell : Giant Cell Ratio") +
  theme(text = element_text(size = 20),axis.title.x=element_blank(),axis.title.y=element_blank(), legend.position = "none") +
  #Add Sample Sizes
  stat_n_text(size = 6)
#Exported as 6 inch by 6 inch pdf