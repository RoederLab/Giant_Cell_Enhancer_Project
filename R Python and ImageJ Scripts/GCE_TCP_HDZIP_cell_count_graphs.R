### TCP and Hdzip Cell Count Boxplots ###
### Load in the dataset GCE_Image_processing_data_pBRs_jawcrosses_T6_parsed.csv

#rename df

TCP_hdzip_counts = `GCE_Image_processing_data_pBRs_jawcrosses_T6_parsed`

#Subset the data by cell type (one df for small cells, one df for giant cells)
x = split(TCP_hdzip_counts, TCP_hdzip_counts$Cell_type)
list2env(x, envir = .GlobalEnv)

#Load in required libraries
library(dplyr)
library(ggplot2)
library(agricolae)
library(EnvStats)

#Subset the genotypes of interest based on which TF is targeted
Giant_sub_TCP = subset(Giant, Genotype %in% c("pAR111", "pAR254", "pAR111xCol0","pAR254xCol0","pAR111xjaw1d","pAR254xjaw1d","pBR63","pBR67"))
Small_sub_TCP = subset(Small, Genotype %in% c("pAR111", "pAR254", "pAR111xCol0","pAR254xCol0","pAR111xjaw1d","pAR254xjaw1d","pBR63","pBR67"))

Giant_sub_HDzip = subset(Giant, Genotype %in% c("pAR111", "pAR254","pBR65","pBR69"))
Small_sub_HDzip = subset(Small, Genotype %in% c("pAR111", "pAR254","pBR65","pBR69"))


#Perform an Anova and Tukey HSD on the Small Cell data for TCPs
Signif_small = aov(Count ~ Genotype, data = Small_sub_TCP)
summary(Signif_small)
hsd = HSD.test(Signif_small, "Genotype")

#Perform an Anova and Tukey HSD on the Giant Cell data for TCPs
Signif_giant = aov(Count ~ Genotype, data = Giant_sub_TCP)
summary(Signif_giant)
hsd2 = HSD.test(Signif_giant, "Genotype")

#group data by genotype and extract the significance letters from HSD test
Small_summary = Small_sub_TCP %>% group_by(Genotype) %>% summarize(Mean.ratio=mean(Count))
signif_letters_S = data.frame(row.names(hsd$groups),hsd$groups$groups)
#reorder the the significance letters by alphabetical order of the genotypes so the order is the same as the plot
signif_letters_S <- signif_letters_S[order(signif_letters_S$row.names.hsd.groups.),]

#group data by genotype and extract the significance letters from HSD test
Giant_summary = Giant_sub_TCP %>% group_by(Genotype) %>% summarize(Mean.ratio=mean(Count))
signif_letters_G = data.frame(row.names(hsd2$groups),hsd2$groups$groups)
#reorder the the significance letters by alphabetical order of the genotypes so the order is the same as the plot
signif_letters_G <- signif_letters_G[order(signif_letters_G$row.names.hsd2.groups.),]

#Rearrange the boxplots in an order that makes more sense
Giant_sub_TCP %>%
  arrange(Genotype) %>% 
  mutate(Genotype = factor(Genotype, levels = c("pAR111","pAR111xCol0","pAR111xjaw1d","pBR63","pAR254","pAR254xCol0","pAR254xjaw1d","pBR67"))) %>%

#Plot the giant cell counts
ggplot(aes(x=Genotype, y=Count, color = Genotype)) + 
  geom_boxplot(lwd = 1.5) +
  theme_classic() +
  scale_color_manual(values = c("blue","blue", "red","red","purple","purple", "red","red")) +
  geom_text(data=Giant_summary,aes(x=Genotype,y=67,label =signif_letters_G$hsd2.groups.groups),vjust=0, color = "black", size = 8) +
  labs(x ="Genotype", y = "Giant Cell Count") +
  theme(text = element_text(size = 20), legend.position = "none") +
  stat_n_text(size = 6)
#Export as a 1400 x 1200 Tif

#Rearrange the boxplots in an order that makes more sense
Small_sub_TCP %>%
  arrange(Genotype) %>% 
  mutate(Genotype = factor(Genotype, levels = c("pAR111","pAR111xCol0","pAR111xjaw1d","pBR63","pAR254","pAR254xCol0","pAR254xjaw1d","pBR67"))) %>%

#plot the small cell counts
ggplot(aes(x=Genotype, y=Count, color = Genotype)) + 
  geom_boxplot(lwd = 1.5) +
  theme_classic() +
  scale_color_manual(values = c("blue","blue", "red","red","purple","purple", "red","red")) +
  geom_text(data=Small_summary,aes(x=Genotype,y=1100,label =signif_letters_S$hsd.groups.groups),vjust=0, color = "black", size = 8) +
  labs(x ="Genotype", y = "Small Cell Count") +
  theme(text = element_text(size = 20), legend.position = "none") +
  stat_n_text(size = 6)
#Export as a 1400 x 1200 pixel tif

###  Figure 6 ####

#Perform an Anova and Tukey HSD on the Small Cell data for HDzips
Signif_small = aov(Count ~ Genotype, data = Small_sub_HDzip)
summary(Signif_small)
hsd = HSD.test(Signif_small, "Genotype")

#Perform an Anova and Tukey HSD on the Giant Cell data for HDzips
Signif_giant = aov(Count ~ Genotype, data = Giant_sub_HDzip)
summary(Signif_giant)
hsd2 = HSD.test(Signif_giant, "Genotype")

#group data by genotype and extract the significance letters from HSD test
Small_summary = Small_sub_HDzip %>% group_by(Genotype) %>% summarize(Mean.ratio=mean(Count))
signif_letters_S = data.frame(row.names(hsd$groups),hsd$groups$groups)
#reorder the the significance letters by alphabetical order of the genotypes so the order is the same as the plot
signif_letters_S <- signif_letters_S[order(signif_letters_S$row.names.hsd.groups.),]

#group data by genotype and extract the significance letters from HSD test
Giant_summary = Giant_sub_HDzip %>% group_by(Genotype) %>% summarize(Mean.ratio=mean(Count))
signif_letters_G = data.frame(row.names(hsd2$groups),hsd2$groups$groups)
#reorder the the significance letters by alphabetical order of the genotypes so the order is the same as the plot
signif_letters_G <- signif_letters_G[order(signif_letters_G$row.names.hsd2.groups.),]

#Plot the giant cell counts
ggplot(subset(Giant, Genotype %in% c("pAR111", "pAR254", "pBR65","pBR69")), aes(x=Genotype, y=Count, color = Genotype)) + 
  geom_boxplot() +
  theme_classic() +
  scale_color_manual(values = c("blue","purple", "red","red")) +
  geom_text(data=Giant_summary,aes(x=Genotype,y=45,label =signif_letters_G$hsd2.groups.groups),vjust=0, color = "black", size = 8) +
  labs(title="Giant Cell Count",
       x ="Genotype", y = "Count") +
  theme(text = element_text(size = 20), legend.position = "none") +
  stat_n_text(size = 6)
#Export as 400x525 pixel tif
#Export as 6 inch by 6 inch PDF

#plot the small cell counts
ggplot(subset(Small, Genotype %in% c("pAR111", "pAR254", "pBR65","pBR69")), aes(x=Genotype, y=Count, color = Genotype)) + 
  geom_boxplot() +
  theme_classic() +
  scale_color_manual(values = c("blue","purple", "red","red")) +
  geom_text(data=Small_summary,aes(x=Genotype,y=1100,label =signif_letters_S$hsd.groups.groups),vjust=0, color = "black", size = 8) +
  labs(title="Small Cell Count",
       x ="Genotype", y = "Count") +
  theme(text = element_text(size = 20), legend.position = "none") +
  stat_n_text(size = 6)
#Export as 6 inch by 6 inch PDF

