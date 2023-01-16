#manually load in the files: 
#Intensity_data_ATML1_crosses_Giantv3_Smallv4_parsed.csv

#Load libraries
library(dplyr)
library(ggplot2)
library(agricolae)
library(EnvStats)
#install.packages('EnvStats')

#Split dataset by cell type (giant and small cells)
x = split(Intensity_data_ATML1_crosses_Giantv3_Smallv4_parsed, Intensity_data_ATML1_crosses_Giantv3_Smallv4_parsed$Cell_type)
list2env(x, envir = .GlobalEnv)
#Now we have pAR254 Small cells (v4) in the df "Small", and pAR111 and pAR254 Giant
#Cells (v3) in the df "Giant". We still have to separate 111 and 254 giant cells.
Giant254 = subset(Giant, Genotype %in% c("pAR254xPDFflagATML1ox","pAR254xCol0","pAR254xatml14"))
Giant111 = subset(Giant, Genotype %in% c("pAR111xPDFflagATML1ox","pAR111xCol0","pAR111xatml14ko"))



#Perform ANOVA and Tukey HSD Statistical test on pAR254 small cell intensity
Signif_small = aov(Mean_intensity ~ Genotype, data = Small)
summary(Signif_small)
hsd = HSD.test(Signif_small, "Genotype")

#Try using t-tests instead
t.test(Mean_intensity ~ Genotype, data = subset(Small, Genotype %in% c("pAR254xatml14","pAR254xPDFflagATML1ox"))) #.73
t.test(Mean_intensity ~ Genotype, data = subset(Small, Genotype %in% c("pAR254xatml14","pAR254xCol0"))) #.08
t.test(Mean_intensity ~ Genotype, data = subset(Small, Genotype %in% c("pAR254xPDFflagATML1ox","pAR254xCol0"))) #.11

#Perform ANOVA and Tukey HSD Statistical test on pAR254 giant cell intensity
Signif_giant = aov(Mean_intensity ~ Genotype, data = Giant254)
summary(Signif_giant)
hsd2 = HSD.test(Signif_giant, "Genotype")

#Try using t-tests instead
t.test(Mean_intensity ~ Genotype, data = subset(Giant254, Genotype %in% c("pAR254xatml14","pAR254xPDFflagATML1ox"))) #1e-9
t.test(Mean_intensity ~ Genotype, data = subset(Giant254, Genotype %in% c("pAR254xatml14","pAR254xCol0"))) #8e-5
t.test(Mean_intensity ~ Genotype, data = subset(Giant254, Genotype %in% c("pAR254xPDFflagATML1ox","pAR254xCol0"))) #2e-5

#Perform ANOVA and Tukey HSD Statistical test on pAR111 giant cell intensity
Signif_giant = aov(Mean_intensity ~ Genotype, data = Giant111)
summary(Signif_giant)
hsd3 = HSD.test(Signif_giant, "Genotype")

#Try using t-tests instead
t.test(Mean_intensity ~ Genotype, data = subset(Giant111, Genotype %in% c("pAR111xatml14ko","pAR111xPDFflagATML1ox"))) #7e-5
t.test(Mean_intensity ~ Genotype, data = subset(Giant111, Genotype %in% c("pAR111xatml14ko","pAR111xCol0"))) #0.18
t.test(Mean_intensity ~ Genotype, data = subset(Giant111, Genotype %in% c("pAR111xPDFflagATML1ox","pAR111xCol0"))) #.0003

#group data by genotype and extract the significance letters from HSD test
Small_summary = Small %>% group_by(Genotype) %>% summarize(Mean.ratio=mean(Mean_intensity))
signif_letters_S = data.frame(row.names(hsd$groups),hsd$groups$groups)
#reorder the the significance letters by alphabetical order of the genotypes so the order is the same as the plot
signif_letters_S <- signif_letters_S[order(signif_letters_S$row.names.hsd.groups.),]

#group data by genotype and extract the significance letters from HSD test
Giant_summary254 = Giant254 %>% group_by(Genotype) %>% summarize(Mean.ratio=mean(Mean_intensity))
signif_letters_G254 = data.frame(row.names(hsd2$groups),hsd2$groups$groups)
#reorder the the significance letters by alphabetical order of the genotypes so the order is the same as the plot
signif_letters_G254 <- signif_letters_G254[order(signif_letters_G254$row.names.hsd2.groups.),]

#group data by genotype and extract the significance letters from HSD test
Giant_summary111 = Giant111 %>% group_by(Genotype) %>% summarize(Mean.ratio=mean(Mean_intensity))
signif_letters_G111 = data.frame(row.names(hsd3$groups),hsd3$groups$groups)
#reorder the the significance letters by alphabetical order of the genotypes so the order is the same as the plot
signif_letters_G111 <- signif_letters_G111[order(signif_letters_G111$row.names.hsd3.groups.),]

#Plot pAR254 Small Cell Intensity
ggplot(Small, aes(x=Genotype, y=Mean_intensity, color = Genotype)) + 
  geom_boxplot(lwd = 1.5) +
  theme_classic() +
  scale_color_manual( values = c("green", "light green", "dark green")) +
  scale_x_discrete(labels = c('pAR254 X atml1-4','pAR254 X WT','pAR254 X PDF-flag::ATML1 OX')) +
  geom_text(data=Small_summary,aes(x=Genotype,y=35,label=signif_letters_S$hsd.groups.groups),vjust=0, color = "black", size = 8) +
  labs(title="pAR254 Small Cell Intensity") +
  theme(text = element_text(size = 20),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
  stat_n_text(size = 6)

#Plot pAR254 Giant Cell Intensity
ggplot(Giant254, aes(x=Genotype, y=Mean_intensity, color = Genotype)) + 
  geom_boxplot(lwd = 1.5) +
  theme_classic() +
  scale_color_manual( values = c("green", "light green", "dark green")) +
  scale_x_discrete(labels = c('pAR254 X atml1-4','pAR254 X WT','pAR254 X PDF-flag::ATML1 OX')) +
  geom_text(data=Giant_summary254,aes(x=Genotype,y=80,label =signif_letters_G254$hsd2.groups.groups),vjust=0, color = "black", size = 8) +
  #labs(title="pAR254 Giant Cell Intensity") +
  theme(text = element_text(size = 20), axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
  stat_n_text(size = 6)
#Export as 360x525 pixel tif
#Export as 6 inch by 6 inch PDF

#Plot pAR111 Giant Cell Intensity
ggplot(Giant111, aes(x=Genotype, y=Mean_intensity, color = Genotype)) + 
  geom_boxplot(lwd = 1.5) +
  theme_classic() +
  scale_color_manual( values = c("green", "light green", "dark green")) +
  scale_x_discrete(labels = c('pAR111 X atml1-4','pAR111 X WT','pAR111 X PDF-flag::ATML1 OX')) +
  geom_text(data=Giant_summary111,aes(x=Genotype,y=80,label =signif_letters_G111$hsd3.groups.groups),vjust=0, color = "black", size = 8) +
  #labs(title="pAR111 Giant Cell Intensity") +
  theme(text = element_text(size = 20), axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
  stat_n_text(size = 6)
#Export as 360x525 pixel tif




############### ALL NUCLEI ANALYSIS ###################
#load Intensity_data_ATML1_crosses_Giantv3_Smallv4_allNuclei.csv

#Split dataset by cell type (giant and small cells)
x = split(Intensity_data_ATML1_crosses_Giantv3_Smallv4_allNuclei, Intensity_data_ATML1_crosses_Giantv3_Smallv4_allNuclei$cell_type)
list2env(x, envir = .GlobalEnv)
#Now we have pAR254 Small cells (v4) in the df "Small", and pAR111 and pAR254 Giant
#Cells (v3) in the df "Giant". We still have to separate 111 and 254 giant cells.
Giant254 = subset(Giant, genotype %in% c("pAR254xOX","pAR254xCol","pAR254xKO"))
Giant111 = subset(Giant, genotype %in% c("pAR111xOX","pAR111xCol","pAR111xKO"))

##assess normality and transform data
#install.packages("ggpubr")
library(ggpubr)
ggqqplot(1/(Giant111$Mean))
ggqqplot(1/subset(Giant111, genotype %in% c("pAR111xOX"))$Mean)
ggqqplot(1/subset(Giant111, genotype %in% c("pAR111xCol"))$Mean)
ggqqplot(1/subset(Giant111, genotype %in% c("pAR111xKO"))$Mean)

ggqqplot(1/subset(Giant254, genotype %in% c("pAR254xOX"))$Mean)
ggqqplot(1/subset(Giant254, genotype %in% c("pAR254xCol"))$Mean)
ggqqplot(1/subset(Giant254, genotype %in% c("pAR254xKO"))$Mean)


#Perform ANOVA and Tukey HSD Statistical test on pAR254 small cell intensity
Signif_small = aov(1/Mean ~ genotype, data = Small)
summary(Signif_small)
hsd = HSD.test(Signif_small, "genotype")

#Perform ANOVA and Tukey HSD Statistical test on pAR254 giant cell intensity
Signif_giant = aov(1/Mean ~ genotype, data = Giant254)
summary(Signif_giant)
hsd2 = HSD.test(Signif_giant, "genotype")

#Perform ANOVA and Tukey HSD Statistical test on pAR111 giant cell intensity
Signif_giant = aov(1/Mean ~ genotype, data = Giant111)
summary(Signif_giant)
hsd3 = HSD.test(Signif_giant, "genotype")


#group data by genotype and extract the significance letters from HSD test
Small_summary = Small %>% group_by(genotype) %>% summarize(Mean.ratio=mean(Mean))
signif_letters_S = data.frame(row.names(hsd$groups),hsd$groups$groups)
#reorder the the significance letters by alphabetical order of the genotypes so the order is the same as the plot
signif_letters_S <- signif_letters_S[order(signif_letters_S$row.names.hsd.groups.),]

#group data by genotype and extract the significance letters from HSD test
Giant_summary254 = Giant254 %>% group_by(genotype) %>% summarize(Mean.ratio=mean(Mean))
signif_letters_G254 = data.frame(row.names(hsd2$groups),hsd2$groups$groups)
#reorder the the significance letters by alphabetical order of the genotypes so the order is the same as the plot
signif_letters_G254 <- signif_letters_G254[order(signif_letters_G254$row.names.hsd2.groups.),]

#group data by genotype and extract the significance letters from HSD test
Giant_summary111 = Giant111 %>% group_by(genotype) %>% summarize(Mean.ratio=mean(Mean))
signif_letters_G111 = data.frame(row.names(hsd3$groups),hsd3$groups$groups)
#reorder the the significance letters by alphabetical order of the genotypes so the order is the same as the plot
signif_letters_G111 <- signif_letters_G111[order(signif_letters_G111$row.names.hsd3.groups.),]


#Plot pAR254 Small Cell Intensity
ggplot(Small, aes(x=reorder(genotype,Mean), y=Mean, color = genotype)) + 
  geom_boxplot(lwd = 1) +
  theme_classic() +
  scale_color_manual( values = c("green", "light green", "dark green")) +
  scale_x_discrete(labels = c('pAR254 X atml1-4','pAR254 X WT','pAR254 X PDF-flag::ATML1 OX')) +
  geom_text(data=Small_summary,aes(x=genotype,y=35,label=signif_letters_S$hsd.groups.groups),vjust=0, color = "black", size = 8) +
  labs(title="pAR254 Small Cell Intensity") +
  theme(text = element_text(size = 20),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
  stat_n_text(size = 6)

#Plot pAR254 Giant Cell Intensity
ggplot(Giant254, aes(x=reorder(genotype,Mean), y=Mean, fill = genotype)) + 
  geom_violin() +
  geom_boxplot(width=0.07, outlier.shape = NA) +
  theme_classic() +
  scale_fill_manual( values = c("green", "light green", "dark green")) +
  scale_x_discrete(labels = c('','','')) +
  geom_text(data=Giant_summary254,aes(x=genotype,y=250,label =signif_letters_G254$hsd2.groups.groups),vjust=0, color = "black", size = 8) +
  #labs(title="pAR254 Giant Cell Intensity") +
  theme(text = element_text(size = 20), axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
  stat_n_text(size = 6)
#Export as PDF

#Plot pAR111 Giant Cell Intensity
ggplot(Giant111, aes(x=reorder(genotype,Mean), y=Mean, fill = genotype)) + 
  geom_violin() + #, outlier.shape = NA
  geom_boxplot(width=0.07, outlier.shape = NA) +
  #geom_jitter(size = 0.2, width = 0.2) +
  #geom_point(data = function(x) dplyr::filter_(x, ~ outlier), position = 'jitter') +
  theme_classic() +
  scale_fill_manual( values = c("green", "light green", "dark green")) +
  scale_x_discrete(labels = c('pAR111 X atml1-4','pAR111 X WT','pAR111 X PDF-flag::ATML1 OX')) +
  geom_text(data=Giant_summary111,aes(x=genotype,y=250,label =signif_letters_G111$hsd3.groups.groups),vjust=0, color = "black", size = 8) +
  #labs(title="pAR111 Giant Cell Intensity") +
  theme(text = element_text(size = 20), axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
  stat_n_text(size = 6)
#Export as pdf





##### Make Histograms

#plottings all of the giant cells across all genotypes w/pAR254 together
ggplot(Giant254, aes(x=Mean)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=mean(Mean)),
             color="blue", linetype="dashed", size=1)

ggplot(Giant254, aes(x=Mean, color=genotype)) +
  geom_histogram(fill="white", position="dodge")+
  #geom_vline(data=Giant254, aes(xintercept=grp.mean, color=genotype),
   #          linetype="dashed")+
  theme(legend.position="top")

#install.packages('ggridges')
library(ggridges)
theme_set(theme_minimal())

ggplot(Giant111, aes(x = Mean, y = genotype, group = genotype, fill = stat(x))) + 
  #geom_density_ridges(fill = "blue", rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 0.5)
  geom_density_ridges_gradient(scale = 1, size = 0.3, rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_gradient(name = "A.U", low = "white", high = "limegreen")

ggplot(Giant254, aes(x = Mean, y = genotype, group = genotype, fill = stat(x))) + 
  #geom_density_ridges(fill = "blue", rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 0.5)
  geom_density_ridges_gradient(scale = 1, size = 0.3, rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_gradient(name = "A.U", low = "white", high = "limegreen")

ggplot(Small, aes(x = Mean, y = genotype, group = genotype, fill = stat(x))) + 
  #geom_density_ridges(fill = "blue", rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 0.5)
  geom_density_ridges_gradient(scale = 1, size = 0.3, rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 0.5) +
  scale_fill_gradient(name = "A.U", low = "white", high = "limegreen")
