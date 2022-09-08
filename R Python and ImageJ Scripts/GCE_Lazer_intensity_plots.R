### Figure 4 ###

#Load in data manually (pAR.intensity.csv)

Lazer_intensity = subset(pAR_lazer_intensity, Genotype %in% c("pAR111", "pAR254", "pAR257","pAR307","pLH166","pAR308"))

#Perform an Anova and Tukey HSD on the Small Cell data
Signif_lazer = aov(Lazer_intensity ~ Genotype, data = Lazer_intensity)
summary(Signif_lazer)
hsd = HSD.test(Signif_lazer, "Genotype")

#group data by genotype and extract the significance letters from HSD test
Lazer_summary = Lazer_intensity %>% group_by(Genotype) %>% summarize(Mean.ratio=mean(Lazer_intensity))
signif_letters = data.frame(row.names(hsd$groups),hsd$groups$groups)
#reorder the the significance letters by alphabetical order of the genotypes so the order is the same as the plot
signif_letters <- signif_letters[order(signif_letters$row.names.hsd.groups.),]

#Plot the Lazer_intensity
ggplot(data = Lazer_intensity, aes(x=reorder(Genotype, -Lazer_intensity), y=Lazer_intensity, color = Genotype)) + 
  geom_boxplot(lwd = 1.5) +
  theme_classic() +
  scale_color_manual(values = c("blue","purple","#01af50","#01af50","purple","purple")) +
  geom_text(data=Lazer_summary,aes(x=Genotype,y=11,label =signif_letters$hsd.groups.groups),vjust=0, color = "black", size = 8) +
  labs(y = "Laser power (%)") +
  theme(text = element_text(size = 20), axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") +
  stat_n_text(size = 6)
#Exported as 1350 x 900 Pixel Tif

### Figure 3 ###

Lazer_intensity3 = subset(pAR_lazer_intensity, Genotype %in% c("pAR111", "pAR254", "pAR257","pAR261"))

#Perform an Anova and Tukey HSD on the Small Cell data
Signif_lazer = aov(Lazer_intensity ~ Genotype, data = Lazer_intensity3)
summary(Signif_lazer)
hsd = HSD.test(Signif_lazer, "Genotype")

#group data by genotype and extract the significance letters from HSD test
Lazer_summary = Lazer_intensity3 %>% group_by(Genotype) %>% summarize(Mean.ratio=mean(Lazer_intensity))
signif_letters = data.frame(row.names(hsd$groups),hsd$groups$groups)
#reorder the the significance letters by alphabetical order of the genotypes so the order is the same as the plot
signif_letters <- signif_letters[order(signif_letters$row.names.hsd.groups.),]

#Plot the Lazer_intensity
ggplot(data = Lazer_intensity3, aes(x=reorder(Genotype, -Lazer_intensity), y=Lazer_intensity, color = Genotype)) + 
  geom_boxplot(lwd = 1.5) +
  theme_classic() +
  scale_color_manual(values = c("blue","purple","#01af50","purple")) +
  geom_text(data=Lazer_summary,aes(x=Genotype,y=11,label =signif_letters$hsd.groups.groups),vjust=0, color = "black", size = 8) +
  labs(y = "Laser power (%)", x = "Genotype") +
  theme(text = element_text(size = 20), legend.position = "none") +
  stat_n_text(size = 6)
#Export as a 760 x 1400 pixel Tif
