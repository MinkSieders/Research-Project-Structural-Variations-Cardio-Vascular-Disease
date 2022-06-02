
---
title: "PCA analysis"
authors: "Mink Sieders, Daan van der Giezen"
date: "2022/05/08"
---
  
  
#1 preparation
  
  #1.1 import functions
  
library(pacman)
p_load(ggplot2, corrplot, tidyr, Hmisc, psych, reshape2, dplyr, tidyverse, rstatix, ggpubr, ggrepel, ape, vegan, ggfortify)


  #1.2 import input files

all_dsv_df   =  read.table('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Data\\lld_dsv.tsv',
                           header = TRUE,
                           sep = '\t',
                           fill = TRUE)

all_vsv_df   =  read.table('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Data\\lld_vsv.tsv',
                           header = TRUE,
                           sep = '\t',
                           fill = TRUE)

all_prot_df  =  read.table('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Data\\lld_prot.tsv',
                           header = TRUE,
                           sep = '\t',
                           fill = TRUE)



  #1.3 combine sv dataframes

d1 = all_dsv_df
d2 = all_vsv_df

d3 = cbind(d1, d2)


#2 statistical analysis of the data

  #2.1 code that calculates adonis R values of a list of species and a list of proteins

proteins_of_interest = c('FAS','PON3','Ep.CAM','IL.1RT2', 'GDF.15', 'CCL24')
species_of_interest = c('Blautia.wexlerae', 'Ruminococcus.lactaris', 'Blautia.obeum', 'Dorea.formicigenerans', 'Collinsella.sp', 'Holdemanella.biformis', 'Coprococcus.comes', 'Oscillibacter.sp')
R_value_table = data.frame()
R_value_table2 = data.frame()
for(i in 1:length(species_of_interest)){
  for(j in 1:length(proteins_of_interest)){
    
    sv = d3[ , grepl(species_of_interest[i], names(d3))]
    pm = as.matrix(all_prot_df[ , grepl(proteins_of_interest[j], names(all_prot_df))])

    pos1 = c()
    pos2 = c()
    for(k in 1:length(rownames(sv))){
      if(rowSums(is.na(sv[k,])) != length(colnames(sv))){
        pos1 = c(pos1, k)} }
    sv = sv[pos1,]
    pm = pm[pos1]
    for(l in 1:length(pm)){
      if(is.na(pm[l]) == FALSE){
        pos2 = c(pos2, l)} }
    sv = sv[pos2,]
    pm = pm[pos2]

    dm = vegdist(sv, method = "canberra", na.rm = TRUE) %>% as.matrix
    ad = adonis2(dm ~ pm)

    adr = data.frame(ad$R2[1], 
                     proteins_of_interest[j], 
                     gsub('.', ' ', species_of_interest[i], fixed = TRUE),
                     ad$`Pr(>F)`[1])

    R_value_table = rbind(R_value_table, adr)

  }
}

colnames(R_value_table) = c('R','prot','spec_name','p')

R_value_table = R_value_table %>%
  mutate(correlation = case_when(p <= .01 ~ '  ***',
                                 p <= .05 ~ '  **',
                                 p <= .1 ~ '  *',
                                 TRUE ~ ' '))

colnames(R_value_table) = c('R','prot','spec_name','p', 'significance_level')


  #2.3 PCA analysis between species and the best associated proteins gathered and from the table in '2.1'

filtered_data_sv = d3[ , grepl( 'Dorea.formicigenerans', names(d3))]
prot_mat = as.matrix(all_prot_df$Ep.CAM)

pos1 = c()
for(i in 1:length(rownames(filtered_data_sv))){
  if(rowSums(is.na(filtered_data_sv[i,])) != length(colnames(filtered_data_sv))){
    pos1 = c(pos1, i)} }

filtered_data_sv = filtered_data_sv[pos1,]
prot_mat = prot_mat[pos1]

pos2 = c()
for(i in 1:length(prot_mat)){
  if(is.na(prot_mat[i]) == FALSE){
    pos2 = c(pos2, i)} }

filtered_data_sv = filtered_data_sv[pos2,]
prot_mat = prot_mat[pos2]

dist_mat = vegdist(filtered_data_sv, method = "canberra", na.rm = TRUE) %>% as.matrix
adonis = adonis2(dist_mat ~ prot_mat)

PCA = pcoa(dist_mat)
PCA_vecs = PCA$vectors
PCA_vecs_protein_levels = cbind(PCA_vecs[, 1:2], prot_mat)
colnames(PCA_vecs_protein_levels) = c('PCA1','PCA2','prot_level')
eigen_titles_1 = paste('PCA1 (', round(PCA$values$Eigenvalues[1], digits = 2), '%)', sep = '') 
eigen_titles_2 = paste('PCA2 (', round(PCA$values$Eigenvalues[2], digits = 2), '%)', sep = '')
R2_PCA_plot = round(adonis$R2[1], digits = 4)


#3 visualization of the PCA analysis

  #3.1 barplot for adonis R values for protein and species (each combination is calculated using '2.1')

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\PCA barplot sv ~ prot.png', 
    height = 700, width = 800) 

ggplot(R_value_table, aes(fill = prot, y = R, x = spec_name)) + 
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(label = significance_level, angle = 90), position = position_dodge(0.9), vjust = .6, hjust = 0.1, size = 6) +
  scale_y_continuous(limits = c(0, 0.016)) +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 17),
        axis.text.x = element_text(size = 12, angle = 60, hjust = .95, vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=25)) +
  labs(fill = 'protein',
       y = 'variance explained by protein levels')
  
dev.off()


  #3.2 PCA plot for best associated proteins and species calculated in '2.3'

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\PCoA sv ~ prot.png', 
    height = 700, width = 700) 

ggplot(PCA_vecs_protein_levels, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = prot_level), size = 5, alpha = 0.8) +
  scale_color_distiller(palette = 'Spectral') +
  labs(x = eigen_titles_1,
       y = eigen_titles_2) +
  geom_vline(xintercept = 0, 
             linetype = 'dotted') +
  geom_hline(yintercept = 0, 
             linetype = 'dotted') +
  ggtitle(paste('R2 = ', R2_PCA_plot, ',             D. formicigenerans - Ep.CAM', sep = '')) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        title = element_text(size = 20),
        legend.key.height = unit(3, 'cm'),
        legend.key.width = unit(1.2, 'cm'),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 17)) +
  labs(colour = 'protein \n levels')
  
dev.off()


