
---
title: "PCA CVD analysis"
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

all_cvd_df   =  read.table('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Data\\lld_cvd.tsv',
                           header = TRUE,
                           sep = '\t',
                           fill = TRUE)


  #1.3 combine sv dataframes

d1 = all_dsv_df
d2 = all_vsv_df

d3 = cbind(d1, d2)
  
  
  #1.4 change 1's and 0's to male and female and filter for 1 specific gender

all_cvd_df = all_cvd_df %>%
  mutate(Sex = case_when(Sex == 1 ~ 'male',
                         Sex == 0 ~ 'female'))

analyzing = 'male' # choose either 'male' or 'female' to do analysis

filter_pos = c()
for(i in 1:length(rownames(all_cvd_df))){
  if(all_cvd_df$Sex[i] == analyzing){
    filter_pos = c(filter_pos, i)
  }
}

all_cvd_df = all_cvd_df[filter_pos,]
d3 = d3[filter_pos,]


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
    cm = as.matrix(all_cvd_df)
    
    pos3 = c()
    pos1 = c()
    pos2 = c()
    for(k in 1:length(rownames(sv))){
      if(rowSums(is.na(sv[k,])) != length(colnames(sv))){
        pos1 = c(pos1, k)} }
    sv = sv[pos1,]
    pm = pm[pos1]
    cm = cm[pos1,]
    for(l in 1:length(pm)){
      if(is.na(pm[l]) == FALSE){
        pos2 = c(pos2, l)} }
    sv = sv[pos2,]
    pm = pm[pos2]
    cm = cm[pos2,]
    for(m in 1:length(rownames(cm))){
      if(is.na(cm[m, 1]) == FALSE){
        pos3 = c(pos3, m)} }
    sv = sv[pos3,]
    pm = pm[pos3]
    cm = cm[pos3,]
    
    dm = vegdist(sv, method = "canberra", na.rm = TRUE) %>% as.matrix
    ad = adonis2(dm ~ pm + as.numeric(cm[,1]))

    adr = data.frame(ad$R2[1], 
                     proteins_of_interest[j], 
                     gsub('.', ' ', species_of_interest[i], fixed = TRUE),
                     ad$`Pr(>F)`[1])
    adr2 = data.frame(ad$R2[2], 
                     proteins_of_interest[j], 
                     gsub('.', ' ', species_of_interest[i], fixed = TRUE),
                     ad$`Pr(>F)`[2])
    R_value_table = rbind(R_value_table, adr)
    R_value_table2 = rbind(R_value_table2, adr2)
  }
}

colnames(R_value_table) = c('R','prot','spec_name','p')
R_value_table['factor'] = 'protein_level'
colnames(R_value_table2) = c('R','prot','spec_name','p')
R_value_table2['factor'] = 'disease_risk'

Com_R_val_tab = rbind(R_value_table, R_value_table2)

Com_R_val_tab = Com_R_val_tab %>%
  mutate(correlation = case_when(p <= .01 ~ '***',
                                 p <= .05 ~ '**',
                                 p <= .1 ~ '*',
                                 TRUE ~ ' '))

save(Com_R_val_tab, file = 'C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Results\\male sv~cvd+prot.RData')


  #2.3 PCA analysis between species, disease risk, and the best associated proteins gathered and from the table in '2.1'

filtered_data_sv = d3[ , grepl( 'Blautia.wexlerae', names(d3))]
prot_mat = as.matrix(all_prot_df$CCL24)
cvd_mat = all_cvd_df

pos1 = c()
for(i in 1:length(rownames(filtered_data_sv))){
  if(rowSums(is.na(filtered_data_sv[i,])) != length(colnames(filtered_data_sv))){
    pos1 = c(pos1, i)} }

filtered_data_sv = filtered_data_sv[pos1,]
prot_mat = prot_mat[pos1]
cvd_mat = cvd_mat[pos1,]

pos2 = c()
for(i in 1:length(prot_mat)){
  if(is.na(prot_mat[i]) == FALSE){
    pos2 = c(pos2, i)} }

filtered_data_sv = filtered_data_sv[pos2,]
prot_mat = prot_mat[pos2]
cvd_mat = cvd_mat[pos2,]

pos3 = c()
for(i in 1:length(rownames(cvd_mat))){
  if(is.na(cvd_mat[i,1]) == FALSE){
    pos3 = c(pos3, i)} }

filtered_data_sv = filtered_data_sv[pos3,]
prot_mat = prot_mat[pos3]
cvd_mat = cvd_mat[pos3,]

dist_mat = vegdist(filtered_data_sv, method = "canberra", na.rm = TRUE) %>% as.matrix
adonis = adonis2(dist_mat ~ prot_mat + as.numeric(cvd_mat[,1]))

PCA = pcoa(dist_mat)
PCA_vecs = PCA$vectors
PCA_vecs_protein_levels = cbind(PCA_vecs[, 1:2], prot_mat)
colnames(PCA_vecs_protein_levels) = c('PCA1','PCA2','prot_level')
eigen_titles_1 = paste('PCA1 (', round(PCA$values$Eigenvalues[1], digits = 2), '%)', sep = '') 
eigen_titles_2 = paste('PCA2 (', round(PCA$values$Eigenvalues[2], digits = 2), '%)', sep = '')
R2_PCA_plot = round(adonis$R2[1], digits = 4)

R2_PCA_plot2 = round(adonis$R2[2], digits = 4)
PCA_vecs_disease_risk = cbind(PCA_vecs[, 1:2], cvd_mat[, 1])
colnames(PCA_vecs_disease_risk) = c('PCA1','PCA2','disease_risk')


#3 visualization of the PCA analysis

  #3.1 barplot for adonis R values for protein and species (each combination is calculated using '2.1')

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\PCA barplot.png', 
    height = 700, width = 800) 

ggplot(Com_R_val_tab,                         
       aes(x = spec_name,
           y = R,
           fill = factor)) + 
  geom_bar(stat = "identity",
           position = "stack") +
  scale_fill_manual(values = c("lightblue3", "darkcyan")) +
  geom_text(aes(label = correlation, angle = 90), position = position_stack(vjust = 0.5), hjust = .8, size = 6) +
  facet_grid(~ prot) +
  labs(y = 'variance explained by 2 factors') +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 17),
        axis.text.x = element_text(size = 12, angle = 60, hjust = .95, vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20))
  
dev.off()


  #3.2 PCA plot for best associated proteins and species calculated in '2.1' - disease risk

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\PCoA1.png', 
    height = 700, width = 700) 

ggplot(PCA_vecs_disease_risk, aes(x = PCA1, y = PCA2)) +
  geom_point(aes(colour = disease_risk), size = 5, alpha = 0.8) +
  scale_color_distiller(palette = 'Spectral') +
  labs(x = eigen_titles_1,
       y = eigen_titles_2) +
  geom_vline(xintercept = 0, 
             linetype = 'dotted') +
  geom_hline(yintercept = 0, 
             linetype = 'dotted') +
  ggtitle(paste('R2 = ', R2_PCA_plot2, ',             B. wexlerae - CCL24', sep = '')) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        title = element_text(size = 20),
        legend.key.height = unit(3, 'cm'),
        legend.key.width = unit(1.2, 'cm'),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 17)) +
  labs(colour = 'disease \n risk')

dev.off()


#3.3 PCA plot for best associated proteins and species calculated in '2.1' - protein_level

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\PCoA2.png', 
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
  ggtitle(paste('R2 = ', R2_PCA_plot, ',             B. wexlerae - CCL24', sep = '')) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        title = element_text(size = 20),
        legend.key.height = unit(3, 'cm'),
        legend.key.width = unit(1.2, 'cm'),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 17)) +
  labs(colour = 'protein \n level')

dev.off()
