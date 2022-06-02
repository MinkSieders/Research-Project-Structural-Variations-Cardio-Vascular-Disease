
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

all_prot_df  =  read.table('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Data\\lld_prot.tsv',
                           header = TRUE,
                           sep = '\t',
                           fill = TRUE)

all_cvd_df   =  read.table('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Data\\lld_cvd.tsv',
                           header = TRUE,
                           sep = '\t',
                           fill = TRUE)


  #1.2 change 1's and 0's to male and female and filter for 1 specific gender

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
all_prot_df = all_prot_df[filter_pos,]


#2 statistical analysis of the data

  #2.1 PCA analysis between overall protein level makeup and disease risk

prot_mat = all_prot_df
cvd_mat = all_cvd_df
pos1 = c()
for(i in 1:length(rownames(prot_mat))){
  if(rowSums(is.na(prot_mat[i,])) != length(colnames(prot_mat))){
    pos1 = c(pos1, i)} }
prot_mat = prot_mat[pos1,]
cvd_mat = cvd_mat[pos1,]
pos2 = c()
for(i in 1:length(rownames(cvd_mat))){
  if(is.na(cvd_mat[i,1]) == FALSE){
    pos2 = c(pos2, i)} }
prot_mat = prot_mat[pos2,]
cvd_mat = cvd_mat[pos2,]

dist_mat = vegdist(prot_mat, method = "canberra", na.rm = TRUE) %>% as.matrix
adonis = adonis2(dist_mat ~ as.numeric(cvd_mat[,1]))

PCA = pcoa(dist_mat)
PCA_vecs = PCA$vectors
PCA_vecs_disease_risk = cbind(PCA_vecs[, 1:2], cvd_mat[1])
colnames(PCA_vecs_disease_risk) = c('PCA1','PCA2','disease_risk')
eigen_titles_1 = paste('PCA1 (', round(PCA$values$Eigenvalues[1], digits = 2), '%)', sep = '') 
eigen_titles_2 = paste('PCA2 (', round(PCA$values$Eigenvalues[2], digits = 2), '%)', sep = '')
R2_PCA_plot = round(adonis$R2[1], digits = 4)
p_PCA_plot = adonis$`Pr(>F)`[1]

#3 visualization of the PCA analysis

  #3.1 PCA plot for protein profile and disease risk

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\PCoA prot ~ cvd.png', 
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
  ggtitle(paste('R2 = ', R2_PCA_plot, ', p < 0.002', sep = '')) +
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

