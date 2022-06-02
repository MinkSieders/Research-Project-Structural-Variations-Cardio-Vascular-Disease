
---
title: "Statistics dSV data"
authors: "Mink Sieders, Daan van der Giezen"
date: "2022/04/23"
---
  
  
#1 preparation

  #1.1 import functions
  
library(pacman)
p_load(ggplot2, ggpubr, rstatix, gtools, corrplot, ggcorrplot, tidyr, Hmisc, psych, reshape2, dplyr, tidyverse, rstatix, ggpubr, ggrepel)


#1.2 import input files

all_dsv_df   =  read.table('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Data\\lld_dsv.tsv',
                           header = TRUE,
                           sep = '\t',
                           fill = TRUE)

all_prot_df   =  read.table('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Data\\lld_prot.tsv',
                            header = TRUE,
                            sep = '\t',
                            fill = TRUE)

d1 = all_dsv_df
d2 = all_prot_df


#2 statistical analysis of the data

  #2.1 correlation analysis between dSV and protein levels

dsv_corr_data = data.frame()  
  
for(i in 1:length(colnames(d2))){
    for(j in 1:length(colnames(d1))){
    
        if(length(which(all_dsv_df[, j] == 0)) >= 10 & length(which(all_dsv_df[, j] == 1)) >= 10){
        pos_0 = which(d1[, j] == 0)
        pos_1 = which(d1[, j] == 1)
    
        protein_level_0 = d2[pos_0, i]
        protein_level_1 = d2[pos_1, i]  
      
        corr_test = wilcox.test(protein_level_0, protein_level_1, paired = FALSE)
        corr_data = data.frame(colnames(d1)[j], colnames(d2)[i], length(c(protein_level_0, protein_level_1)), length(protein_level_0), 
                               length(protein_level_1), corr_test$p.value, (mean(protein_level_0, na.rm = TRUE) - 
                               mean(protein_level_1, na.rm = TRUE)))
        
        dsv_corr_data = rbind(dsv_corr_data, corr_data)
        
        }else{
        corr_data = data.frame(colnames(d1)[j], colnames(d2)[i], length(c(protein_level_0, protein_level_1)), length(protein_level_0), 
                               length(protein_level_1), NA, NA)
        names(corr_data)[names(corr_data) == 'NA.'] = colnames(dsv_corr_data)[6]
        names(corr_data)[names(corr_data) == 'NA..1'] = colnames(dsv_corr_data)[7]
        
        dsv_corr_data = rbind(dsv_corr_data, corr_data)}
    } 
  print(i)
  print(j)
}

p_adjust = p.adjust(dsv_corr_data[, 6], method = 'BH')
dsv_corr_data = cbind(dsv_corr_data, p_adjust)

colnames(dsv_corr_data) = c('dsv','protein', 'n', 'n1', 'n0', 'p', 'mean_difference', 'FDR')


  #2.2 save statistical analysis file containing the P values for later usage

save(dsv_corr_data, file = 'C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Results\\dsv_stat_anal_data.RData')


#3 visualisation

  #3.1 get the file with the statistical analysis data performed in '2'

dsv_corr_data = get(load('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Results\\dsv_stat_anal_data.RData'))
  

  #3.2 calculate absolute log2() values of the fold change

n_corr_data = dsv_corr_data

  #3.3 gathering the corresponding species names and making a list of top 3 positive and negative associations

n_corr_data = n_corr_data %>%
  mutate(effect = case_when(mean_difference >= .1 & p <= 0.05 ~ 'positive',
                            mean_difference <= -.1 & p <= 0.05 ~ 'negative',
                            TRUE ~ 'not significant'))

n_corr_data = n_corr_data[order(n_corr_data[,6],decreasing=FALSE),]

n_corr_data['name'] = NA

pos_rel = c()
neg_count = 0 
labs_rel = c()
for(i in 1:1000){
  x = sapply(str_split(n_corr_data[i, 1], pattern = '[.]'), `[`, 1)
    y = sapply(str_split(n_corr_data[i, 1], pattern = '[.]'), `[`, 2)
      name = paste(x, y, sep = ' ')
  if(neg_count >= 3){
    break
  }
  if(n_corr_data$effect[i] == 'negative'){ 
     pos_rel = c(pos_rel, i) 
     neg_count = neg_count + 1
     labs_rel = c(labs_rel, paste(n_corr_data$protein[i], name, sep = ' - ')) } }

n_corr_data[pos_rel, 10] = labs_rel

labs_rel2 = c()
pos_rel2 = c()
pos_count = 0
for(i in 1:1000){
  x = sapply(str_split(n_corr_data[i, 1], pattern = '[.]'), `[`, 1)
    y = sapply(str_split(n_corr_data[i, 1], pattern = '[.]'), `[`, 2)
      name = paste(x, y, sep = ' ')
  if(pos_count == 3){
    break
  }
  if(n_corr_data$effect[i] == 'positive'){
    pos_count = pos_count + 1
    pos_rel2 = c(pos_rel2, i)
    labs_rel2 = c(labs_rel2, paste(n_corr_data$protein[i], name, sep = ' - '))
  }
}

n_corr_data[pos_rel2, 10] = labs_rel2


  #3.4 volcano plot - visualisation of the correlation analysis between dSV and protein levels

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\dSV p-value volcano map2.png', 
    height = 2000, width = 1500) 

line_pos = c(-0.8234977, -.1, .1, 0.2143361)

ggplot(n_corr_data, 
       aes(x = mean_difference, 
           y = -log10(p), 
           fill = effect,    
           size = effect,
           alpha = effect)) + 
  geom_point(shape = 21, colour = 'black') + 
  geom_vline(xintercept = line_pos, 
             linetype = 'dotdash') +
  geom_hline(yintercept = -log10(0.05), 
             linetype = 'dashed') +
  geom_hline(yintercept = -log10(1e-5), 
             linetype = 'dotted') +
  geom_label_repel(show.legend = FALSE, aes(label = name, fontface = 'bold'), size = 7,
                   box.padding = unit(2, 'lines')) +
  scale_x_continuous(breaks = line_pos, labels = c('-0.82', '-0.1', '0.1', '0.21')) +
  scale_y_continuous(breaks = c(-log10(1e-5), -log10(0.05)), labels = c('p < 1e-5', 'p < 0.05')) +
  scale_fill_manual(values = c('positive' = 'chocolate1', 'negative' = 'darkcyan', 'not significant' = 'grey')) + 
  scale_size_manual(values = c('positive' = 4, 'negative' = 4, 'not significant' = 2)) + 
  scale_alpha_manual(values = c('positive' = .8, 'negative' = .8, 'not significant' = .4)) +
  labs(x = 'mean difference',
       y = '-log10(p)') +  
  guides(fill = guide_legend('correlation', override.aes = list(size=10)), size = 'none', alpha = 'none') + 
  theme_bw() +   
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.key.size = unit(4, 'cm'),
        legend.key.height = unit(4, 'cm'),
        legend.key.width = unit(4, 'cm'), 
        legend.title = element_text(size=30, face = 'bold'), 
        legend.text = element_text(size=18),
        plot.title = element_text(size=40, face = 'bold'),
        axis.text = element_text(size=20),
        axis.title = element_text(size=40))

dev.off()


  #3.5 heatmap of the r-value - visualisation of the correlation analysis between vSV and protein levels

t200_ncd = n_corr_data

matrix_r = pivot_wider(t200_ncd, id_cols = dsv, names_from = protein, values_from = mean_difference)
matrix_r = column_to_rownames(matrix_r, 'dsv')

matrix_FDR = pivot_wider(t200_ncd, id_cols = dsv, names_from = protein, values_from = FDR)
matrix_FDR = column_to_rownames(matrix_FDR, 'dsv')

matrix_p = pivot_wider(t200_ncd, id_cols = dsv, names_from = protein, values_from = p)
matrix_p = column_to_rownames(matrix_p, 'dsv')

col_pos = c()
for(i in 1:length(colnames(matrix_p))){
  if(length(which(matrix_p[, i] < 0.00006))){
    col_pos = c(col_pos, i)} }

row_pos = c()
for(i in 1:length(rownames(matrix_p))){
  if(length(which(matrix_p[i, ] < 0.00004))){
    row_pos = c(row_pos, i)} }

fm_FDR = as.matrix(matrix_FDR[row_pos, col_pos])
fm_r = as.matrix(matrix_r[row_pos, col_pos])  

fm_FDR_flip = fm_FDR
fm_FDR_flip[fm_FDR_flip < .1] = 2 
fm_FDR_flip[fm_FDR_flip > .1 & fm_FDR_flip <= 1] = 0.01  

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\dSV r-value heatmap.png', 
    height = 2000, width = 1600)

ggcorrplot(fm_r, 
           p.mat = fm_FDR_flip, 
           lab = TRUE, 
           outline.color = "white",
           ggtheme = ggplot2::theme_gray,
           sig.level = .1,
           pch = 1,
           pch.cex = 17) +
  scale_fill_gradient2(low = 'darkcyan', high = 'chocolate1', breaks=c(-1, 1), limit=c(-1, 1)) +
  theme(legend.key.height = unit(3, 'cm'),
        legend.key.width = unit(1.2, 'cm'),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 17),
        axis.text.x = element_text(size = 12, angle = 70)) +
  labs(fill = 'r value')

dev.off()


  #3.6 boxplot of best association

     #3.6.1 gather the relevant data

rel_data = dsv_corr_data[dsv_corr_data$dsv == "Blautia.wexlerae.DSM.19850.3737_3738", ]
rel_data = rel_data[rel_data$protein == "Ep.CAM", ]


rel_combo = paste('Blautia.wexlerae.DSM.19850.3737_3738', 'Ep.CAM', sep= ' - ')
fdr_val = rel_data$FDR
n_val = paste('n = ', rel_data$n, sep = '')
p_manual = data.frame('deletion', 'no.deletion', paste('FDR = ', format(fdr_val, digits = 3), sep = ''), 5.5)
colnames(p_manual) = c('group1','group2','p.adj', 'y.position')


pos_0 = which(d1[, 'Blautia.wexlerae.DSM.19850.3737_3738'] == 0)
pos_1 = which(d1[, 'Blautia.wexlerae.DSM.19850.3737_3738'] == 1)

protein_level_0 = d2[pos_0, 'Ep.CAM']
protein_level_1 = d2[pos_1, 'Ep.CAM']  

x = protein_level_0
y = protein_level_1

max_ln <- max(c(length(x), length(y)))
gfg_data<- data.frame('deletion' = c(x,rep(NA, max_ln - length(x))),
                      'no deletion' = c(y,rep(NA, max_ln - length(y))))

gfg_data = pivot_longer(gfg_data, cols = c('deletion','no.deletion'))


     #3.6.2 visualisation

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\dSV boxplot.png', 
    height = 700, width = 700)

ggboxplot(
  gfg_data, x = 'name', y = 'value', title = paste(rel_combo, n_val,sep = '\n'),
  ylab = 'protein level', xlab = 'groups') +
  stat_pvalue_manual(p_manual) +
  theme(axis.title = element_text(size = 16))


dev.off()
