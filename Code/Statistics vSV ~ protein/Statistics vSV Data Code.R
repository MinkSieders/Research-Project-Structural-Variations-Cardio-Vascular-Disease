
---
title: "Statistics vSV data"
authors: "Mink Sieders, Daan van der Giezen"
date: "2022/04/23"
---
  
  
#1 preparation

  #1.1 import functions

library(pacman)
p_load(ggplot2, corrplot, tidyr, Hmisc, psych, reshape2, dplyr, tidyverse, ggcorrplot, ggrepel, ggpubr)


  #1.2 import input files

all_vsv_df    =  read.table('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Data\\lld_vsv.tsv',
                           header = TRUE,
                           sep = '\t',
                           fill = TRUE)

all_prot_df   =  read.table('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Data\\lld_prot.tsv',
                            header = TRUE,
                            sep = '\t',
                            fill = TRUE)

vsv_info      =  read.delim('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Data\\vsgv_info_anno.tsv',
                           header = TRUE,
                           sep = '\t',
                           fill = TRUE)


#2 statistical analysis of the data

  #2.1 correlation analysis between vSV and protein levels

d1 = data.frame(all_vsv_df)
d2 = data.frame(all_prot_df)

vsv_corr_data = data.frame()

for(i in 1:length(colnames(d2))){
  for(j in 1:length(colnames(d1))){
    
    corr_test = psych::corr.test(d1[, j], d2[, i], method= 'spearman')
    corr_data = data.frame(colnames(d1)[j], colnames(d2)[i], corr_test$r, (corr_test$r)^2, corr_test$p)
    vsv_corr_data = rbind(vsv_corr_data, corr_data) } } 

p_adjust = p.adjust(vsv_corr_data[, 5], method = 'BH')
vsv_corr_data = cbind(vsv_corr_data, p_adjust)

colnames(vsv_corr_data) = c('vsv','protein','r', 'R2', 'p', 'FDR')


  #2.2 save statistical analysis file containing the statistical data values for later usage

save(vsv_corr_data, file = 'C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Results\\vsv_stat_anal_data.RData')


#3 visualisation

  #3.1 get the file with the statistical analysis data performed in '2'

vsv_corr_data = get(load('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Results\\vsv_stat_anal_data.RData'))

  
  #3.2 gathering the corresponding species names and making labels for relevant names (top 3 negative/positive corrolations)

x = data.frame(vsv_info$Taxonomy_Name)
spec_name = data.frame()
for(i in 1:92){
  spec_name = rbind(spec_name, x)}

n_corr_data = data.frame(cbind(vsv_corr_data$protein, vsv_corr_data$vsv, spec_name, vsv_corr_data$r, vsv_corr_data$FDR, vsv_corr_data$p))

colnames(n_corr_data) = c('protein', 'vsv', 'spec_name', 'r_value', 'FDR', 'p')

n_corr_data = transform(n_corr_data, r_value = as.numeric(r_value))
n_corr_data = transform(n_corr_data, FDR = as.numeric(FDR))
n_corr_data = transform(n_corr_data, p = as.numeric(p))

n_corr_data = n_corr_data %>%
  mutate(correlation = case_when(r_value >= .1 & p <= 0.05 ~ 'positive',
                                 r_value <= -.1 & p <= 0.05 ~ 'negative',
                                 TRUE ~ 'not significant'))

n_corr_data = n_corr_data[order(n_corr_data[,6],decreasing=FALSE),]

n_corr_data['name'] = NA

pos_rel = c()
neg_count = 0
labs_rel = c()
for(i in 1:1000){
  if(neg_count == 3){
    break
  }
  if(n_corr_data$correlation[i] == 'negative'){
     neg_count =  neg_count + 1  
     pos_rel = c(pos_rel, i)
     labs_rel = c(labs_rel, paste(n_corr_data$protein[i], n_corr_data$spec_name[i], sep = ' - ')) } }

n_corr_data[pos_rel, 8] = labs_rel

labs_rel2 = c()
pos_rel2 = c()
pos_count = 0
for(i in 1:1000){
  if(pos_count == 3){
    break
  }
  if(n_corr_data$correlation[i] == 'positive'){
    pos_count = pos_count + 1
    pos_rel2 = c(pos_rel, i)
    labs_rel2 = c(labs_rel2, paste(n_corr_data$protein[i], n_corr_data$spec_name[i], sep = ' - '))
  }
}

n_corr_data[pos_rel2, 8] = labs_rel2


  #3.3 volcano map of p values - visualisation of the correlation analysis between vSV and protein levels

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\vSV p-value volcano map.png', 
    height = 2000, width = 1500) 

line_pos = c(-0.2330253, -.1, .1, 0.2844537)

ggplot(n_corr_data, 
       aes(x = r_value, 
           y = -log10(p), 
           fill = correlation,    
           size = correlation,
           alpha = correlation)) + 
       geom_point(shape = 21, colour = 'black') + 
       geom_vline(xintercept = line_pos, 
                  linetype = 'dotdash') +
       geom_hline(yintercept = -log10(0.05), 
                  linetype = 'dashed') +
       geom_hline(yintercept = -log10(1e-5), 
                  linetype = 'dotted') +
       geom_label_repel(show.legend = FALSE, aes(label = name, fontface = 'bold'), size = 7,
                        box.padding = unit(2, 'lines')) +
       scale_x_continuous(breaks = line_pos, labels = c('-0.23', '-0.1', '0.1', '0.28')) +
       scale_y_continuous(breaks = c(-log10(1e-5), -log10(0.05)), labels = c('p < 1e-5', 'p < 0.05')) +
       scale_fill_manual(values = c('positive' = 'chocolate1', 'negative' = 'darkcyan', 'not significant' = 'grey')) + 
       scale_size_manual(values = c('positive' = 4, 'negative' = 4, 'not significant' = 2)) + 
       scale_alpha_manual(values = c('positive' = .8, 'negative' = .8, 'not significant' = .4)) +
       labs(x = 'r-value',
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


  #3.3 heatmap of the r-value - visualisation of the correlation analysis between vSV and protein levels

t200_ncd = n_corr_data

matrix_r = pivot_wider(t200_ncd, id_cols = vsv, names_from = protein, values_from = r_value)
matrix_r = column_to_rownames(matrix_r, 'vsv')

matrix_FDR = pivot_wider(t200_ncd, id_cols = vsv, names_from = protein, values_from = FDR)
matrix_FDR = column_to_rownames(matrix_FDR, 'vsv')

matrix_p = pivot_wider(t200_ncd, id_cols = vsv, names_from = protein, values_from = p)
matrix_p = column_to_rownames(matrix_p, 'vsv')

col_pos = c()
for(i in 1:length(colnames(matrix_p))){
  if(length(which(matrix_p[, i] < 0.0001))){
    col_pos = c(col_pos, i)} }

row_pos = c()
for(i in 1:length(rownames(matrix_p))){
  if(length(which(matrix_p[i, ] < 0.00005))){
    row_pos = c(row_pos, i)} }

fm_FDR = as.matrix(matrix_FDR[row_pos, col_pos])
fm_r = as.matrix(matrix_r[row_pos, col_pos])  

fm_FDR_flip = fm_FDR
fm_FDR_flip[fm_FDR_flip < .1] = 2 
fm_FDR_flip[fm_FDR_flip > .1 & fm_FDR_flip < 1] = 0.01  

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\vSV r-value heatmap.png', 
    height = 2000, width = 1600)

ggcorrplot(fm_r, 
           p.mat = fm_FDR_flip, 
           lab = TRUE, 
           outline.color = "white",
           ggtheme = ggplot2::theme_gray,
           sig.level = .1,
           pch = 1,
           pch.cex = 17) +
   scale_fill_gradient2(low = 'darkcyan', high = 'chocolate1', breaks=c(-.3, .3), limit=c(-.3, .3)) +
   theme(legend.key.height = unit(3, 'cm'),
         legend.key.width = unit(1.2, 'cm'),
         legend.title = element_text(size = 20),
         legend.text = element_text(size = 17),
         axis.text.x = element_text(size = 12, angle = 70)) +
   labs(fill = 'r value')
          
dev.off()


  #3.4 graph of best association

     #3.4.1 get the relevant data

data_vsv = all_vsv_df$Blautia.wexlerae.DSM.19850.3394_3400.3400_3404
data_prot = all_prot_df$Ep.CAM

my_data = data.frame(data_prot, data_vsv)

rel_data = vsv_corr_data[vsv_corr_data$vsv == 'Blautia.wexlerae.DSM.19850.3394_3400.3400_3404', ]
rel_data = rel_data[rel_data$protein == 'Ep.CAM', ]

fdr_val = format(rel_data$FDR, digits = 3)
n_val = length(na.omit(my_data[,2]))
r_val = format(rel_data$r, digits = 3)
rel_combo = paste('Blautia.wexlerae.DSM.19850.3394_3400.3400_3404', 'Ep.CAM', sep= ' - ')

stat_info = paste('R = ', r_val, ', n = ', n_val, ', FDR = ', fdr_val, sep = '')
     

     #3.4.2 visualize

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\vSV graph.png', 
    height = 700, width = 700)

ggscatter(my_data, x = 'data_vsv', y = 'data_prot', 
          add = 'reg.line', cor.method = 'spearman',conf.int = TRUE, 
          cor.coef = TRUE,
          xlab = "variable SV's", ylab = "protein levels", show.legend.text = FALSE, parse = FALSE, font.label = c(1, "bold", "white"),
          title = paste(rel_combo, stat_info,  sep = '\n')) +
          theme(axis.title = element_text(size = 16),
                legend.position="none")

dev.off()
