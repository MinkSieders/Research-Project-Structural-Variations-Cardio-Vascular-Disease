
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

all_dsv_df    =  read.table('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Data\\lld_dsv.tsv',
                           header = TRUE,
                           sep = '\t',
                           fill = TRUE)

all_cvd_df    =  read.table('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Data\\lld_cvd.tsv',
                            header = TRUE,
                            sep = '\t',
                            fill = TRUE)

dsv_corr_data = get(load('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Results\\dsv_stat_anal_data.RData'))


#2 statistical analysis of the data

  #2.1 correlation analysis between vSV and protein levels

d1 = data.frame(all_dsv_df)
d2 = data.frame(all_cvd_df)

d2 = d2 %>%
  mutate(Sex = case_when(Sex == 1 ~ 'male',
                                 Sex == 0 ~ 'female'))

analyzing = 'male' # choose either 'male' or 'female' to do analysis

filter_pos = c()
for(i in 1:length(rownames(d2))){
  if(d2$Sex[i] == analyzing){
    filter_pos = c(filter_pos, i)
  }
}

Sex_d2 = d2[filter_pos,]
Sex_d1 = d1[filter_pos,]

dsv_corr_data_cvd = data.frame()

  for(j in 1:length(colnames(Sex_d1))){
    
    if(length(which(all_dsv_df[, j] == 0)) >= 10 & length(which(all_dsv_df[, j] == 1)) >= 10){
      pos_0 = which(Sex_d1[, j] == 0)
      pos_1 = which(Sex_d1[, j] == 1)
      
      protein_level_0 = Sex_d2[pos_0, 1]
      protein_level_1 = Sex_d2[pos_1, 1]  
      
      corr_test = wilcox.test(protein_level_0, protein_level_1, paired = FALSE)
      corr_data = data.frame(colnames(Sex_d1)[j], corr_test$p.value, (mean(protein_level_0, na.rm = TRUE) - 
                                                                            mean(protein_level_1, na.rm = TRUE)))
      
      dsv_corr_data_cvd = rbind(dsv_corr_data_cvd, corr_data)
      
    }else{
      corr_data = data.frame(colnames(Sex_d1)[j], NA, NA)
      colnames(corr_data) = colnames(dsv_corr_data_cvd)
      
      dsv_corr_data_cvd = rbind(dsv_corr_data_cvd, corr_data)}
  } 


p_adjust = p.adjust(dsv_corr_data_cvd[, 2], method = 'BH')
dsv_corr_data_cvd = cbind(dsv_corr_data_cvd, p_adjust)

colnames(dsv_corr_data_cvd) = c('dsv', 'p', 'mean_difference', 'FDR')

dsv_corr_data_cvd = dsv_corr_data_cvd %>%
  mutate(significant = case_when(FDR <= .01 ~ '  ***',
                                 FDR <= .05 ~ '  **',
                                 FDR <= .1 ~ '  *',
                                 TRUE ~ ' '))


#3 visualisation

  #3.1 gathering the corresponding species names and making labels for relevant names (top 3 negative/positive corrolations)

n_corr_data = dsv_corr_data_cvd 

n_corr_data = n_corr_data %>%
  mutate(correlation = case_when(mean_difference >= .1 & p <= 0.05 ~ 'positive',
                                 mean_difference <= -.1 & p <= 0.05 ~ 'negative',
                                 TRUE ~ 'not significant'))

n_corr_data = n_corr_data[order(n_corr_data[,2],decreasing=FALSE),]

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
     labs_rel = c(labs_rel, n_corr_data$dsv[i]) } }

n_corr_data[pos_rel, 7] = labs_rel

labs_rel2 = c()
pos_rel2 = c()
pos_count = 0
for(i in 1:1000){
  if(pos_count == 3){
    break
  }
  if(n_corr_data$correlation[i] == 'positive'){
    pos_count = pos_count + 1
    pos_rel2 = c(pos_rel2, i)
    labs_rel2 = c(labs_rel2, n_corr_data$dsv[i])
  }
}

n_corr_data[pos_rel2, 7] = labs_rel2


fn_corr_data = n_corr_data #use this one when you are analyzing 'female'
mn_corr_data = n_corr_data #use this one when you are analyzing 'male'
fn_corr_data['gender'] = 'female'
mn_corr_data['gender'] = 'male'

n_combined_corr_data = rbind(fn_corr_data, mn_corr_data) #combine the two for the barplots


  #3.2 volcano map of p values - visualisation of the correlation analysis between vSV and disease risk

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\dSV cvd p-value volcano map.png', 
    height = 2000, width = 1500) 

line_pos = c(-3.166149, -.1, .1, 3.433475)

ggplot(n_corr_data, 
       aes(x = mean_difference, 
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
       scale_x_continuous(breaks = line_pos, labels = c('-3.1', '-.1', '.1', '3.4'), limits = c(-4, 4)) +
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


  #3.3 barplot of r2 - visualisation of the correlation analysis between vsvs and disease riks using vsvs that are associated with proteins (FDR cutoff at .1)

dsv_corr_data = dsv_corr_data[order(dsv_corr_data[,8],decreasing=FALSE),]

associatednames = c()
for(i in 1:1000){
  if(dsv_corr_data[i, 8] >= 0.1){
    break   }
  uno = dsv_corr_data[i, 1]
  associatednames = c(associatednames, uno)
  
}

bp_corr_data = data.frame()
for(i in 1:length(associatednames)){
  uno = n_corr_data[n_corr_data$dsv == associatednames[i], ]
  bp_corr_data = rbind(bp_corr_data, uno)
}

bp_corr_data = bp_corr_data %>%
mutate(correlation2 = case_when(mean_difference >= .1 & FDR <= 0.1 ~ 'positive',
                                mean_difference <= -.1 & FDR <= 0.1 ~ 'negative',
                                TRUE ~ 'not significant'))

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\dsv cvd barplot1.png', 
    height = 700, width =900) 

ggplot(bp_corr_data, aes(y = abs(mean_difference), x = reorder(dsv, -(abs(mean_difference))), fill = correlation2)) + 
  geom_bar(position = "dodge", stat = 'identity') +
  scale_fill_manual(values = c('grey')) +
  geom_text(aes(label = significant, angle = 90), position = position_dodge(0.9), vjust = .5, hjust = 0.04, size = 6) +
  #scale_y_continuous( limits = c(0, 0.17)) +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 17),
        axis.text.x = element_text(size = 12, angle = 75, hjust = .95, vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20)) +
  labs(fill = 'protein',
       y = 'mean difference in disease risk')

dev.off()


fn_corr_data = (fn_corr_data[order(fn_corr_data[,1],decreasing=FALSE),])
mn_corr_data = (mn_corr_data[order(mn_corr_data[,1],decreasing=FALSE),])

fn_corr_data = fn_corr_data[!is.na(fn_corr_data$FDR),]
mn_corr_data = mn_corr_data[!is.na(mn_corr_data$FDR),]

indexes = c()
f_count = 0
m_count = 0
for(i in 1:length(row.names(fn_corr_data))){
  if(fn_corr_data$FDR[i] <= .1){
    indexes = c(indexes, i)
    f_count = f_count + 1
  }
  if(mn_corr_data$FDR[i] <= .1){
    indexes = c(indexes, i)
    m_count = m_count + 1
  }
  
}

bp2_corr_data = distinct(rbind(fn_corr_data[indexes,], mn_corr_data[indexes,]))


png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\dsv cvd barplot top20.png', 
    height = 700, width = 1300) 

ggplot(bp2_corr_data, aes(y = mean_difference, x = reorder(dsv, -mean_difference), fill = gender)) + 
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = c('darkcyan', 'chocolate1')) +
  #geom_text(aes(label = significant, angle = 90), position = position_dodge(0.9), vjust = .5, hjust = 0.04, size = 6) +
  #scale_y_continuous( limits = c(0, 7.5)) +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 17),
        axis.text.x = element_text(size = 15, angle = 75, hjust = .95, vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15)) +
  labs(fill = 'gender',
       y = 'mean difference in disease risk')

dev.off()
