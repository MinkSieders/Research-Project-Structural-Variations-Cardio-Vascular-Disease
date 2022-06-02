
---
title: "Statistics vSV data"
authors: "Mink Sieders, Daan van der Giezen"
date: "2022/04/23"
---
  
  
#1 preparation

  #1.1 import functions

library(pacman)
p_load(ggplot2, corrplot, tidyr, Hmisc, psych, reshape2, dplyr, tidyverse, ggcorrplot, ggrepel, ggpubr, scales)


  #1.2 import input files

all_vsv_df    =  read.table('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Data\\lld_vsv.tsv',
                           header = TRUE,
                           sep = '\t',
                           fill = TRUE)

vsv_info      =  read.delim('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Data\\vsgv_info_anno.tsv',
                           header = TRUE,
                           sep = '\t',
                           fill = TRUE)

all_cvd_df    =  read.table('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Data\\lld_cvd.tsv',
                            header = TRUE,
                            sep = '\t',
                            fill = TRUE)

vsv_corr_data = get(load('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Results\\vsv_stat_anal_data.RData'))


#2 statistical analysis of the data

  #2.1 correlation analysis between vSV and protein levels

d1 = data.frame(all_vsv_df)
d2 = data.frame(all_cvd_df)

d2 = d2 %>%
  mutate(Sex = case_when(Sex == 1 ~ 'male',
                                 Sex == 0 ~ 'female'))

analyzing = 'female' # choose either 'male' or 'female' to do analysis

filter_pos = c()
for(i in 1:length(rownames(d2))){
  if(d2$Sex[i] == analyzing){
    filter_pos = c(filter_pos, i)
  }
}

Sex_d2 = d2[filter_pos,]
Sex_d1 = d1[filter_pos,]

vsv_corr_data_cvd = data.frame()

  for(j in 1:length(colnames(Sex_d1))){
    corr_test = psych::corr.test(Sex_d1[, j], Sex_d2[1], method= 'spearman')
    corr_data = data.frame(colnames(d1)[j], corr_test$r, (corr_test$r)^2, corr_test$p)
    vsv_corr_data_cvd = rbind(vsv_corr_data_cvd, corr_data) }

p_adjust = p.adjust(vsv_corr_data_cvd[, 4], method = 'BH')
vsv_corr_data_cvd = cbind(vsv_corr_data_cvd, p_adjust)

colnames(vsv_corr_data_cvd) = c('vsv','r', 'R2', 'p', 'FDR')

vsv_corr_data_cvd = vsv_corr_data_cvd %>%
  mutate(significant = case_when(FDR <= .01 ~ '  ***',
                                 FDR <= .05 ~ '  **',
                                 FDR <= .1 ~ '  *',
                                 TRUE ~ ' '))


#3 visualisation

  #3.1 gathering the corresponding species names and making labels for relevant names (top 3 negative/positive corrolations)

n_corr_data = vsv_corr_data_cvd
n_corr_data = n_corr_data %>%
  mutate(correlation = case_when(r >= .1 & p <= 0.05 ~ 'positive',
                                 r <= -.1 & p <= 0.05 ~ 'negative',
                                 TRUE ~ 'not significant'))

n_corr_data = n_corr_data[order(n_corr_data[,5],decreasing=FALSE),]

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
     labs_rel = c(labs_rel, n_corr_data$vsv[i]) } }

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
    pos_rel2 = c(pos_rel2, i)
    labs_rel2 = c(labs_rel2, n_corr_data$vsv[i])
  }
}

n_corr_data[pos_rel2, 8] = labs_rel2

fn_corr_data = n_corr_data #use this one when you are analyzing 'female'
mn_corr_data = n_corr_data #use this one when you are analyzing 'male'
fn_corr_data['gender'] = 'female'
mn_corr_data['gender'] = 'male'

n_combined_corr_data = rbind(fn_corr_data, mn_corr_data) #combine the two for the barplots

  #3.2 volcano map of p values - visualisation of the correlation analysis between vSV and disease risk

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\vSV cvd p-value volcano map.png', 
    height = 2000, width = 1500) 

line_pos = c(-0.2853201, -.1, .1, 0.29443687)

ggplot(n_corr_data, #use either the fn_corr_data or the mn_corr_data depending on what gender you want to visualize
       aes(x = r, 
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
                        box.padding = unit(5, 'lines')) +
       scale_x_continuous(breaks = line_pos, labels = c('-0.28', '-0.1', '0.1', '0.29')) +
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


  #3.3 barplot of r2 - visualisation of the correlation analysis between vsvs and disease riks using vsvs that are associated with proteins (FDR cutoff at .1)

vsv_corr_data = vsv_corr_data[order(vsv_corr_data[,6],decreasing=FALSE),]

associatednames = c()
for(i in 1:1000){
  if(vsv_corr_data[i, 6] >= 0.1){
    break   }
  uno = vsv_corr_data[i, 1]
  associatednames = c(associatednames, uno)

}

bp_corr_data = data.frame()
for(i in 1:length(associatednames)){
  uno = n_combined_corr_data[n_combined_corr_data$vsv == associatednames[i], ]
  bp_corr_data = rbind(bp_corr_data, uno)
}
 

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\cvd barplot1.png', 
    height = 700, width = 850) 

ggplot(bp_corr_data, aes(y = r, x = reorder(vsv, -r), fill = gender)) + 
  geom_bar(position = "dodge", stat = 'identity') +
  scale_fill_manual(values = c('darkcyan', 'chocolate1')) +
  geom_text(aes(label = significant, angle = 90), position = position_dodge(0.9), vjust = .5, hjust = 0.04, size = 6) +
  #scale_y_continuous( limits = c(0, 0.17)) +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 17),
        axis.text.x = element_text(size = 12, angle = 70, hjust = .95, vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20)) +
  labs(fill = 'gender',
       y = 'variance explained by vSV')

dev.off()


    #get the top 10 associations from both females and males and put them in a dataframe. 
      #Also finds the corresponding female or male associations for each sex so we can compare them

fn_corr_data = fn_corr_data[order(fn_corr_data[,1],decreasing=FALSE),]
mn_corr_data = mn_corr_data[order(mn_corr_data[,1],decreasing=FALSE),]

indexes = c()
f_count = 0
m_count = 0
for(i in 1:length(row.names(fn_corr_data))){
  if(fn_corr_data$FDR[i] <= .05){
    indexes = c(indexes, i)
    f_count = f_count + 1
  }
  if(mn_corr_data$FDR[i] <= .05){
    indexes = c(indexes, i)
    m_count = m_count + 1
  }
  
}

bp2_corr_data = distinct(rbind(fn_corr_data[indexes,], mn_corr_data[indexes,]))

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\cvd barplot top.png', 
    height = 700, width = 1300) 

ggplot(bp2_corr_data, aes(y = r, x = reorder(vsv, -r), fill = gender)) + 
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = c('darkcyan', 'chocolate1')) +
  geom_text(aes(label = significant, angle = 90), position = position_dodge(0.9), size = 6) +
  scale_y_continuous( limits = c(0, 0.17)) +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 17),
        axis.text.x = element_text(size = 10, angle = 70, hjust = .95, vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20)) +
  labs(fill = 'gender',
       y = 'r value')

dev.off()

