
---
title: "SV summary statistics"
authors: "Mink Sieders, Daan van der Giezen"
date: "2022/04/23"
---


#1 preparation
  
  #1.1 import functions
  
library(pacman)
p_load(ggplot2, tidyr)


  #1.2 input files

info_species   =  read.table('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Data\\Informative_species_information.tsv',
                             header = TRUE,
                             sep = '\t')


#2 SVs totals

  #2.1 calculate SV totals

dSV_sum = sum(info_species$Deletion_SVs_number)
vSV_sum = sum(info_species$Variable_SVs_number)

slices <- c(dSV_sum, vSV_sum)
lbs <- c("Deletion SV", "Variable SV")
pct_num <- paste(lbs, '\n', round(100*slices/sum(slices), 1), 
                 '% - ', slices, sep = '')

pie_df = data.frame(slices, lbs, pct_num)


  #2.2 visualize SV totals

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\ratio_total_SV.png', 
    height = 900, width = 900)

ggplot(pie_df, aes(x = '', y = slices, fill = lbs)) +
  geom_col(color = 'white', size = 2) +
  scale_fill_manual(values = c('lightblue3', 'dark cyan')) +
  coord_polar(theta = 'y') + 
  theme_void() + 
  theme(legend.position="none") +
  geom_text(aes(label = pct_num), color = "white", size=12, 
            position = position_stack(vjust = .5)) 

dev.off()


#3 SVs totals per species 

  #3.1 calculate SV totals per species

dSV_num = info_species$Deletion_SVs_number
vSV_num = info_species$Variable_SVs_number
spec_name = info_species$Short_name

info_df = data.frame(dSV_num, vSV_num, spec_name)
info_df_reshaped = tidyr::pivot_longer(info_df, cols = c(dSV_num, vSV_num), names_to = 'variable', 
                                      values_to = 'value')


  #3.2 visualize SV totals per species

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\Total SV distribution for each bacterial spiecies.png', 
    height = 900, width = 1300)

ggplot(info_df_reshaped, aes(x = reorder(spec_name, value), y = value, fill = variable)) + 
  geom_bar(width = 0.5, position='stack', stat = 'identity') +
  scale_fill_manual(name = "", labels = c("Deletion SV", "Variable SV"), values = c('lightblue3', 'darkcyan')) +
  ylab('Number of SVs') +
  theme(axis.text.x = element_text(size = 20, angle = 60, hjust = .95, vjust = 1), 
        axis.title = element_text(size = 35),
        legend.text = element_text(size = 30), 
        legend.title = element_text(size = 35, face = 'bold')) +
  labs(x = '', fill = '')      
  
dev.off()


#4 number of samples for each species

  #4.1 calculate number of samples for each species

spec_name = info_species$Short_name
sample_num = info_species$Total_samples_number

size_df = data.frame(spec_name, sample_num)


  #4.1 visualize number of samples for each species

png('C:\\Users\\minks\\OneDrive\\Bureaublad\\UNI\\Bachelor\\y4\\Semester II\\Research Project\\Figures\\sample size for each bacterial spiecies.png', 
    height = 900, width = 900)

ggplot(size_df, aes(x = reorder(spec_name, -sample_num), y = sample_num)) + 
       geom_bar(stat = "identity", color = "dark cyan", fill = "dark cyan", width = 0.5) + 
       ylab('Number of samples') +
       theme(axis.text.y = element_text(size = 20), 
             axis.title = element_text(size = 35), 
             axis.title.y = element_blank()) +
       coord_flip()

dev.off()

