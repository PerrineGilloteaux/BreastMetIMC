####PART 2 ABUNDANCE

#Won Jin Ho 
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 18.04.5 LTS

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
work<-getwd()

####REQUIRED LIBRARIES####
library(readxl)
library(stringr)
library(matrixcalc)
library(matrixStats)
library(Hmisc)
library(reshape2)
library(dplyr)
library(plotrix)
library(multcomp)
library(flowCore)
library(sf)
library(clusterSim)
library(limma)
library(corrplot)
library(packcircles)
library(cowplot)
library(autoimage)
library(ggplot2)
library(ggpubr)
library(ggiraphExtra)
library(ggvoronoi)
library(ggridges)
library(gridExtra)
library(igraph)
library(qgraph)
library(circlize)
library(scales)
library(RColorBrewer)
library(ComplexHeatmap)
library(pheatmap)
library(pals)
library(plot3D)
library(akima)
library(basetheme)

####LOAD CLUSTERED DATA####
output<-readRDS('..\\Input_Files\\backup_output.rds')
clusterMergeFile = "..\\Input_Files\\merge.xlsx" #create dummy merger numbers prior to annotation
cluster_merging <- read_excel(clusterMergeFile)

##Levels

samplevels=c(paste0("P_",1:342))
samplestoexclude <- c(output$meta_data$sample_id[output$meta_data$Tumor=="CTRL"]) # CTRL 
typelevels=c("CTRL","LUM","TNC","HER2")
sitelevels=c("PBC","PBC_1YR","PBC_tx","PBC_recur","GI","PANC","LIVER","BRAIN","OVARY","PLEURA","SPINE","LUNG","LN","CTRL")
tumorlevels=c("IDC","ILC","IMC","CTRL")
TMAlevels=c("TMA788","TMA789","TMA801","TMA974","TMA975")
SPClevels=c("SPC01","SPC02","SPC03","SPC04","SPC06","SPC07","SPC09","SPC10",
            "SPC11","SPC12","SPC13","SPC14","SPC15","SPC16","SPC17","SPC18",
            "SPC19","SPC20","SPC21","SPC22","SPC23","SPC24","SPC25","SPC26",
            "CTRL")
clusterlevels=c("Tc", "ThN", "ThEM", "Treg",
                "NK_T", "NK", "B_I", "B_II",
                "Gran", "DC",
                "Mac_I", "Mac_II", "Mac_III", "Mac_IV",
                "Str_PS", "Str_V", "Str_VC", "Str_VS",
                "CK8", "CK8_E", "CK8_V",
                "CK14", "CK14_E", "CK14_EV",
                "CK8_14_E", "CK8_14_EV",
                "CKlo", "CKlo_EV", "CKother",
                "Neuron", "NA")
clusterlevels_macs=c("CK8", "CK8_E", "CK8_V",
                     "CK14", "CK14_E", "CK14_EV",
                     "CK8_14_E", "CK8_14_EV",
                     "CKlo", "CKlo_EV", "CKother", 
                     "Neuron", "NA",
                     "Str_PS", "Str_V", "Str_VC", "Str_VS",
                     "NK_T", "NK", "B_I", "B_II",                     
                     "Tc", "ThN", "ThEM", "Treg", 
                     "Gran", "DC",
                     "Mac_I", "Mac_II", "Mac_III", "Mac_IV")

##colors

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
hex <- hue_pal()(9)
colorassignedbroad <- c(rep(hex[1],3), #ck8
                        rep(hex[2],3), #ck14
                        rep(hex[3],2), #ck8and14
                        rep(hex[4],3), #ckother
                        rep(hex[5],2), #other
                        rep(hex[6],4), #str
                        rep(hex[7],8), #lym
                        rep(hex[8],2), #myeloid1
                        rep(hex[9],4)) #myeloid2
clusternames<-clusterlevels
names(colorassigned)<-clusternames


####ABUNDANCE PLOTS####

## Proportion calculations
counts_table <- table(output$cell_clustering1m, output$sample_ids)
counts_table_noNA <- counts_table[clusterlevels[clusterlevels != "NA"],]
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
props_table_noNA <- t(t(counts_table_noNA) / colSums(counts_table_noNA)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)
props_noNA <- as.data.frame.matrix(props_table_noNA)

#areas <- read_xlsx(paste0(work,'/Config/areas.xlsx'))
#densities <- t(t(counts)/areas$TotalArea)

write.csv(counts,'Results_counts.csv')
write.csv(props,'Results_props.csv')
write.csv(props_noNA, "Results_props_noNA.csv")
#write.csv(densities, 'Results_densities.csv')

## Set up the data frame for proportional plotting
ggdf <- melt(data.frame(cluster = rownames(props_noNA),props_noNA, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdf$SPC <- factor(output$meta_data$Case[match(ggdf$sample_id,output$meta_data$sample_id)], levels = SPClevels)
ggdf$tumor <- factor(output$meta_data$Tumor[match(ggdf$sample_id,output$meta_data$sample_id)], levels=tumorlevels)
ggdf$type <- factor(output$meta_data$Phenotype[match(ggdf$sample_id,output$meta_data$sample_id)], levels=typelevels)
ggdf$site <- factor(output$meta_data$Anno[match(ggdf$sample_id,output$meta_data$sample_id)], levels=sitelevels)
ggdf$TMA <- factor(output$meta_data$TMA[match(ggdf$sample_id,output$meta_data$sample_id)], levels=TMAlevels)
ggdf$met <- factor(output$meta_data$Met[match(ggdf$sample_id,output$meta_data$sample_id)])

ggdf<-ggdf[ggdf$sample_id %nin% samplestoexclude,]

#plot box plots

#% CELLS
ggp2<-ggplot(ggdf,aes(x=SPC,y=proportion,fill=met))+
  geom_boxplot(outlier.shape=NA, lwd=0.25, color="black")+
  geom_jitter(width=0.2, size=0.5)+
  facet_wrap(~cluster,ncol=4,scales="free")+
  ylab("% of Cells")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(linewidth=0.25, color="black"),
        axis.line.y = element_line(linewidth=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=7, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white")  
  )

pdf('Abundance_box.pdf',width=16,height=8)
ggp2
dev.off()

## Stacked bars

#to get out the sample levels based on ranking the % of total macrophages 
ggdf_noNA <- ggdf[ggdf$cluster %nin% c("NA"),]
ggdf_noNA_macs <- ggdf_noNA[ggdf_noNA$cluster == "Mac_I",]
ggdf_noNA_macs1 <- ggdf_noNA[ggdf_noNA$cluster == "Mac_I",]
ggdf_noNA_macs1 <- ggdf_noNA_macs1$proportion
ggdf_noNA_macs2 <- ggdf_noNA[ggdf_noNA$cluster == "Mac_II",]
ggdf_noNA_macs2 <- ggdf_noNA_macs2$proportion
ggdf_noNA_macs3 <- ggdf_noNA[ggdf_noNA$cluster == "Mac_III",]
ggdf_noNA_macs3 <- ggdf_noNA_macs3$proportion
ggdf_noNA_macs4 <- ggdf_noNA[ggdf_noNA$cluster == "Mac_IV",]
ggdf_noNA_macs4 <- ggdf_noNA_macs4$proportion
ggdf_noNA_macstot<-ggdf_noNA_macs1+ggdf_noNA_macs2+ggdf_noNA_macs3+ggdf_noNA_macs4
ggdf_noNA_macstot<-data.frame(sample_id=ggdf_noNA_macs$sample_id, tot=ggdf_noNA_macstot)
ggdf_noNA_macstot <- arrange(ggdf_noNA_macstot,desc(tot))
macptlevels<-ggdf_noNA_macstot$sample_id

ggdf_noNA$sample_id <- factor(ggdf_noNA$sample_id, levels = macptlevels)
ggdf_noNA$cluster <- factor(ggdf_noNA$cluster, levels = clusterlevels_macs)
ggdf_noNA$site2 <- ggdf_noNA$site
ggdf_noNA$site2[str_detect(ggdf_noNA$site2,"PBC")] <- "PBC"

bp <- ggplot(ggdf_noNA, aes(x = sample_id, y = proportion, fill=cluster)) +
  geom_bar(stat = "identity", position="stack", width=0.85) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
  ) +
  ylab("% of Total Cells")+
  facet_grid(.~site2, space="free", scales="free_x")+
  #scale_fill_manual(values = colorassigned, breaks = clusterlevels, labels = clusterlevels)+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=2))

bp2 <- ggplot(ggdf_noNA, aes(x = sample_id, y = proportion, fill=cluster)) +
  geom_bar(stat = "identity", position="stack", width=0.85) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
  ) +
  ylab("% of Total Cells")+
  facet_grid(.~type, space="free", scales="free_x")+
  #scale_fill_manual(values = colorassigned, breaks = clusterlevels, labels = clusterlevels)+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=2))

pdf('Abundance_stackedbar_noNA.pdf', width=20, height=4); bp; bp2; dev.off()

## Compact stacked bar plots

bp <- ggplot(ggdf_noNA, aes(x = sample_id, y = proportion, fill=cluster)) +
  geom_col(position="stack", width=0.85) +
  ggtitle('All Samples')+
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
  ) +
  ylab("% of Total Cells")+
  facet_grid(.~site2, space="free", scales="free_x")+
  scale_fill_manual(values = colorassignedbroad, breaks = clusterlevels_macs, labels = clusterlevels_macs)+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=2))

bp2 <- ggplot(ggdf_noNA, aes(x = sample_id, y = proportion, fill=cluster)) +
  geom_col(position="stack", width=0.85) +
  ggtitle('All Samples')+
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.ticks.x = element_blank(),        
        axis.ticks.y = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
  ) +
  ylab("% of Total Cells")+
  facet_grid(.~type, space="free", scales="free_x")+
  scale_fill_manual(values = colorassignedbroad, breaks = clusterlevels_macs, labels = clusterlevels_macs)+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=2))

bp3 <- ggplot(ggdf_noNA[ggdf_noNA$type=="TNC",], aes(x = sample_id, y = proportion, fill=cluster)) +
  geom_col(position="stack", width=0.85) +
  ggtitle('TNC') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.ticks.x = element_blank(),        
        axis.ticks.y = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
  ) +
  ylab("% of Total Cells")+
  facet_grid(.~site2, space = "free", scales = "free_x")+
  scale_fill_manual(values = colorassignedbroad, breaks = clusterlevels_macs, labels = clusterlevels_macs)+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=2))

bp4 <- ggplot(ggdf_noNA[ggdf_noNA$type=="LUM",], aes(x = sample_id, y = proportion, fill=cluster)) +
  geom_col(position="stack", width=0.85) +
  ggtitle('LUM') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.ticks.x = element_blank(),        
        axis.ticks.y = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white")
  ) +
  ylab("% of Total Cells")+
  facet_grid(.~site2, space = "free", scales = "free_x")+
  scale_fill_manual(values = colorassignedbroad, breaks = clusterlevels_macs, labels = clusterlevels_macs)+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=2))

pdf('Abundance_stackedbar_compact_noNA.pdf', width=9, height=3); bp; bp2; dev.off()

pdf('Abundance_stackedbar_compact2_noNA.pdf', width=9, height=6); 
ggdraw()+
  draw_plot(bp3, x=0, y=0.5, width = 0.47, height = 0.5)+
  draw_plot(bp4, x=0, y=0, width = 1, height = 0.5)
dev.off()

## Stacked bar plots by patients

ggdf_byspc <- ggdf_noNA %>% group_by(cluster, SPC, met, type, site2) %>% summarise(avg = mean(proportion))

ggdf_noNA_macstot <- arrange(ggdf_noNA_macstot,desc(tot))
ggdf_noNA_macstot$SPC <- factor(output$meta_data$Case[match(ggdf_noNA_macstot$sample_id,output$meta_data$sample_id)])
ggdf_noNA_macstot$site <- factor(output$meta_data$Anno[match(ggdf_noNA_macstot$sample_id,output$meta_data$sample_id)])
ggdf_noNA_macstot$site2 <- ggdf_noNA_macstot$site
ggdf_noNA_macstot$site2[str_detect(ggdf_noNA_macstot$site2,"PBC")] <- "PBC"

ggdf_noNA_macstot_byspc <- ggdf_noNA_macstot %>% group_by(SPC, site2) %>% summarise(avg = mean(tot))
ggdf_noNA_macstot_byspc <- arrange(ggdf_noNA_macstot_byspc[ggdf_noNA_macstot_byspc$site2=="PBC",],desc(avg))
macagglevels<-ggdf_noNA_macstot_byspc$SPC

ggdf_byspc$SPC <- factor(ggdf_byspc$SPC, levels = macagglevels)

bp3 <- ggplot(ggdf_byspc[ggdf_byspc$type=="TNC",], aes(x = SPC, y = avg, fill=cluster)) +
  geom_col(position="stack", width=0.75) +
  ggtitle('TNC') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.ticks.x = element_blank(),        
        axis.ticks.y = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
  ) +
  ylab("% of Total Cells")+
  facet_grid(.~site2, space = "free", scales = "free_x")+
  scale_fill_manual(values = colorassignedbroad, breaks = clusterlevels_macs, labels = clusterlevels_macs)+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=2))

bp4 <- ggplot(ggdf_byspc[ggdf_byspc$type=="LUM",], aes(x = SPC, y = avg, fill=cluster)) +
  geom_col(position="stack", width=0.75) +
  ggtitle('LUM') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.ticks.x = element_blank(),        
        axis.ticks.y = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white")
  ) +
  ylab("% of Total Cells")+
  facet_grid(.~site2, space = "free", scales = "free_x")+
  scale_fill_manual(values = colorassignedbroad, breaks = clusterlevels_macs, labels = clusterlevels_macs)+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=2))

pdf('Abundance_stackedbar_compact_noNA_bySPC.pdf', width=9, height=6); 
ggdraw()+
  draw_plot(bp3, x=0, y=0.5, width = 0.25, height = 0.5)+
  draw_plot(bp4, x=0, y=0, width = .49, height = 0.5)
dev.off()

#Stacked bar plots only for macs

ggdf_noNA_macs_sub <- ggdf_noNA[ggdf_noNA$cluster %in% c("Mac_I","Mac_II","Mac_III","Mac_IV"),]
ggdf_noNA_macs_sub$sample_id <- factor(ggdf_noNA_macs_sub$sample_id, levels = macptlevels)
ggdf_noNA_macs_sub$site2 <- ggdf_noNA_macs_sub$site
ggdf_noNA_macs_sub$site2[str_detect(ggdf_noNA_macs_sub$site2,"PBC")] <- "PBC"

bp5 <- ggplot(ggdf_noNA_macs_sub[ggdf_noNA_macs_sub$type=="LUM",], aes(x = sample_id, y = proportion, fill=cluster)) +
  geom_col(position="stack", width=0.85) +
  ggtitle('LUM Macs') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.ticks.x = element_blank(),        
        axis.ticks.y = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white")
  ) +
  ylab("% of Total Cells")+
  facet_grid(.~site2, space = "free", scales = "free_x")+
  scale_fill_manual(values = brewer.dark2(4))+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=1))

bp6 <- ggplot(ggdf_noNA_macs_sub[ggdf_noNA_macs_sub$type=="TNC",], aes(x = sample_id, y = proportion, fill=cluster)) +
  geom_col(position="stack", width=0.85) +
  ggtitle('TNC Macs') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.ticks.x = element_blank(),        
        axis.ticks.y = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white")
  ) +
  ylab("% of Total Cells")+
  facet_grid(.~site2, space = "free", scales = "free_x")+
  scale_fill_manual(values = brewer.dark2(4))+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=1))

pdf('Abundance_stackedbar_compact3_noNA.pdf', width=9, height=6); 
ggdraw()+
  draw_plot(bp6, x=0, y=0.5, width = 0.6, height = 0.5)+
  draw_plot(bp5, x=0, y=0, width = 1, height = 0.5)
dev.off()

#Stacked bar plots only for macs by patients

ggdf_noNA_macs_sub <- ggdf_noNA_macs_sub %>% group_by(cluster, SPC, met, type, site2) %>% summarise(avg = mean(proportion))
ggdf_noNA_macs_sub$SPC <- factor(ggdf_noNA_macs_sub$SPC, levels = macagglevels)

bp5 <- ggplot(ggdf_noNA_macs_sub[ggdf_noNA_macs_sub$type=="LUM",], aes(x = SPC, y = avg, fill=cluster)) +
  geom_col(position="stack", width=0.75) +
  ggtitle('LUM Macs') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.ticks.x = element_blank(),        
        axis.ticks.y = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white")
  ) +
  ylab("% of Total Cells")+
  facet_grid(.~site2, space = "free", scales = "free_x")+
  scale_fill_manual(values = brewer.dark2(4))+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=1))

bp6 <- ggplot(ggdf_noNA_macs_sub[ggdf_noNA_macs_sub$type=="TNC",], aes(x = SPC, y = avg, fill=cluster)) +
  geom_col(position="stack", width=0.75) +
  ggtitle('TNC Macs') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.ticks.x = element_blank(),        
        axis.ticks.y = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white")
  ) +
  ylab("% of Total Cells")+
  facet_grid(.~site2, space = "free", scales = "free_x")+
  scale_fill_manual(values = brewer.dark2(4))+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=1))

pdf('Abundance_stackedbar_compact3_noNA_bySPC.pdf', width=9, height=6); 
ggdraw()+
  draw_plot(bp6, x=0, y=0.5, width = 0.34, height = 0.5)+
  draw_plot(bp5, x=0, y=0, width = 0.49, height = 0.5)
dev.off()

#Stacked bar plots only for stromal cells

ggdf_noNA <- ggdf[ggdf$cluster %nin% c("NA"),]
ggdf_noNA_str <- ggdf_noNA[ggdf_noNA$cluster == "Str_PS",]
ggdf_noNA_str1 <- ggdf_noNA[ggdf_noNA$cluster == "Str_PS",]
ggdf_noNA_str1 <- ggdf_noNA_str1$proportion
ggdf_noNA_str2 <- ggdf_noNA[ggdf_noNA$cluster == "Str_V",]
ggdf_noNA_str2 <- ggdf_noNA_str2$proportion
ggdf_noNA_str3 <- ggdf_noNA[ggdf_noNA$cluster == "Str_VC",]
ggdf_noNA_str3 <- ggdf_noNA_str3$proportion
ggdf_noNA_str4 <- ggdf_noNA[ggdf_noNA$cluster == "Str_VS",]
ggdf_noNA_str4 <- ggdf_noNA_str4$proportion
ggdf_noNA_strtot<-ggdf_noNA_str1+ggdf_noNA_str2+ggdf_noNA_str3+ggdf_noNA_str4
ggdf_noNA_strtot<-data.frame(sample_id=ggdf_noNA_str$sample_id, tot=ggdf_noNA_strtot)
ggdf_noNA_strtot <- arrange(ggdf_noNA_strtot,desc(tot))
strptlevels<-ggdf_noNA_strtot$sample_id

ggdf_noNA_strtot$SPC <- factor(output$meta_data$Case[match(ggdf_noNA_strtot$sample_id,output$meta_data$sample_id)])
ggdf_noNA_strtot$site <- factor(output$meta_data$Anno[match(ggdf_noNA_strtot$sample_id,output$meta_data$sample_id)])
ggdf_noNA_strtot$site2 <- ggdf_noNA_strtot$site
ggdf_noNA_strtot$site2[str_detect(ggdf_noNA_strtot$site2,"PBC")] <- "PBC"

ggdf_noNA_str_sub <- ggdf_noNA[ggdf_noNA$cluster %in% c("Str_PS","Str_V","Str_VC","Str_VS"),]
ggdf_noNA_str_sub$sample_id <- factor(ggdf_noNA_str_sub$sample_id, levels = strptlevels)
ggdf_noNA_str_sub$site2 <- ggdf_noNA_str_sub$site
ggdf_noNA_str_sub$site2[str_detect(ggdf_noNA_str_sub$site2,"PBC")] <- "PBC"

ggdf_noNA_strtot_byspc <- ggdf_noNA_strtot %>% group_by(SPC, site2) %>% summarise(avg = mean(tot))
ggdf_noNA_strtot_byspc <- arrange(ggdf_noNA_strtot_byspc[ggdf_noNA_strtot_byspc$site2=="PBC",],desc(avg))
stragglevels<-ggdf_noNA_strtot_byspc$SPC

bp5 <- ggplot(ggdf_noNA_str_sub[ggdf_noNA_str_sub$type=="LUM",], aes(x = sample_id, y = proportion, fill=cluster)) +
  geom_col(position="stack", width=0.85) +
  ggtitle('LUM Str') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.ticks.x = element_blank(),        
        axis.ticks.y = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white")
  ) +
  ylab("% of Total Cells")+
  facet_grid(.~site2, space = "free", scales = "free_x")+
  scale_fill_manual(values = brewer.set2(4))+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=1))

bp6 <- ggplot(ggdf_noNA_str_sub[ggdf_noNA_str_sub$type=="TNC",], aes(x = sample_id, y = proportion, fill=cluster)) +
  geom_col(position="stack", width=0.85) +
  ggtitle('TNC Str') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.ticks.x = element_blank(),        
        axis.ticks.y = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white")
  ) +
  ylab("% of Total Cells")+
  facet_grid(.~site2, space = "free", scales = "free_x")+
  scale_fill_manual(values = brewer.set2(4))+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=1))

pdf('Abundance_stackedbar_compact3_noNA_str.pdf', width=9, height=6); 
ggdraw()+
  draw_plot(bp6, x=0, y=0.5, width = 0.6, height = 0.5)+
  draw_plot(bp5, x=0, y=0, width = 1, height = 0.5)
dev.off()

#Stacked bar plots only for str by patients

ggdf_noNA_str_sub <- ggdf_noNA_str_sub %>% group_by(cluster, SPC, met, type, site2) %>% summarise(avg = mean(proportion))
ggdf_noNA_str_sub$SPC <- factor(ggdf_noNA_str_sub$SPC, levels = stragglevels)

bp5 <- ggplot(ggdf_noNA_str_sub[ggdf_noNA_str_sub$type=="LUM",], aes(x = SPC, y = avg, fill=cluster)) +
  geom_col(position="stack", width=0.75) +
  ggtitle('LUM Str') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.ticks.x = element_blank(),        
        axis.ticks.y = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white")
  ) +
  ylab("% of Total Cells")+
  facet_grid(.~site2, space = "free", scales = "free_x")+
  scale_fill_manual(values = brewer.set2(4))+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=1))

bp6 <- ggplot(ggdf_noNA_str_sub[ggdf_noNA_str_sub$type=="TNC",], aes(x = SPC, y = avg, fill=cluster)) +
  geom_col(position="stack", width=0.75) +
  ggtitle('TNC Str') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.ticks.x = element_blank(),        
        axis.ticks.y = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white")
  ) +
  ylab("% of Total Cells")+
  facet_grid(.~site2, space = "free", scales = "free_x")+
  scale_fill_manual(values = brewer.set2(4))+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=1))

pdf('Abundance_stackedbar_compact3_noNA_strbySPC.pdf', width=9, height=6); 
ggdraw()+
  draw_plot(bp6, x=0, y=0.5, width = 0.34, height = 0.5)+
  draw_plot(bp5, x=0, y=0, width = 0.49, height = 0.5)
dev.off()

##Stacked bar plots only for lym

ggdf_noNA$site2 <- ggdf_noNA$site
ggdf_noNA$site2[str_detect(ggdf_noNA$site2, "PBC")] <- "PBC"
ggdf_noNA_lym<-ggdf_noNA[ggdf_noNA$cluster %in% c("Tc","ThN","ThEM","Treg","B_I","B_II","NK_T","NK"),]
bp7 <- ggplot(ggdf_noNA_lym[ggdf_noNA_lym$type=="LUM",], aes(x = sample_id, y = proportion, fill=cluster)) +
  geom_col(position="stack", width=0.85) +
  ggtitle('LUM Lym') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.ticks.x = element_blank(),        
        axis.ticks.y = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white")
  ) +
  ylab("% of Total Cells")+
  facet_grid(.~site2, space = "free", scales = "free_x")+
  scale_fill_manual(values = brewer.dark2(8))+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=1))

bp8 <- ggplot(ggdf_noNA_lym[ggdf_noNA_lym$type=="TNC",], aes(x = sample_id, y = proportion, fill=cluster)) +
  geom_col(position="stack", width=0.85) +
  ggtitle('TNC Lym') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.ticks.x = element_blank(),        
        axis.ticks.y = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white")
  ) +
  ylab("% of Total Cells")+
  facet_grid(.~site2, space = "free", scales = "free_x")+
  scale_fill_manual(values = brewer.dark2(8))+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=1))

pdf('Abundance_stackedbar_compact4_noNA.pdf', width=9, height=6); 
ggdraw()+
  draw_plot(bp8, x=0, y=0.5, width = 0.6, height = 0.5)+
  draw_plot(bp7, x=0, y=0, width = 1, height = 0.5)
dev.off()

## Line plots - proportion taken by celltypes (CK, Str, Others)
# cell counts = 0 samples are removed. 

# proportions for CK 
counts_table_CK <- counts_table[grepl("CK", rownames(counts_table)),]
allCKcounts <- colSums(counts_table_CK)
counts_table_CK <- counts_table_CK[,allCKcounts>0]

allCKcounts <- allCKcounts[allCKcounts != 0]
props_table_CK<-t(t(counts_table_CK) / allCKcounts) * 100
props_CK<- as.data.frame.matrix(props_table_CK)

# proportions for Str
counts_table_Str <- counts_table[grepl("Str", rownames(counts_table)),]
allStrcounts <- colSums(counts_table_Str)
counts_table_Str <- counts_table_Str[,allStrcounts>0]

allStrcounts <- allStrcounts[allStrcounts != 0]
props_table_Str<-t(t(counts_table_Str) / allStrcounts) * 100
props_Str<- as.data.frame.matrix(props_table_Str)

# proportions for others 
counts_table_other <- counts_table[!grepl("Str|CK|Neuron|NA", rownames(counts_table)), ]
allothercounts <- colSums(counts_table_other)
counts_table_other <- counts_table_other[,allothercounts>0]

allothercounts <- allothercounts[allothercounts != 0]
props_table_other<-t(t(counts_table_other) / allothercounts) * 100
props_other<- as.data.frame.matrix(props_table_other)


# CK 
ggdf <- melt(data.frame(cluster = rownames(props_CK),props_CK, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")

ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
# levels not based on ranking the % of total macrophages 
ggdf$SPC <- factor(output$meta_data$Case[match(ggdf$sample_id,output$meta_data$sample_id)], levels = SPClevels)
ggdf$tumor <- factor(output$meta_data$Tumor[match(ggdf$sample_id,output$meta_data$sample_id)], levels=tumorlevels)
ggdf$type <- factor(output$meta_data$Phenotype[match(ggdf$sample_id,output$meta_data$sample_id)], levels=typelevels)
ggdf$site <- factor(output$meta_data$Anno[match(ggdf$sample_id,output$meta_data$sample_id)], levels=sitelevels)
ggdf$TMA <- factor(output$meta_data$TMA[match(ggdf$sample_id,output$meta_data$sample_id)], levels=TMAlevels)
ggdf$met <- factor(output$meta_data$Met[match(ggdf$sample_id,output$meta_data$sample_id)])
ggdf$celltype <- "CK"

ggdf<-ggdf[ggdf$sample_id %nin% samplestoexclude,]

ggdf$site2 <- ggdf$site
ggdf$site2[str_detect(ggdf$site2,"PBC")] <- "PBC"


ggdf_byspc_CK <- ggdf %>% group_by(cluster, SPC, met, type, site2,celltype) %>% summarise(avg = mean(proportion))

# Str
ggdf <- melt(data.frame(cluster = rownames(props_Str),props_Str, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")

ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
# levels not based on ranking the % of total macrophages 
ggdf$SPC <- factor(output$meta_data$Case[match(ggdf$sample_id,output$meta_data$sample_id)], levels = SPClevels)
ggdf$tumor <- factor(output$meta_data$Tumor[match(ggdf$sample_id,output$meta_data$sample_id)], levels=tumorlevels)
ggdf$type <- factor(output$meta_data$Phenotype[match(ggdf$sample_id,output$meta_data$sample_id)], levels=typelevels)
ggdf$site <- factor(output$meta_data$Anno[match(ggdf$sample_id,output$meta_data$sample_id)], levels=sitelevels)
ggdf$TMA <- factor(output$meta_data$TMA[match(ggdf$sample_id,output$meta_data$sample_id)], levels=TMAlevels)
ggdf$met <- factor(output$meta_data$Met[match(ggdf$sample_id,output$meta_data$sample_id)])
ggdf$celltype <- "Str"

ggdf<-ggdf[ggdf$sample_id %nin% samplestoexclude,]

ggdf$site2 <- ggdf$site
ggdf$site2[str_detect(ggdf$site2,"PBC")] <- "PBC"


ggdf_byspc_Str <- ggdf %>% group_by(cluster, SPC, met, type, site2,celltype) %>% summarise(avg = mean(proportion))


# Others
ggdf <- melt(data.frame(cluster = rownames(props_other),props_other, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")

ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
# levels not based on ranking the % of total macrophages 
ggdf$SPC <- factor(output$meta_data$Case[match(ggdf$sample_id,output$meta_data$sample_id)], levels = SPClevels)
ggdf$tumor <- factor(output$meta_data$Tumor[match(ggdf$sample_id,output$meta_data$sample_id)], levels=tumorlevels)
ggdf$type <- factor(output$meta_data$Phenotype[match(ggdf$sample_id,output$meta_data$sample_id)], levels=typelevels)
ggdf$site <- factor(output$meta_data$Anno[match(ggdf$sample_id,output$meta_data$sample_id)], levels=sitelevels)
ggdf$TMA <- factor(output$meta_data$TMA[match(ggdf$sample_id,output$meta_data$sample_id)], levels=TMAlevels)
ggdf$met <- factor(output$meta_data$Met[match(ggdf$sample_id,output$meta_data$sample_id)])
ggdf$celltype <- "Other"

ggdf<-ggdf[ggdf$sample_id %nin% samplestoexclude,]

ggdf$site2 <- ggdf$site
ggdf$site2[str_detect(ggdf$site2,"PBC")] <- "PBC"


ggdf_byspc_Other <- ggdf %>% group_by(cluster, SPC, met, type, site2,celltype) %>% summarise(avg = mean(proportion))


# combine all 
ggdf_byspc<- rbind(ggdf_byspc_CK, ggdf_byspc_Str, ggdf_byspc_Other)
ggdf_byspc$cluster <- factor(ggdf_byspc$cluster, levels = clusterlevels_macs)

library(cowplot)
plots_list <- list()
ggdf_byspc$met<- as.character(ggdf_byspc$met)
ggdf_byspc$met[grepl("0", ggdf_byspc$met)]<-"Primary"
ggdf_byspc$met[grepl("1", ggdf_byspc$met)]<-"Metastatic"
ggdf_byspc$met<- factor(ggdf_byspc$met, level=c('Primary', 'Metastatic'))

for(celltype in unique(ggdf_byspc$celltype)) {
  # Subset the data for the current celltype
  subset_data <- ggdf_byspc[ggdf_byspc$celltype == celltype & ggdf_byspc$type %in% c("TNC", "LUM"),]
  
  p <- ggplot(subset_data, aes(x=type, y=avg, color=met)) +
    geom_boxplot(alpha=0.66, outlier.fill=NA, outlier.color=NA, lwd=0.25) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.25), size=1, alpha=0.5) +
    theme_bw() + 
    theme(axis.text.y = element_text(color="black"),
          axis.ticks.y = element_line(linewidth=0.25),
          axis.line = element_line(linewidth=0.25),
          strip.text=element_text(size=8),
          strip.background = element_rect(fill=NA, color=NA),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.position = "none")+
    facet_wrap(~cluster, scales="free", ncol=6) +
    scale_y_continuous(expand = expansion(mult = c(0.15,0.33))) +
    stat_compare_means(paired=F, 
                       label = "p.signif", 
                       label.x = 1.4, 
                       label.y.npc="top",
                       hide.ns = T,
                       size=5) +
    xlab("") +
    ylab("")

  plots_list[[celltype]] <- p
}


p2<- plot_grid(plots_list[[2]]+ylab("% Stromal cells"),NULL, rel_widths = c(4,2))

combined_plot <- plot_grid(
  plots_list[[1]]+ylab("% CK cells"),
  p2,
  plots_list[[3]]+theme( legend.position = "bottom",
                         legend.title = element_blank(),
                         legend.key.size = unit(.75,'lines'),
                         legend.text = element_text(size=8),
                         legend.key = element_rect(fill="white")) +
    ylab("% Other cells") +
    xlab("type"), 
  nrow = 3,
  axis = 'l',
  rel_heights = c(2,1,3)
) 


pdf('Abundance_boxes_lines.pdf',width=8,height=9);combined_plot;dev.off()


ggdfsub <- ggdf[ggdf$SPC %nin% c("SPC06","SPC09","SPC16","SPC17"),]

ggdfavg <- ggdfsub %>% group_by(cluster, SPC, met) %>% summarise(avg = mean(proportion))

ggp<-ggplot(ggdfavg, aes(x=met, y=avg))+
  geom_line(aes(group=SPC), alpha=0.33, linewidth=0.25)+
  geom_jitter(width=0, size=1, alpha=0.5)+
  geom_boxplot(aes(fill=met), alpha=0.66, outlier.fill=NA, lwd=0.25)+
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color="black"),
        axis.ticks.x = element_blank(),        
        axis.ticks.y = element_line(linewidth=0.25),
        axis.line = element_line(linewidth=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"))+
  facet_wrap(~cluster, ncol=6, scales="free")+
  scale_y_continuous(expand = expansion(mult = c(0.15,0.33)))+
  stat_compare_means(paired=T, 
                     label = "p.signif", 
                     label.x = 1.4, 
                     label.y.npc="top",
                     hide.ns = T,
                     size=5)
pdf('Abundance_lines.pdf',width=8,height=8);ggp;dev.off()

##Lollipop Plots

# for CK+ 
# CK subset
counts_table_CK <- counts_table[grepl("CK", rownames(counts_table)),]
allCKcounts <- colSums(counts_table_CK)
counts_table_CK <- counts_table_CK[,allCKcounts>0]

allCKcounts <- allCKcounts[allCKcounts != 0]
props_table_CK<-t(t(counts_table_CK) / allCKcounts) * 100
props_CK<- as.data.frame.matrix(props_table_CK)


ggdf <- melt(data.frame(cluster = rownames(props_CK),props_CK, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")

ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdf$SPC <- factor(output$meta_data$Case[match(ggdf$sample_id,output$meta_data$sample_id)], levels = SPClevels)
ggdf$tumor <- factor(output$meta_data$Tumor[match(ggdf$sample_id,output$meta_data$sample_id)], levels=tumorlevels)
ggdf$type <- factor(output$meta_data$Phenotype[match(ggdf$sample_id,output$meta_data$sample_id)], levels=typelevels)
ggdf$site <- factor(output$meta_data$Anno[match(ggdf$sample_id,output$meta_data$sample_id)], levels=sitelevels)
ggdf$TMA <- factor(output$meta_data$TMA[match(ggdf$sample_id,output$meta_data$sample_id)], levels=TMAlevels)
ggdf$met <- factor(output$meta_data$Met[match(ggdf$sample_id,output$meta_data$sample_id)])

ggdf<-ggdf[ggdf$sample_id %nin% samplestoexclude,]

ggdfsub <- ggdf[ggdf$met == 0,]
ggdfavg <- ggdfsub %>% group_by(cluster, type) %>% summarise(avg = mean(proportion))
ggdfavg_L <- ggdfavg[ggdfavg$type=="LUM",]
ggdfavg <- ggdfavg[ggdfavg$type=="TNC",]
ggdfavg <- data.frame(cluster=ggdfavg$cluster, LUM=ggdfavg_L$avg, TNC=-ggdfavg$avg)
ggdfavg$diff <- ggdfavg$LUM+ggdfavg$TNC
ggdfavg <- arrange(ggdfavg,desc(diff))
ggdfavg$cluster <- factor(ggdfavg$cluster, levels=ggdfavg$cluster)

lollipop_ck_prim <- ggplot(ggdfavg) +
  ggtitle('CK+ Cell Types
          (Primary)') +
  geom_segment( aes(x=cluster, xend=cluster, y=TNC, yend=LUM), color="grey") +
  geom_point( aes(x=cluster, y=TNC), size=3 , color=hue_pal()(2)[1]) +
  geom_point( aes(x=cluster, y=LUM), size=3 , color=hue_pal()(2)[2]) +
  ylab("-TNC to +LUM") +
  ylim(-25,65) + 
  coord_flip() +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank()
  )

ggdfsub <- ggdf[ggdf$met == 1,]
ggdfavg <- ggdfsub %>% group_by(cluster, type) %>% summarise(avg = mean(proportion))
ggdfavg_L <- ggdfavg[ggdfavg$type=="LUM",]
ggdfavg <- ggdfavg[ggdfavg$type=="TNC",]
ggdfavg <- data.frame(cluster=ggdfavg$cluster, LUM=ggdfavg_L$avg, TNC=-ggdfavg$avg)
ggdfavg$diff <- ggdfavg$LUM+ggdfavg$TNC
ggdfavg <- arrange(ggdfavg,desc(diff))
ggdfavg$cluster <- factor(ggdfavg$cluster, levels=ggdfavg$cluster)

lollipop_ck_met <- ggplot(ggdfavg) +
  ggtitle('CK+ Cell Types
          (Metastatic)') +
  geom_segment( aes(x=cluster, xend=cluster, y=TNC, yend=LUM), color="grey") +
  geom_point( aes(x=cluster, y=TNC), size=3 , color=hue_pal()(2)[1]) +
  geom_point( aes(x=cluster, y=LUM), size=3 , color=hue_pal()(2)[2]) +
  ylab("-TNC to +LUM") +
  ylim(-26,70) + 
  coord_flip() +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank()
  )

pdf('Abundance_lolli_ck_subsetCK.pdf',width=5.5, height=3.5)
ggarrange(lollipop_ck_prim, lollipop_ck_met, common.legend = T)
dev.off()


#for macs
ggdfsub <- ggdf[ggdf$cluster %in% clusterlevels[str_detect(clusterlevels,"Mac")],]
ggdfsub <- ggdfsub[ggdfsub$met == 0,]
ggdfavg <- ggdfsub %>% group_by(cluster, type) %>% summarise(avg = mean(proportion))
ggdfavg_L <- ggdfavg[ggdfavg$type=="LUM",]
ggdfavg <- ggdfavg[ggdfavg$type=="TNC",]
ggdfavg <- data.frame(cluster=ggdfavg$cluster, LUM=ggdfavg_L$avg, TNC=-ggdfavg$avg)
ggdfavg$diff <- ggdfavg$LUM+ggdfavg$TNC
ggdfavg <- arrange(ggdfavg,desc(diff))
ggdfavg$cluster <- factor(ggdfavg$cluster, levels=ggdfavg$cluster)

lollipop_mac_prim <- ggplot(ggdfavg) +
  ggtitle('Macrophage Subtypes
          (Primary)') +
  geom_segment( aes(x=cluster, xend=cluster, y=TNC, yend=LUM), color="grey") +
  geom_point( aes(x=cluster, y=TNC), size=3 , color=hue_pal()(2)[1]) +
  geom_point( aes(x=cluster, y=LUM), size=3 , color=hue_pal()(2)[2]) +
  ylab("-TNC to +LUM") +
  ylim(-5,2.75) + 
  coord_flip() +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank()
  )

ggdfsub <- ggdf[ggdf$cluster %in% clusterlevels[str_detect(clusterlevels,"Mac")],]
ggdfsub <- ggdfsub[ggdfsub$met == 1,]
ggdfavg <- ggdfsub %>% group_by(cluster, type) %>% summarise(avg = mean(proportion))
ggdfavg_L <- ggdfavg[ggdfavg$type=="LUM",]
ggdfavg <- ggdfavg[ggdfavg$type=="TNC",]
ggdfavg <- data.frame(cluster=ggdfavg$cluster, LUM=ggdfavg_L$avg, TNC=-ggdfavg$avg)
ggdfavg$diff <- ggdfavg$LUM+ggdfavg$TNC
ggdfavg <- arrange(ggdfavg,desc(diff))
ggdfavg$cluster <- factor(ggdfavg$cluster, levels=ggdfavg$cluster)

lollipop_mac_met <- ggplot(ggdfavg) +
  ggtitle('Macrophage Subtypes
          (Metastatic)') +
  geom_segment( aes(x=cluster, xend=cluster, y=TNC, yend=LUM), color="grey") +
  geom_point( aes(x=cluster, y=TNC), size=3 , color=hue_pal()(2)[1]) +
  geom_point( aes(x=cluster, y=LUM), size=3 , color=hue_pal()(2)[2]) +
  ylab("-TNC to +LUM") +
  ylim(-5,2.75) + 
  coord_flip() +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank()
  )

pdf('Abundance_lolli_mac.pdf',width=5.2, height=3.5)
ggarrange(lollipop_mac_prim, lollipop_mac_met, common.legend = T)
dev.off()

ggdf_lvt <- ggdf[ggdf$type != "HER2",]
ggdf_lvt <- ggdf_lvt[ggdf_lvt$met==0,]
primaryres <- compare_means(formula= proportion ~ type, 
                            data=ggdf_lvt, 
                            method = "wilcox.test", 
                            paired=F, 
                            group.by = "cluster")
write.csv(primaryres,"Abundance_primary_LvT.csv")

ggdf_lvt <- ggdf[ggdf$type != "HER2",]
ggdf_lvt <- ggdf_lvt[ggdf_lvt$met==1,]
metres <- compare_means(formula= proportion ~ type, 
                        data=ggdf_lvt, 
                        method = "wilcox.test", 
                        paired=F, 
                        group.by = "cluster")
write.csv(metres,"Abundance_met_LvT.csv")

####CELL GEOMETRIC FEATURES PLOTS####

exprall<-fsApply(output$fcs1,exprs)
rng <- colQuantiles(exprall[,c("Area","MinAxis","MajAxis")], probs = c(0.01, 0.95))

exprall<-as.data.frame(exprall)
exprall$sample_ids<-output$sample_ids
exprall$cluster<-factor(output$cell_clustering1m, levels=clusterlevels_macs)

exprall_area<-exprall[exprall$Area<rng["Area","95%"],]
ggpcell_area<-ggplot(exprall_area, aes(x=Area, y=cluster))+
  geom_density_ridges(aes(fill=cluster), scale=5)+
  theme_bw()+
  theme(axis.ticks.y=element_blank(), 
        panel.border=element_blank(), 
        panel.grid = element_blank())+
  scale_fill_manual(values=colorassignedbroad)

exprall_majaxis<-exprall[exprall$MajAxis<rng["MajAxis","95%"],]
ggpcell_majaxis<-ggplot(exprall_majaxis, aes(x=MajAxis, y=cluster))+
  geom_density_ridges(aes(fill=cluster), scale=5)+
  theme_bw()+
  theme(axis.ticks.y=element_blank(), 
        panel.border=element_blank(), 
        panel.grid=element_blank())+
  scale_fill_manual(values=colorassignedbroad)

pdf("Cellgeometricfeatures.pdf",width=5,height=5)
ggpcell_area
ggpcell_majaxis
dev.off()

ggpcellgeom<- exprall %>% group_by(cluster) %>% summarise_at(., vars(Area, MajAxis), list(mean=mean, med=median, sd=sd, se=std.error))

ggpcell_area_dot<-ggplot(ggpcellgeom, aes(x=cluster))+
  geom_errorbar(aes(ymin=Area_med-Area_sd, ymax=Area_med+Area_sd), linewidth=0.25)+
  geom_point(aes(y=Area_med, color=cluster), size=3, show.legend = F)+
  ylab("Area of Each Cell (um)")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        panel.border=element_rect(linewidth=0.25),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(color="black", size=7),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth=0.25))+
  scale_color_manual(values=colorassignedbroad)

ggpcell_majaxis_dot<-ggplot(ggpcellgeom, aes(x=cluster))+
  geom_errorbar(aes(ymin=MajAxis_med-MajAxis_sd, ymax=MajAxis_med+MajAxis_sd), linewidth=0.25)+
  geom_point(aes(y=MajAxis_med, color=cluster), size=3, show.legend = F)+
  ylab("Major Axis of Each Cell (um)")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        panel.border=element_rect(linewidth=0.25),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text = element_text(color="black", size=7),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth=0.25))+
  scale_color_manual(values=colorassignedbroad)

pdf("Cellgeometricfeatures_dots.pdf",width=4.5,height=4)
ggpcell_area_dot
ggpcell_majaxis_dot
dev.off()


