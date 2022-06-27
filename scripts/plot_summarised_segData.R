library(tidyverse)

if(!dir.exists("data/cellPlotsComp/")){
	dir.create("data/cellPlotsComp/")
}
plot_seg_table <- function(segTable = NULL){
  if(is.null(segTable)){
    stop("SegTable not provided")
  } else if(any(!colnames(segTable) %in% c("chromosome","start","end","segVal","sample","group"))){
    stop("incorrect column names - columns should be 'chromosome','start','end','segVal', and 'sample'")
  }
  segTable$chromosome <- factor(segTable$chromosome,levels = str_sort(unique(segTable$chromosome),numeric = T))
  segTable$segVal[segTable$segVal > 10] <- 10
  for(i in unique(segTable$sample)){
    as.pr <- round(combine_cellLine_stats$purity.x[combine_cellLine_stats$sample == i],digits = 2)
    fix.pr <- round(combine_cellLine_stats$purity.y[combine_cellLine_stats$sample == i],digits = 2)
    p <- ggplot(segTable[segTable$sample == i,]) +
      geom_linerange(aes(xmin = start,xmax = end,y = segVal,color=group),alpha = 0.75,size = 2) +
      scale_y_continuous(breaks = seq.int(0,10,1),
                         labels = seq.int(0,10,1),
                         limits = c(0,10),expand = c(0,0)) +
      labs(title = paste0(i," | ascat pu: ",as.pr," | fixed pu: ",fix.pr)) +
      theme_bw() +
      theme(legend.position = "bottom",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.background = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(colour = "grey90"),
            panel.grid.minor.y = element_blank(),
            panel.spacing.x = unit(0, "lines"),
            strip.background = element_blank()) +
      facet_grid(cols = vars(chromosome),scales = "free_x",space = "free_x",switch = "x")
    png(filename = paste0("data/cellPlotsComp/",i,"_tCN_profile.png"),width = 12,height = 5,units = "in",res = 300)
    print(p)
    dev.off()
  }
}

## Plots below were used during method development of cell line-specific fitting - not required
tcga <- read.table("resources/ASCAT_metadata_TCGA_PCAWG_with_dCIN_extended.txt",header = T)
cellLine <- read.table("data/ascat_results/ascat_profile_statistics.tsv",sep = "\t",header = F)
cellLine_fixed <- read.table("data/ascat_results/ascat_sc_profile_statistics.tsv",sep = "\t",header = F)

tcga <- tcga[,c("patient","purity","ploidy")]
colnames(tcga) <- c("sample","purity","ploidy")
tcga$group <- rep("TCGA",times=nrow(tcga))

colnames(cellLine) <- c("sample","purity","ploidy","GoF")
cellLine <- cellLine[,c("sample","purity","ploidy")]
cellLine$group <- rep("cellLine",times=nrow(cellLine))

colnames(cellLine_fixed) <- c("sample","purity","ploidy")
cellLine_fixed$group <- rep("cellLine_fixed",times=nrow(cellLine_fixed))

plot.d <- rbind(tcga,cellLine,cellLine_fixed)

#pl_pur <- ggplot(plot.d) + geom_point(aes(ploidy,purity,color = group),alpha = 0.5) +
#  geom_smooth(aes(ploidy,purity,color = group),method = "lm") + theme_bw()

#png(filename = "data/ascat_ascatsc_plpur.png",width = 10,height = 8,units = "in",res = 300)
#print(pl_pur)
#dev.off()

#pl_plot <- ggplot(plot.d) +
#  geom_density(aes(ploidy,fill = group),alpha = 0.5) +
#  theme_bw()

#png(filename = "data/ascat_ascatsc_ploidy.png",width = 10,height = 8,units = "in",res = 300)
#print(pl_plot)
#dev.off()

#pur_plot <- ggplot(plot.d) +
#  geom_density(aes(purity,fill = group),alpha = 0.5) +
#  theme_bw()

#png(filename = "data/ascat_ascatsc_purity.png",width = 10,height = 8,units = "in",res = 300)
#print(pur_plot)
#dev.off()

## plot and compare profiles

combine_cellLine_stats <- left_join(cellLine,cellLine_fixed,"sample")
seg_data <- read.table("data/cellLine_segment_ascat_sc_fixed_purity_tCN.tsv",
                          header = TRUE,
                          sep = "\t")
seg_data$group <- rep("ascat.fixed",times=nrow(seg_data))

plot_seg_table(segTable = seg_data)

## Clonality diffs

#seg_data_diffs <- seg_data %>%
#  group_by(sample) %>%
#  mutate(rCN = round(segVal)) %>%
#  mutate(diffs = abs(segVal - rCN)) %>%
#  summarise(across(diffs,mean)) %>%
#  mutate(group = rep("ascat.fixed",times=nrow(.)))

#comb_seg_data_diffs <- seg_data_diffs

#cl_plot <- ggplot(comb_seg_data_diffs) +
  #geom_point(aes(group,diffs),position = "jitter") +
#  labs(title = "Clonality function values",subtitle = "Error function measure of distance from integer state") +
#  geom_violin(aes(group,diffs,fill=group),alpha = 0.5) +
  #ggsignif::geom_signif(aes(group,diffs),comparisons = list(c("ascat","ascat.fixed"),c("ascat.sc","ascat.fixed"))) +
#  theme_bw()

#png(filename = "data/ascat_ascatsc_clonality.png",width = 10,height = 8,units = "in",res = 300)
#print(cl_plot)
#dev.off()

## Construct QC table
plot.list <- list.files(path = "data/cellPlotsComp/",pattern = "*.png")
print(plot.list)
sample.list <- gsub(x = plot.list,pattern = "_tCN_profile.png",replacement = "")
qc.len <- length(sample.list)
qc_table <- data.frame(sample=sample.list,plot=plot.list,use=rep(NA,times=qc.len),notes=rep(NA,times=qc.len))
write.table(x = qc_table,file = "data/cell_fit_qc_table.tsv",quote = F,sep = "\t",row.names = F,col.names = T)
