library(Seurat)
library(ggplot2);library(ggthemes) #for scale_color_tableau
library(ggpubr)
library(dplyr);library(readxl)
library(rstatix)
library(stringr)
library(tidyverse)
library(clusterProfiler);library(enrichplot)


###seurat object can be downloaded from GEO with the accession number: GSE269669
avm.aml=readRDS("./avml_aml.rds")
avm.t=readRDS("./avm_t.rds")
######Main figures
##Figure 4B: dotplot for cell type in all cells
feat_AML=c("CD38","KIT","CD34")
feat_prog=c("CFD","AZU1")
feat_T=c("CD3E","CD8A");feat_B=c("CD19","CD79A","JCHAIN","SDC1");feat_NK=c("NCAM1","NKG7")
feat_mono=c("CD14","ITGAX","CD68","CSF1R","CD1C")
feat_ery=c("HBA1","HBB","GATA1");feat_dc=c("PLD4")
dot_plot_x(avm.combined, rev(features), group = "class2")
ggsave("",width=10,height =15,unit="cm")

allcell_dot <- group_dot_plot(avm.combined, group = "class2",
                              features=list(feat_AML,feat_prog,feat_T,feat_B,feat_NK,feat_mono,feat_ery,feat_dc), 
                              lim = c("AML","Progenitor","T","B/Plasma","NK","Myeloid","Erythroid","pDC"))
leg <- get_legend(dot_plot_leg(avm.combined, feat_AML))

#######Figure 4C: UMAP of all cells
avm.combined=merge(avm.aml,y=avm.tme,add.cell.ids=NULL,merge.data=T)
avm.combined$class2=factor(avm.combined$class2,levels = c("AML","Progenitor","T","B/Plasma","NK","Myeloid","Erythroid", "pDC"))
DimPlot(avm.combined, group.by = "class2",  shuffle=T, label = T) + 
  scale_color_tableau(palette = "Tableau 10") +
  NoAxes()  

#######Figure 4D: Celltype fraction
df=avm.aml@meta.data[,c("orig.ident","response_time","Symphony_Granja2019")]
fraction=df %>% group_by(orig.ident,Symphony_Granja2019)  %>% summarise(count=n()) 
fraction=fraction %>% group_by(orig.ident) %>% mutate(Percent = round(count/sum(count),4))
fraction=as.data.frame(fraction)
fraction$orig.ident=factor(fraction$orig.ident,levels = c("PT33A","PT33B","PT33C","PT35A","PT35B", "PT35C","PT36A","PT36B","PT36C","PT37A","PT37B","PT38A","PT38B","PT39A","PT39B","PT39C","PT41A","PT41B", "PT42A","PT42B","PT42C","PT34A","PT34B","PT40A","PT40B","PT43A","PT43B"))
ggplot(fraction, aes(fill=Celltype, y=Percent, x=orig.ident)) + theme_classic() + xlab("") +
  geom_bar(position="stack",stat="identity",width = 0.8) + labs(fill='Cell type') +
  scale_fill_tableau(palette = "Tableau 10")+
  theme(plot.title = element_text(hjust = 0.5),axis.title.y = element_text(size = 12),
        axis.text.y =element_text(size=12), axis.text.x =element_text(angle = 90,hjust = 0.5,size=12))



#######Figure 4F: LSC17 score in Relapse patients (PT35 and PT42)
df1=avm.aml@meta.data[,c("orig.ident","LSC17")]
df1=df1[df1$orig.ident %in% c("PT34A","PT34B","PT40A","PT40B","PT43A","PT43B","PT35A","PT35C","PT42A","PT42C"),]
df1$orig.ident=factor(df1$orig.ident,levels = c("PT34A","PT34B","PT40A","PT40B","PT43A","PT43B","PT35A","PT35C","PT42A","PT42C"))
p=ggplot(df1,aes(x=orig.ident,y=LSC17))+
  geom_boxplot(aes(fill=orig.ident),width=0.5,color='gray40',outlier.shape=NA,size=0.2)+
  theme_classic() + #ylim(0,0.43) +
  scale_fill_manual(values=c("#dbdb8d", "#dbdb8d", "#c7c7c7", "#c7c7c7", "#9edae5", "#9edae5", "#c5b0d5", "#c5b0d5", "#c49c94", "#c49c94"))+
  theme(legend.position="none",
        axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12,angle = 30,hjust = 1),axis.text.y = element_text(size = 12))
wilcox_stat=df1 %>% group_by("orig.ident") %>% pairwise_wilcox_test(LSC17~orig.ident, p.adjust.method = "BH") %>% add_significance("p.adj") %>% add_xy_position(x="orig.ident", fun="max", step.increase = 0.1)
wilcox_stat=wilcox_stat[c(1,18,31,40,45),];wilcox_stat$y.position=c(0.25,0.2,0.25,0.45,0.35)
p + stat_pvalue_manual(wilcox_stat, label="p.adj.signif", tip.length = 0.01,hide.ns = F,bracket.size = 0.5, size = 4.2)



#######Figure 4H: Pronomoblast percentage 
clinical_table <- read_excel("./AVM Sheet for Yousuf_JS_June 2024.xlsx")

BM_Diff <- read_excel("./BM_Diff_(KS)1128-2023v2_Yousuf_add_2_patients.xlsx")

merged_table <- left_join(BM_Diff,clinical_table, by = "MRN")

merged_table$days_from_RX_date <- as.integer(difftime(merged_table$Date,
                                                      merged_table$RX_Date,
                                                      units = "days"))
merged_table$days_from_DXDate <- as.integer(difftime(merged_table$Date,
                                                     merged_table$DXDate,
                                                     units = "days"))
merged_table$days_from_best_response <- as.integer(difftime(merged_table$Date,
                                                            merged_table$BestRespDate,
                                                            units = "days"))
merged_table$days_from_LossDate <- as.integer(difftime(merged_table$Date,
                                                       merged_table$LossDate,
                                                       units = "days"))
merged_table$days_from_OffDate <- as.integer(difftime(merged_table$Date,
                                                      merged_table$OffDate,
                                                      units = "days"))
merged_table$timeline_sample <- NA
merged_table$patient_group <- NA
for (mrn in merged_table$MRN) {
  baseline_days_from_RX_date <- merged_table %>%
    filter(MRN == mrn & days_from_RX_date <= 0) %>%
    filter(days_from_RX_date == max(days_from_RX_date)) %>%
    select(days_from_RX_date) %>% dplyr::first() %>% as.integer()
  baseline_days_from_DXDate <- merged_table %>%
    filter(MRN == mrn & days_from_DXDate <= 0) %>%
    filter(days_from_DXDate == max(days_from_DXDate)) %>%
    select(days_from_DXDate) %>% dplyr::first() %>% as.integer()
  best_response_days_sample <- merged_table %>%
    filter(MRN == mrn & days_from_best_response <= 0) %>%
    filter(days_from_best_response == max(days_from_best_response)) %>%
    select(days_from_best_response) %>% dplyr::first() %>% as.integer()
  if(!is.na(merged_table %>%
            filter(MRN == mrn)  %>% select(LossDate))[1]){
    LossDate_days_sample <- merged_table %>%
      filter(MRN == mrn)  %>%
      filter(days_from_LossDate <= 0) %>%
      filter(days_from_LossDate == max(days_from_LossDate)) %>%
      select(days_from_LossDate) %>% dplyr::first() %>% as.integer()
  }else {
    LossDate_days_sample <- merged_table %>%
      filter(MRN == mrn)  %>%
      filter(days_from_OffDate <= 0) %>%
      filter(days_from_OffDate == max(days_from_OffDate)) %>%
      select(days_from_OffDate) %>% dplyr::first() %>% as.integer()
  }
  Last_sample <- merged_table %>%
    filter(MRN == mrn) %>%
    filter(days_from_RX_date < 730) %>%
    filter(days_from_RX_date == max(days_from_RX_date)) %>%
    select(days_from_RX_date) %>% dplyr::first() %>% as.integer()
  # mark samples
  merged_table[(merged_table$days_from_RX_date < baseline_days_from_RX_date) &
                 (merged_table$MRN == mrn),]$timeline_sample <- "before_baseline"
  merged_table[(merged_table$days_from_RX_date == baseline_days_from_RX_date) &
                 (merged_table$MRN == mrn),]$timeline_sample <- "baseline" ##Rx
  merged_table[(merged_table$days_from_RX_date > baseline_days_from_RX_date) &
                 (merged_table$MRN == mrn),]$timeline_sample <- "after_baseline"
  merged_table[(merged_table$days_from_RX_date == Last_sample) &
                 (merged_table$days_from_RX_date != baseline_days_from_RX_date) &
                 (merged_table$MRN == mrn),]$timeline_sample <- "last_sample"
  if(unique(merged_table %>%
            filter(MRN == mrn) %>% select(BestResponse)) %in%
     c("NR", "Died")){
    merged_table[(merged_table$MRN == mrn),]$patient_group <- "NR/Died"
    merged_table[(merged_table$days_from_best_response == best_response_days_sample) &
                   (merged_table$MRN == mrn),]$timeline_sample <- "best_response"
    merged_table[(merged_table$days_from_best_response > best_response_days_sample) &
                   (merged_table$MRN == mrn),]$timeline_sample <- "after_best_response"
  }else if(unique(merged_table %>%
                  filter(MRN == mrn) %>% select(BestResponse)) %in%
           c("CR", "CRi", "HI", "PR", "MLFS")){
    if(is.na((merged_table[(merged_table$MRN == mrn),]$days_from_LossDate)[1])){
      if(is.na((merged_table[(merged_table$MRN == mrn),]$days_from_OffDate)[1])){
        merged_table[(merged_table$MRN == mrn),]$patient_group <- "OnProtocol_CR"
        merged_table[(merged_table$days_from_best_response == best_response_days_sample) &
                       (merged_table$MRN == mrn),]$timeline_sample <- "best_response"
        merged_table[(merged_table$days_from_best_response > best_response_days_sample) &
                       (merged_table$MRN == mrn),]$timeline_sample <- "after_best_response"
      }else{
        merged_table[(merged_table$MRN == mrn),]$patient_group <- "CR/CRi/HI/PR/MLFS_OffDate"
        merged_table[(merged_table$MRN == mrn) &
                       (merged_table$days_from_OffDate == LossDate_days_sample),]$timeline_sample <- "relapse"
        merged_table[(merged_table$MRN == mrn) &
                       (merged_table$days_from_OffDate > LossDate_days_sample),]$timeline_sample <- "after_relapse"
      }
    }else{
      merged_table[(merged_table$MRN == mrn),]$patient_group <- "CR/CRi/HI/PR/MLFS_LossDate"
      merged_table[(merged_table$MRN == mrn) &
                     (merged_table$days_from_LossDate == LossDate_days_sample),]$timeline_sample <- "relapse"
      merged_table[(merged_table$MRN == mrn) &
                     (merged_table$days_from_LossDate > LossDate_days_sample),]$timeline_sample <- "after_relapse"
    }
  }
}
merged_table$unique_collection_day <- paste0(merged_table$MRN, "_", merged_table$days_from_RX_date)
merged_table <- merged_table %>%
  distinct(unique_collection_day, .keep_all = TRUE) %>% 
  arrange(MRN, days_from_RX_date)
write.table(merged_table,"./merged_table_cleaned_20240612.txt", row.names = FALSE, sep = "\t")


# Trend of BM blasts % for time points between baseline and relapse/best_response ----
merged_table <- merged_table %>%
  filter(timeline_sample %in% c("baseline", "after_baseline", "relapse", "best_response"))

merged_table[merged_table$TP53 %in% c("Y", "y"),]$TP53 <- "Mutated"
merged_table[merged_table$TP53 %in% c("N"),]$TP53 <- "WT"
merged_table <- as.data.frame(merged_table)

plot_trend <- function(y_axis = "BM_Blasts", group = "TP53"){
  data <- merged_table[!is.na(merged_table[,y_axis]),]
  data <- data[!is.na(data[,"patient_group"]),]
  ggplot(data=data, aes(x=days_from_RX_date,
                        y=data[,y_axis],
                        group= data[,"MRN"])) +
    geom_point(aes(colour=data[,group]), alpha = 0.6) +
    geom_line(aes(colour=data[,group]), alpha = 0.5, linewidth=0.4) +
    # geom_smooth(aes(group=data[,group], colour=data[,group]),
    #             alpha = 0) + #  method = "lm" "gam"
    #facet_wrap(~data[,group], ncol = 1) +
    ylim(min(data[,y_axis],na.rm = TRUE), max(data[,y_axis],na.rm = TRUE)) +
    scale_color_manual(values=c("#CD534CFF", "#0073C2FF")) +
    labs(y= gsub("_", " ", y_axis), x = "Days from RX date") +
    guides(color = guide_legend(title = group)) +
    theme_bw()
}

plot_trend("BM_Blasts")+ theme(legend.position = "none") +
  plot_trend("BM_Pronormoblast") + plot_layout(ncol = 2)



#####################################
##Figure 5A:proportion of each subtype at baseline
df=avm.t@meta.data[,c("orig.ident","BestResponse","time","class7","class8")]
df=df[df$time=="Pre",]
fraction <- table(data.frame(df$orig.ident, df$class7))
fraction <- as.data.frame(fraction/rowSums(fraction))
colnames(fraction) <- c("orig.ident", "class7", "Percent")
fraction=merge(fraction,distinct(df[,c("orig.ident","BestResponse")],orig.ident,.keep_all = TRUE),by="orig.ident")

p=ggplot(fraction,aes(x=class7,y=Percent))+
  geom_boxplot(aes(fill=BestResponse),width=0.5,color='gray40',outlier.shape=NA,size=0.2)+
  #geom_jitter(fill='gray80',color='black',shape=21,width =0,size=2,stroke = 0)+
  theme_classic() + #facet_grid(. ~ class4) +
  scale_fill_manual(values=c("#98df8a","#ff9896"))+
  theme(legend.position="right",
        axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12,angle = 0,hjust = 0.5),axis.text.y = element_text(size = 12))
wilcox_stat_data = fraction %>% group_by(class7) %>% pairwise_wilcox_test(Percent ~ BestResponse, p.adjust.method = "BH") %>% add_significance("p.adj") %>% add_xy_position(x="class7", fun="max", step.increase = 0.1)
p+stat_pvalue_manual(wilcox_stat_data, label="p.adj", tip.length = 0.01)


##Figure 5B:T cell signatures in CD4 and CD8 effector cells
df=avm.t@meta.data[,c("orig.ident","BestResponse","time","class7","class8","CD8_dysfunction","senescence","exhaustion", "CD8_Tumor_reactivity","tumor_specific","NeoTCR_CD8","NeoTCR_CD4")]
df=df[df$time=="Pre",]
df=df[df$class7=="CD4",]

df=df[,c(2,6:12)] %>% pivot_longer(cols = -"BestResponse", names_to="Signature", values_to="values") %>%
  mutate(BestResponse=factor(BestResponse)) %>%
  mutate(Signature=factor(Signature, levels=c("CD8_dysfunction","senescence","exhaustion","CD8_Tumor_reactivity","tumor_specific","NeoTCR_CD8","NeoTCR_CD4"))) %>% as.data.frame()

ggplot(df,aes(x=Signature,y=values)) + 
  geom_boxplot(aes(fill=BestResponse),width=0.5,color='gray40',outlier.shape=NA,size=0.2)+
  theme_classic() + #ylim(0,0.43) +
  scale_fill_manual(values=c("#98df8a","#ff9896"))+ #"#2ca02c", "#d62728"
  theme(legend.position="right",legend.background = element_rect(color="black"),
        axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12,angle = 30,hjust = 1),axis.text.y = element_text(size = 12)) 

##Figure 5C:UMAP. Code same as UMAP for all cells in Figure4

##Figure 5D:Enrichment plot of oxphos
df=tibble::enframe(hallmark) %>%
  dplyr::mutate(value = purrr::map_chr(value, toString)) %>% as.data.frame()
colnames(df)=c("gs_name","gene_symbol")
df=df %>% mutate(gene_symbol = strsplit(as.character(gene_symbol), ", ")) %>% unnest(gene_symbol)
#TERM2GENE <- correlationAnalyzeR::getTERM2GENE(GSEA_Type = c("Hallmark"))
fgseaRes=GSEA(ranks,minGSSize = 1,verbose = FALSE,TERM2GENE = df)
gseaplot2(fgseaRes,geneSetID = "HALLMARK_OXIDATIVE_PHOSPHORYLATION") #,pvalue_table = TRUE

##Figure 5F:Fold change of signatures
df=myeloid@meta.data[,c("orig.ident","response_time","Phago_Zhang","Phagocytosis_regonnition","Phagocytosis_regulation","Phagocytosis_FCGR_signaling","Zhang_m1","Zhang_m2","HLA1","HLA2")]
#df=df[!(df$orig.ident %in% c("PT35C","PT42C")),]
df1=df[df$response_time == "CR-Post", ]#;df1=df1[!(df1$orig.ident %in% c("PT35B","PT42B")),]
df1$relapse=ifelse(df1$orig.ident %in% c("PT35C","PT42C"),"Relapse","Not_relapse")
df1=df1[df1$orig.ident %in% c("PT35B","PT35C","PT42B","PT42C"),]
gm=groupMeans(t(df1[,3:10]),df1$relapse);gm=as.data.frame(gm)
#gm$ratio=gm$Not_relapse/gm$Relapse;gm$name=rownames(gm)
gm$ratio=gm$Relapse/gm$Not_relapse;gm$name=rownames(gm)
gm=gm[order(gm$ratio),]
ggplot(data=gm, aes(x=reorder(name,ratio), y=ratio)) +
  geom_bar(stat="identity",fill="lightblue") +
  ylab("Fold of mean score ") + theme_classic() +
  geom_hline(yintercept = 1,linetype="dotted",linewidth=0.75,color="red") + 
  #geom_errorbar(aes(ymin=nes-error, ymax=nes+error), width=.1) +
  theme(legend.position="none",
        axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12,angle = 90,hjust = 0.5),axis.text.y = element_text(size = 12))


