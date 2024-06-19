library(Seurat)
library(ggplot2);library(ggthemes) #for scale_color_tableau
library(ggpubr)
library(dplyr);library(readxl)
library(rstatix)
library(stringr)
library(tidyverse)

###seurat object can be downloaded from GEO with the accession number: GSE269669
avm.aml=readRDS("./avml.aml.rds")
######Supplemental figures
##oncoprint
tp53=read_excel("Research/manuscript/Poonam_AzaVenMargo/Aza Ven Magro Full group _JS_For Dr. Abbas.xlsx")
tp53=as.data.frame(tp53)
tp53=tp53[tp53$MRN %in% metadata$MRN,]
tp53=merge(tp53,as.data.frame(distinct(metadata[,c("MRN","PatientID")],MRN,.keep_all = TRUE)),by="MRN")
tp53$PatientID=substr(tp53$PatientID,1,4);rownames(tp53)=tp53$PatientID
tp53=tp53[paste0("PT",33:43),]
tp53=tp53[,c(73:78,125:135,137:195)]
colnames(tp53)=unlist(lapply(colnames(tp53),function(x){substr(x,5,nchar(x))}))


tp53[is.na(tp53)]="Wild";tp53[1,5]="II D835H"
for (i in 1: nrow(tp53)){
  for (j in 1:ncol(tp53)){
    if (startsWith(tp53[i,j],"I") | startsWith(tp53[i,j],"II")){tp53[i,j]="Mutation"}
    else {tp53[i,j]="Wild"}
  }
}
tp53=tp53[,apply(tp53,2,function(x){sum(x=="Mutation")>0})]  ##remove negative mutations
tp53=tp53[paste0("PT",c(33,35:39,41,42,34,40,43)),]
Heatmap(as.matrix(tp53), col = c("Mutation" = "#C22A38", "Wild"="lightgray"), #brown4
        column_title = "Patients",column_title_gp = gpar(fontsize = 18), row_title=" ", 
        height=unit(2.5,"in"), rect_gp = gpar(col = "white", lwd = 2), show_heatmap_legend = TRUE)



##correlation aml percentage in single cell with BM blast %
df1 = avm.aml@meta.data %>% group_by(orig.ident) %>% summarise(count=n()) %>% as.data.frame() 
df= avm.tme@meta.data %>% group_by(orig.ident) %>% summarise(count=n()) %>% as.data.frame() 
colnames(df1)[2]="count.aml";colnames(df)[2]="count.tme"
df=merge(df,df1,by="orig.ident",all.x = TRUE)
df[is.na(df)]=0
df$ratio=df$count.aml/(df$count.tme+df$count.aml)
df=merge(df,as.data.frame(metadata[,c("PatientID","BM.Blast")]),by.x="orig.ident",by.y="PatientID")
df$BM.Blast=as.numeric(df$BM.Blast);df$ratio=df$ratio * 100

ggscatter(df, x="ratio", y="BM.Blast", add="reg.line",conf.int = FALSE,shape=1,
          cor.coef = TRUE,cor.coef.size = 4,cor.method = "pearson",
          ylab = "%  by Flow Cytometry", xlab = "% by scRNA-seq",size=2) + 
  theme(legend.position = "right", axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))


######connect lines for erythroid proportion for 3 NR and 2 relapse
df=avm.aml@meta.data[,c("orig.ident","response_time","Symphony_Granja2019")]
fraction=df %>% group_by(orig.ident,Symphony_Granja2019)  %>% summarise(count=n()) 
fraction=fraction %>% group_by(orig.ident) %>% mutate(Percent = round(count/sum(count),4))
fraction=as.data.frame(fraction)
fraction$response_time=factor(fraction$response_time,levels = c("CR-Pre","CR-Post","NR-Pre","NR-Post"))
fraction$Symphony_Granja2019=factor(fraction$Symphony_Granja2019,levels=c("HSC-like","CMP.LMPP-like","CLP-like","GMP-like","Lympho-like", "Mono-like", "Baso-like","Eryth-like","DC-like","Unknown"))
fraction=fraction[fraction$orig.ident %in% c("PT34A","PT34B","PT40A","PT40B","PT43A","PT43B","PT35A","PT35C","PT42A","PT42C"),]
fraction=fraction[fraction$Symphony_Granja2019 == "Eryth-like",]
fraction$timepoint=substr(fraction$orig.ident,5,5);fraction$pt=substr(fraction$orig.ident,1,4)
fraction$response=c("NR","NR","Relapse","Relapse","NR","NR","Relapse","Relapse","NR","NR")
fraction=fraction[,c(4:7)]
fraction=rbind(fraction,c(0,"B","PT35","Relapse"));fraction=rbind(fraction,c(0,"B","PT42","Relapse"))
fraction$Percent=as.numeric(fraction$Percent)
ggplot(fraction, aes(x=timepoint, y=Percent, group=pt, color=response)) + 
  geom_line(linetype=1,linewidth=1.5) +geom_point(size=4) + theme_classic() +
  #scale_y_continuous(n.breaks=0.2)+
  theme(legend.position = "top", legend.text = element_text(size = 12),legend.title = element_blank(),
        axis.title.y = element_text(size = 12), axis.title.x = element_blank(),
        axis.text = element_text(size = 12),axis.text.y = element_text(size = 12))


######Symphony projection
###one needs to build the reference file according to the synphony tutorial
plotBasic = function(sobject,                 #symphony result, used to create umap_labels
                     #umap_labels,            # metadata with UMAP labels in UMAP1 and UMAP2 slots
                     ref_path,                     #symphony reference
                     class,                   #column for Symphony classification
                     celltype_col="BioClassification", #colname for celltype in reference metadata
                     facet.by = NULL,         # (optional) metadata column name for faceting
                     color.mapping = NULL,    # custom color mapping
                     legend.position = 'right') {  # Show cell type legend
  
  reference=readRDS(ref_path) 
  expre<- sobject@assays$RNA@counts
  meta=sobject@meta.data
  query <- mapQuery(expre, # query gene expression (genes x cells)
                    meta,  # query metadata (cells x attributes)
                    reference, do_normalize = TRUE, do_umap = TRUE, verbose=TRUE)
  query <- knnPredict(query, reference, reference$meta_data[[celltype_col]], k = 5)#Predict query cell 
  
  reference$meta_data$ref_query = 'reference'
  query$meta_data$ref_query = 'query'
  # Add the UMAP coordinates to the metadata
  #ref.in = cbind(as.data.frame(reference$meta_data)[,c(celltype_col,"ref_query")], reference$umap$embedding)
  query.in = cbind(query$meta_data[,c(class,"ref_query")],query$umap)
  colnames(query.in)[1] <- celltype_col
  p=ggplot()+geom_point(data=ref.in, aes(x=UMAP1, y=UMAP2),color = 'gray',size=0.3)
  p=p+geom_point(data=query.in %>%dplyr::sample_frac(1L),aes(x = UMAP1, y = UMAP2,col = get(celltype_col)),size = 0.6, stroke = 0.4, shape = 16) 
  if (!is.null(color.mapping)) { p = p + scale_color_manual(values = color.mapping) }
  p=p + theme_classic() +labs(color = "Cell type") + theme(legend.position=legend.position) +
    theme(legend.text = element_text(size=22), legend.title=element_text(size=24),
          axis.title = element_text(size=24), axis.text = element_text(size=24)) +
    guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = 'none')
  if(!is.null(facet.by)) {
    p = p + facet_wrap(~get(facet.by)) +
      theme(strip.text.x = element_text(size = 12))}
  return(p)
}







######T cell proportion
df1=avm.combined@meta.data[,c("orig.ident","class2","response_time")];df1=df1[df1$class2!="AML",]
####merge with metadata for response_time column
fraction=table(df1$orig.ident,df1$class2)
fraction <- as.data.frame(fraction/rowSums(fraction))
colnames(fraction)=c("orig.ident","class2","Percent")
fraction=fraction[fraction$class2=="T cell",]
fraction=merge(fraction,metadata[,c("PatientID","response_time")],by.x="orig.ident",by.y="PatientID")
fraction$response_time=factor(fraction$response_time,levels = c("CR-Pre","CR-Post","NR-Pre","NR-Post"))
p=ggplot(fraction,aes(x=response_time,y=Percent))+#scale_y_continuous(trans=squish_trans(0,6),breaks = seq(6,15))+
  geom_boxplot(aes(fill=response_time),width=0.5,color='gray40',outlier.shape=NA,size=0.2)+
  geom_jitter(fill='gray80',color='black',shape=21,width =0.2,size=1.25,stroke = 0.2)+
  theme_classic() + #ylim(0,0.2) +ylab("HALLMARK IFNG RESPONSE") + 
  scale_fill_manual(values = c("#98df8a","#aec7e8","#ff9896","#ffbb78"))+
  theme(legend.position="none",
        axis.title.y = element_text(size = 12),axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12,angle = 90,hjust = 0),axis.text.y = element_text(size = 12))
kruskal.test(fraction$Percent~fraction$response_time)
p + annotate(geom="text",x=2.5,y=0.4,label="p = 0.4724",cex=4)
fraction %>% group_by("response_time") %>% pairwise_wilcox_test(Percent ~ response_time, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>% add_xy_position(x="response_time", fun="max", step.increase = 0.05)















