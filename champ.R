library(ChAMP)
library(ChAMPdata)
library(stringr)
library(tidyverse)
library(RColorBrewer)
library(plotrix)
# normalizePath("C:/R-3.6.1/etc/Renviron.site", mustWork = FALSE)
rm(list = ls())
options(stringsAsFactors = F)  
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
testDir=system.file("extdata",package="ChAMPdata")
myLoad <- champ.load(testDir) # Or you may separate about code as champ.import(testDir) + champ.filter()
#CpG.GUI()
champ.QC()
#QC.GUI()
myNorm <- champ.norm()
champ.SVD()
#如果上面SVD Plot显示红色块并不在Sample_Group行，而是较多的在Array行，
#那么则需要对这次数据进行批次校正，使用champ.runCombat函数，
#batchname则是根据上述SVD plot出现的红块位置而定，决定哪个batch factor需要被校正；
#这个过程比较耗时，做完后最好再用champ.SVD检查下
#champ.runCombat()    #If Batch detected, run champ.runCombat() here.

myDMP <- champ.DMP()
save(myDMP,file = "myDMP.rdata")
myDMR <- champ.DMR() 
save(myDMR,file = "myDMR.rdata")
gui <- function (DMR = myDMR, beta = myNorm, pheno = myLoad$pd$Sample_Group, runDMP = F, compare.group = NULL, arraytype = "450K") {
  if (runDMP) {
    tmpbeta <- beta
    tmppheno <- pheno
    if (class(pheno) == "numeric") {
      message("  Your pheno parameter is numeric, champ.DMP() function would calculate linear regression for your CpGs.")
    }
    else {
      message("  You pheno is ", class(pheno), " type.")
      message("    Your pheno information contains following groups. >>")
      sapply(unique(pheno), function(x) message("    <", 
                                                x, ">:", sum(pheno == x), " samples."))
      if (length(unique(pheno)) == 2) {
        message("  Your pheno contains EXACTLY two phenotypes, which is good, compare.group is not needed.")
      }
      else {
        message("  Your pheno contains more than 2 phenotypes, please use compare.group to specify only two of them here.")
        if (is.null(compare.group)) {
          stop("  compare.group is needed here, please specify compare.group.")
        }
        else if (sum(compare.group %in% unique(pheno)) == 
                 2) {
          message("  Your compare.group is in accord with your pheno, which is good, now we are about to extract information for your compare.group.")
          tmpbeta <- beta[, which(pheno %in% compare.group)]
          tmppheno <- pheno[which(pheno %in% compare.group)]
        }
        else {
          stop("  Seems your compare.group is not in accord with your pheno, please recheck your pheno and your compare.group.")
        }
      }
    }
    message("Calculating DMP")
    DMP <- champ.DMP(beta = tmpbeta, pheno = tmppheno, adjPVal = 1, 
                     adjust.method = "BH", compare.group = compare.group, 
                     arraytype = arraytype)
    DMP <- DMP[[1]]
  }   
  data(probe.features)
  probe.features <- probe.features[rownames(beta), ]
  DMR[[1]]$seqnames <- as.factor(substr(DMR[[1]]$seqnames, 4, 100))
  index <- apply(DMR[[1]], 1, function(x) which(probe.features$CHR == x[1] & 
                                                  probe.features$MAPINFO >= as.numeric(x[2]) & 
                                                  probe.features$MAPINFO <= as.numeric(x[3])))
  Anno <- data.frame(DMRindex = unname(unlist(sapply(names(index), function(x) rep(x, length(index[[x]]))))), 
                     probe.features[do.call(c, index), 1:8])

}
x <- gui()
# myBlock <- champ.Block()
# # Block.GUI()
# myGSEA <- champ.GSEA()  
# 
# myEpiMod <- champ.EpiMod()
# myCNA <- champ.CNA()
# #champ.refbase()      #If DataSet is Blood samples, run champ.refbase() here.
# myRefbase <- champ.refbase()
#热图代码
dmr <- as.data.frame(myDMR$BumphunterDMR)
met <- function(y){
  if (y=="hyper") hyper <- subset(dmr,value>0)
  else if (y=="hypo") hyper <- subset(dmr,value<0)
  hyperdata <- subset(x,DMRindex %in% rownames(hyper))
  hyperdata$ID <- rownames(hyperdata)
  hyperdata <- hyperdata[,c("DMRindex","ID")]
  hypernormdata <- as.data.frame(myNorm)
  hypernormdata <- hypernormdata[rownames(hyperdata),]
  hypernormdata$ID <- rownames(hypernormdata)
  hypernormdata <- merge(hypernormdata,hyperdata,by="ID")
  hypernormdata$ID <- str_c(hypernormdata$DMRindex,hypernormdata$ID,sep ="_")
  hypernormdata <- hypernormdata[,-ncol(hypernormdata)]
  rownames(hypernormdata) <- hypernormdata$ID
  hypernormdata <- hypernormdata[,-1]
  hypernormdata <- as.matrix(hypernormdata)
}
hyper <- met("hyper")
hypo <- met("hypo")
pheatmapdata <- rbind(hyper,hypo)
anno <- data.frame(group=as.factor(rep(c("C","T"),each=4)))
rownames(anno) <- colnames(pheatmapdata)
pheatmap::pheatmap(pheatmapdata,cluster_cols = F,cluster_rows = T,scale = "none",show_rownames = F,
                   clustering_distance_rows = "euclidean",
                   # clustering_distance_rows = "correlation",
                   angle_col = 45,annotation_col = anno,annotation_names_col = F,)

#feature分布代码
feature <- function(y){
  if (y=="hyper") type <- subset(dmr,value>0)
  else if (y=="hypo") type <- subset(dmr,value<0)
  data <- subset(x,DMRindex %in% rownames(type))
  data <- data %>% group_by(feature) %>% count()
  data <- as.data.frame(data)
  colnames(data) <- c(str_c(y,"methylated_DMRs"),"Numbers")
  data <- data[order(data$Numbers,decreasing=T),]
  data[,1] <- factor( data[,1] ,levels =  data[,1] )
  data
}
hyperfeature <- feature("hyper")
ggplot (hyperfeature,aes(hypermethylated_DMRs,Numbers))+
  geom_bar(stat = 'identity',width = 0.8,color="black",fill="#11bed9",alpha=0.8,size=0.25)+
  ylab("Numbers of DMRs")+
  xlab("")+
  annotate(label="hypermethylated_DMRs",geom = "text",x=4,y=360,size=5)+
  theme_classic()

hypofeature <- feature("hypo")
ggplot (hypofeature,aes(hypomethylated_DMRs,Numbers))+
  geom_bar(stat = 'identity',width = 0.8,color="black",fill="#11bed9",alpha=0.8,size=0.25)+
  ylab("Numbers of DMRs")+
  xlab("")+
  annotate(label="hypomethylated_DMRs",geom = "text",x=3.5,y=70,size=5)+
  theme_classic()

#cgi分布代码
par(mar=c(1,1,1,1))
cgi<- function(y){
  if (y=="hyper") type <- subset(dmr,value>0)
  else if (y=="hypo") type <- subset(dmr,value<0)
  data <- subset(x,DMRindex %in% rownames(type))
  data <- data %>% group_by(cgi) %>% count()
  data <- as.data.frame(data)
  colnames(data) <- c(str_c(y,"methylated_DMRs"),"Numbers")
  data <- data[order(data$Numbers),]
  labs <- paste0(data[,1],"\n",data[,2])
  # pie(data[,2],labels = labs,init.angle = 90,
  #     col = brewer.pal(nrow(data),"Blues"),
  #     border = "black",radius=1,main = colnames(data)[1]
  #     )
  # dev.off()
  pie3D(data[,2],labels = labs,explode = 0.05,main = colnames(data)[1],radius = 1,start = 0,theta = pi/6,
        col = brewer.pal(nrow(data),"Blues"))
  }
cgi("hyper")
cgi("hypo")


# y$gene <- as.character(y$gene)
y1 <- x[x$gene=="",] %>% distinct(DMRindex) %>% nrow()
y <- x[x$gene!="",]
y <- y %>% 
  group_by(DMRindex) %>% 
  distinct(gene) %>% 
  count() %>% 
  ungroup() %>%
  group_by(n) %>% 
  count() %>% 
  ungroup() %>%
  add_row(n=0,nn=y1)
ggplot(y,aes(as.factor(n),nn/sum(nn)))+
  geom_bar(stat = "identity",width = 0.5,color="black",fill="#11bed9",alpha=0.8,size=0.25)+
  theme_classic()+
  xlab("Number of associated genes per region")+
  ylab("Genomic regions")+
  theme(axis.title.y = element_text())+
  scale_y_continuous(labels = scales::percent)

z <- x[x$gene!="",]
z <- z %>% 
  group_by(DMRindex) %>% 
  distinct(gene) %>%
  ungroup()
library(clusterProfiler)
library(biomaRt)
mart <- useMart("ensembl","hsapiens_gene_ensembl",host = "asia.ensembl.org")
# yyy <- listFilters(mart)
description<- getBM(attributes=c("external_gene_name","entrezgene_id"),
                    filters = 'external_gene_name', values = z$gene, mart = mart)
description <- description[!duplicated(description$external_gene_name),]
description <- description[!duplicated(description$entrezgene_id),]

bp <- enrichGO(gene = description$entrezgene_id,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "BP",
               readable = T)
barplot(bp)
cc <- enrichGO(gene = description$entrezgene_id,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "CC",
               readable = T)
enrichplot::goplot(cc,geom = "label")
mf <- enrichGO(gene = description$entrezgene_id,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "MF",
               readable = T)
enrichplot::cnetplot(mf,circular=T)
kk <- enrichKEGG(gene = description$entrezgene_id,
                 organism = "hsa",
                 keyType = "ncbi-geneid")
kk <- setReadable(kk,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
dotplot(kk)

