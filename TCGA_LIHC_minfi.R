library(TCGAbiolinks)
library(TCGAWorkflowData)
library(DT)
library(SummarizedExperiment)
library(stringr)
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(wateRmelon)
library(impute)
library(cluster)
rm(list = ls())
options(stringsAsFactors = F)
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#获取miRNA样本名称
# query.miR <- GDCquery(project = CancerProject, 
#                       data.category = "Gene expression",
#                       data.type = "miRNA gene quantification",
#                       file.type = "hg19.mirna",
#                       legacy = TRUE)
#获取DNA甲基化样本名称    
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# lgg.samples <- matchedMetExp("TCGA-LGG", n = 10)
# gbm.samples <- matchedMetExp("TCGA-GBM", n = 10)
# samples <- c(lgg.samples,gbm.samples)
query <- GDCquery(project = c("TCGA-LIHC"),
                  data.category = "DNA methylation",
                  platform = "Illumina Human Methylation 450",
                  legacy = TRUE,
                  sample.type = c("Primary solid Tumor","Solid Tissue Normal"))
GDCdownload(query)
met <- GDCprepare(query, save = T,save.filename = "lihc_dna_met.rda",summarizedExperiment = T)

load("met.rdata")

met1 <- assay(met)
dim(met1)
met1 <- na.omit(met1)# remove probes with NA (similar to na.omit)
dim(met1)
boxplot(met1,col="blue",xaxt="n",outline=F)
met1 <- betaqn(met1)
boxplot(met1,col="blue",xaxt="n",outline=F)
coldata <- data.frame(samples = rownames(colData(met)),Groups=colData(met)$project_id)
densityPlot(met1, sampGroups = coldata$Groups)
densityBeanPlot(met1,sampGroups = coldata$Groups ,sampNames =coldata$samples )
mdsPlot(met1,sampGroups = coldata$Groups,sampNames =coldata$samples,numPositions = 400000)
keep <- !(colnames(met1) %in% c("TCGA-06-5417-01A-01D-1481-05","TCGA-06-0221-02A-11D-2004-05","TCGA-32-1980-01A-01D-1697-05","TCGA-FG-5963-02A-12D-A29T-05"))
met1 <- met1[,keep]
coldata <- coldata[keep,]
mdsPlot(met1,sampGroups = coldata$Groups,sampNames =coldata$samples,numPositions = 400000)


grset <- makeGenomicRatioSetFromMatrix(met1,what = "Beta")
dim(grset)
grset <- dropLociWithSnps(grset)#删除SNP探针中可能影响CpG甲基化的探针
dim(grset)  
keep <- !(featureNames(grset) %in% ann450k$Name[ann450k$chr %in% c("chrNA","chrX","chrY")])
grset <- grset[keep,]# 如果样本有男有女，移除性染色体上的那些探针
dim(grset)
xReactiveProbes <- read.csv(file=paste(getwd(),"48639-non-specific-probes-Illumina450k.csv",sep="/"))
keep <- !(featureNames(grset) %in% xReactiveProbes$TargetID)#排除具有交叉反应性的探针
grset <- grset[keep,] 
dim(grset)

M=getM(grset)
dmp <- dmpFinder(M,pheno = coldata$Groups,type = 'categorical')
dmpdiff <- subset(dmp,qval<=0.05 & abs(intercept)>=2)
dmpdiff$Name <- rownames(dmpdiff)
ann450kdiff <- ann450k[match(rownames(dmpdiff),ann450k$Name), c(1:4,12:19,24:ncol(ann450k))]
dmpdiff <- as.data.frame(merge(dmpdiff,ann450kdiff,by="Name"))

designmatrix <- model.matrix(~factor(coldata$Groups))
colnames(designmatrix)[2] <- c("TCGA-LGG")
dmrs <- bumphunter(grset,design=designmatrix,cutoff=0.25,b=100,type="Beta")
x <- dmrs$table
x$Ref <- 0
x$Alt <- 0
x <- dplyr::select(x,1:3,Ref,Alt)
write.table(x,row.names = F,quote = F,file = '1.txt')
