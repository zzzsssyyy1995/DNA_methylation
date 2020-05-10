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
rm(list = ls())
options(stringsAsFactors = F)
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#加载illumina注释文件
#Illumina注释文件包含450k阵列上每个CpG探针的所有注释信息
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
baseDir <- system.file("extdata", package="minfiData")
targets <- read.metharray.sheet(baseDir)
rgSet <- read.metharray.exp(targets = targets)
rgSet@colData

detP <- detectionP(rgSet)#计算检测p值
dim(detP)
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]# remove poor quality samples
targets <- targets[keep,]# remove poor quality samples from targets data
detP <- detP[,keep]# remove poor quality samples from detection p-value table
dim(detP)

#归一化
mSetSq <- preprocessQuantile(rgSet) #过滤之后归一化

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 
dim(detP)#未过滤低质量探针之前

keep <- rowSums(detP < 0.01) == ncol(mSetSq) #删除在一个或多个样本中失败的探针
mSetSqFlt <- mSetSq[keep,]
dim(mSetSqFlt)

# 如果样本有男有女，移除性染色体上的那些探针
keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% c("chrNA","chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]

mSetSqFlt <- dropLociWithSnps(mSetSqFlt)# 删除SNP探针中可能影响CpG甲基化的探针
dim(mSetSqFlt)  

xReactiveProbes <- read.csv(file=paste(getwd(),"48639-non-specific-probes-Illumina450k.csv",sep="/"))
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)#排除具有交叉反应性的探针
mSetSqFlt <- mSetSqFlt[keep,] 
dim(mSetSqFlt)


mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)

status <- factor(targets$status,levels = c('normal',"cancer"))
design <- model.matrix(~0+status)
colnames(design) <- c(levels(status))

fit <- lmFit(mVals, design)#拟合线性模型
contMatrix <- makeContrasts(cancer-normal,levels=design)
fit2 <- contrasts.fit(fit, contMatrix)#特定比较的拟合
fit2 <- eBayes(fit2)#经验贝叶斯检验，只能调整p值阈值，不能调整lfc
# look at the numbers of DM CpGs at FDR < 0.05
# et <- decideTests(fit2)#差异表达的探针可以从decideTests的结果中快速提取
summary(decideTests(fit2))
# get the table of results for the first contrast (naive - rTreg)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name), c(1:4,12:19,24:ncol(ann450k))]
# DMPs <- topTable(fit2, number =Inf, coef=1, genelist=ann450kSub,p.value = 0.05,lfc = 1,sort.by="p")#提取所有符合条件的DE基因（按调整p值排序）
DMPs <- topTable(fit2, number =Inf, coef=1, genelist=ann450kSub)
DMPsig <- subset(DMPs,adj.P.Val<=0.05)
write.csv(DMPsig, file="DMPs.csv", row.names=T)


myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M", #获得CpG的相关统计数据
                             analysis.type = "differential", design = design,
                             coef = colnames(design)[1], arraytype = "450K")
# str(myAnnotation)
#endif /* NEWSTUFF */
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)#使用dmrcate函数将它们组合起来以识别差异甲基化区域
# options(timeout= 4000000)
# Sys.setenv(http_proxy="https://127.0.0.1:1080")
# options(download.file.method="internal")
# save(DMRs,file="DMRs.rdata")
# load("DMRs.rdata")
# save(results.ranges,file="results_ranges.rdata")
results.ranges <- extractRanges(DMRs)#网络状况要好，不然运行不了
# load("results_ranges.rdata")
# set up the grouping variables and colours
pal <- brewer.pal(8,"Dark2")
groups <- pal[1:length(unique(targets$Sample_Group))]
names(groups) <- levels(factor(targets$Sample_Group))
cols <- groups[as.character(factor(targets$Sample_Group))]
# draw the plot for the top DMR
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges, dmr = 2, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "450K", genome = "hg19")