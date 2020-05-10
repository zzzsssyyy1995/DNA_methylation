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

#导入数据及低质量样品过滤，归一化
{
  #加载illumina注释文件
  #Illumina注释文件包含450k阵列上每个CpG探针的所有注释信息
  ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  #获取数据目录和加载数据
  dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")#获取数据所在目录
  # list.files(dataDirectory, recursive = TRUE)#列出目录下所有文件（包括子目录中的文件）
  targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheet.csv")# 读入样品信息
  rgSet <- read.metharray.exp(targets=targets)#读r入原始强度数据
  targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
  sampleNames(rgSet) <- targets$ID ##给样本起描述性的名字,即修改样品信息并将新的样品信息附加到原始强度数据中
  #预过滤：过滤低质量样品（注意:未过滤掉低质量的探针)
  detP <- detectionP(rgSet)#计算检测p值
  dim(detP)
  {
    # examine mean detection p-values across all samples to identify any failed samples
    # pal <- brewer.pal(8,"Dark2")
    # par(mfrow=c(1,2))
    # barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
    #         cex.names=0.8, ylab="Mean detection p-values")
    # abline(h=0.05,col="red")
    # legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
    #        bg="white")
    # barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
    #         cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
    # abline(h=0.05,col="red")
    # legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal, 
    #        bg="white")
    # qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group, 
    #          pdf="qcReport.pdf")
    # library(tidyverse)
    # df <- as.data.frame(detP)
    # df <- df %>% mutate(probe_id=rownames(df)) %>% select(probe_id,everything())
    # df_gather <- gather(df,colnames(df)[2:12],key = "sample_type",value = "value")
    # library(ggplot2)
    # ggplot(df_gather,aes(sample_type,value))+
    #   stat_summary(fun.y = "mean",fun.args = list(mult=1),geom = "bar")
    
    # x <- data.frame(colMeans(detP),colSds(detP))
    # colnames(x) <- c("mean","std")
    # x <- x %>% mutate(probe_id=rownames(x)) %>% select(probe_id,everything())
    # ggplot(x,aes(factor(probe_id,levels = probe_id),mean,fill=probe_id))+
    #   geom_bar(stat = "identity")+
    #   theme(axis.text.x = element_text(angle = 44))+
    #   geom_errorbar(data = x,ymin=x$mean,ymax=x$mean+x$std)
  }
  keep <- colMeans(detP) < 0.05
  rgSet <- rgSet[,keep]# remove poor quality samples
  targets <- targets[keep,]# remove poor quality samples from targets data
  detP <- detP[,keep]# remove poor quality samples from detection p-value table
  dim(detP)
  #归一化
  mSetSq <- preprocessQuantile(rgSet) #过滤之后归一化
  {
  # create a MethylSet object from the raw data for plotting
  # mSetRaw <- preprocessRaw(rgSet)#原始数据归一化
  # visualise what the data looks like before and after normalisation
  # par(mfrow=c(1,2))
  # densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
  # legend("top", legend = levels(factor(targets$Sample_Group)), 
  #        text.col=brewer.pal(8,"Dark2"))
  # densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
  #             main="Normalized", legend=FALSE)
  # legend("top", legend = levels(factor(targets$Sample_Group)), 
  #        text.col=brewer.pal(8,"Dark2"))
  }
}


#探索性分析(MDS)
{
  # # MDS plots to look at largest sources of variation
  # par(mfrow=c(1,2))
  # pal <- brewer.pal(8,"Dark2")
  # plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
  #         col=pal[factor(targets$Sample_Group)])
  # legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
  #        bg="white", cex=0.7)
  # plotMDS(getM(mSetSq), top=1000, gene.selection="common",  
  #         col=pal[factor(targets$Sample_Source)])
  # legend("top", legend=levels(factor(targets$Sample_Source)), text.col=pal,
  #        bg="white", cex=0.7)
  # # Examine higher dimensions to look at other sources of variation
  # par(mfrow=c(1,3))
  # plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
  #         col=pal[factor(targets$Sample_Group)], dim=c(1,3))
  # legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal, 
  #        cex=0.7, bg="white")
  # plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
  #         col=pal[factor(targets$Sample_Group)], dim=c(2,3))
  # legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
  #        cex=0.7, bg="white")
  # plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
  #         col=pal[factor(targets$Sample_Group)], dim=c(3,4))
  # legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,
  #        cex=0.7, bg="white")
}


#过滤低质量探针-----探索性分析(MDS)
{
  # ensure probes are in the same order in the mSetSq and detP objects
  detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 
  dim(detP)#未过滤低质量探针之前
  
  keep <- rowSums(detP < 0.01) == ncol(mSetSq) #删除在一个或多个样本中失败的探针
  mSetSqFlt <- mSetSq[keep,]
  dim(mSetSqFlt)
  
  # 如果样本有男有女，移除性染色体上的那些探针
  keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
  table(keep)
  mSetSqFlt <- mSetSqFlt[keep,]
  
  mSetSqFlt <- dropLociWithSnps(mSetSqFlt)# 删除SNP探针中可能影响CpG甲基化的探针
  dim(mSetSqFlt)  
  
  xReactiveProbes <- read.csv(file=paste(dataDirectory,"48639-non-specific-probes-Illumina450k.csv",sep="/"))
  keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)#排除具有交叉反应性的探针
  mSetSqFlt <- mSetSqFlt[keep,] 
  dim(mSetSqFlt)
  
  {
  # par(mfrow=c(1,2))
  # pal <- brewer.pal(8,"Dark2")
  # plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
  #         col=pal[factor(targets$Sample_Group)], cex=0.8)
  # legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
  #        cex=0.65, bg="white")
  # plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
  #         col=pal[factor(targets$Sample_Source)])
  # legend("right", legend=levels(factor(targets$Sample_Source)), text.col=pal,
  #        cex=0.7, bg="white")
  # par(mfrow=c(1,3))
  # # Examine higher dimensions to look at other sources of variation
  # plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
  #         col=pal[factor(targets$Sample_Source)], dim=c(1,3))
  # legend("right", legend=levels(factor(targets$Sample_Source)), text.col=pal,
  #        cex=0.7, bg="white")
  # plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
  #         col=pal[factor(targets$Sample_Source)], dim=c(2,3))
  # legend("topright", legend=levels(factor(targets$Sample_Source)), text.col=pal,
  #        cex=0.7, bg="white")
  # plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
  #         col=pal[factor(targets$Sample_Source)], dim=c(3,4))
  # legend("right", legend=levels(factor(targets$Sample_Source)), text.col=pal,
  #        cex=0.7, bg="white")
  # #calculate M-values for statistical analysis
  # mVals <- getM(mSetSqFlt)
  # head(mVals[,1:5])
  # bVals <- getBeta(mSetSqFlt)
  # head(bVals[,1:5])
  # par(mfrow=c(1,2))
  # densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values",
  #             legend=FALSE, xlab="Beta values")
  # legend("top", legend = levels(factor(targets$Sample_Group)),
  #        text.col=brewer.pal(8,"Dark2"))
  # densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values",
  #             legend=FALSE, xlab="M values")
  # legend("topleft", legend = levels(factor(targets$Sample_Group)),
  #        text.col=brewer.pal(8,"Dark2"))
  }
}


#探针差异甲基化分析,RNA-Seq常见部分
{
  mVals <- getM(mSetSqFlt)
  bVals <- getBeta(mSetSqFlt)
  
  cellType <- factor(targets$Sample_Group)
  individual <- factor(targets$Sample_Source) 
  design <- model.matrix(~0+cellType+individual, data=targets)#创建一个设计矩阵，其中celltype没有参考状态，但individual有参考状态M28
  colnames(design) <- c(levels(cellType),levels(individual)[-1])#去除individual的参考状态M28
  
  fit <- lmFit(mVals, design)#拟合线性模型
  contMatrix <- makeContrasts(naive-rTreg,#为特定的比较创建一个对比矩阵
                              naive-act_naive,
                              rTreg-act_rTreg,
                              act_naive-act_rTreg,
                              levels=design)
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
  DMPsig[duplicated(DMPsig$UCSC_RefGene_Name),]$UCSC_RefGene_Name
  write.csv(DMPsig, file="DMPs.csv", row.names=T)
  
  # plot the top 4 most significantly differentially methylated CpGs 
  # par(mfrow=c(2,2))
  # sapply(rownames(DMPs)[1:4], function(cpg){
  #   plotCpg(bVals, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")
  # })
}


#区域差异甲基化分析
{
  myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M", #获得CpG的相关统计数据
                               analysis.type = "differential", design = design, 
                               contrasts = TRUE, cont.matrix = contMatrix, 
                               coef = "naive - rTreg", arraytype = "450K")
  # str(myAnnotation)
  DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)#使用dmrcate函数将它们组合起来以识别差异甲基化区域
  # options(timeout= 4000000)
  # Sys.setenv(http_proxy="https://127.0.0.1:1080")
  # options(download.file.method="internal")
  # save(DMRs,file="DMRs.rdata")
  # load("DMRs.rdata")
  # save(results.ranges,file="results_ranges.rdata")
  results.ranges <- extractRanges(DMRs)#网络状况要好，不然运行不了
  load("results_ranges.rdata")
  # set up the grouping variables and colours
  pal <- brewer.pal(8,"Dark2")
  groups <- pal[1:length(unique(targets$Sample_Group))]
  names(groups) <- levels(factor(targets$Sample_Group))
  cols <- groups[as.character(factor(targets$Sample_Group))]
  # draw the plot for the top DMR
  par(mfrow=c(1,1))
  DMR.plot(ranges = results.ranges, dmr = 2, CpGs = bVals, phen.col = cols, 
           what = "Beta", arraytype = "450K", genome = "hg19")
}


#自定义甲基化数据的可视化

gen <- "hg19"
dmrIndex <- 1
# extract chromosome number and location from DMR results 
chrom <- as.character(seqnames(results.ranges[dmrIndex]))
start <- as.numeric(start(results.ranges[dmrIndex]))
end <- as.numeric(end(results.ranges[dmrIndex]))
minbase <- start - (0.25*(end-start))
maxbase <- end + (0.25*(end-start))
#添加感兴趣的基因组注释，例如CpG岛的位置和DNAseI超敏感位点；
#CpG岛数据是使用Wu等人发表的方法生成的; 
#DNaseI超敏位点的数据是从UCSC基因组浏览器获得的。
# CpG islands
islandHMM <- read.csv(paste0(dataDirectory,
                             "/model-based-cpg-islands-hg19-chr17.txt"),
                      sep="\t", stringsAsFactors=FALSE, header=FALSE)
head(islandHMM)
islandData <- GRanges(seqnames=Rle(islandHMM[,1]), 
                      ranges=IRanges(start=islandHMM[,2], end=islandHMM[,3]),
                      strand=Rle(strand(rep("*",nrow(islandHMM)))))
islandData
# DNAseI hypersensitive sites
dnase <- read.csv(paste0(dataDirectory,"/wgEncodeRegDnaseClusteredV3chr17.bed"),
                  sep="\t",stringsAsFactors=FALSE,header=FALSE)
head(dnase)
dnaseData <- GRanges(seqnames=dnase[,1],
                     ranges=IRanges(start=dnase[,2], end=dnase[,3]),
                     strand=Rle(rep("*",nrow(dnase))),
                     data=dnase[,5])
dnaseData
iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name="")#设置核型，基因组和RefSeq轨迹，这些轨迹将为我们的甲基化数据提供背景
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq", 
                    from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                    rstarts="exonStarts", rends="exonEnds", gene="name", 
                    symbol="name2", transcript="name", strand="strand", 
                    fill="darkblue",stacking="squish", name="RefSeq", 
                    showId=TRUE, geneSymbol=TRUE)


ann450kOrd <- ann450kSub[order(ann450kSub$chr,ann450kSub$pos),]#确保甲基化数据按染色体和碱基位置排序。
head(ann450kOrd)
bValsOrd <- bVals[match(ann450kOrd$Name,rownames(bVals)),]
head(bValsOrd)
# create genomic ranges object from methylation data
cpgData <- GRanges(seqnames=Rle(ann450kOrd$chr),#使用每种数据类型的适当轨道类型创建数据轨道。
                   ranges=IRanges(start=ann450kOrd$pos, end=ann450kOrd$pos),
                   strand=Rle(rep("*",nrow(ann450kOrd))),
                   betas=bValsOrd)
# extract data on CpGs in DMR
cpgData <- subsetByOverlaps(cpgData, results.ranges[dmrIndex])

# methylation data track
methTrack <- DataTrack(range=cpgData, groups=targets$Sample_Group,genome = gen,
                       chromosome=chrom, ylim=c(-0.05,1.05), col=pal,
                       type=c("a","p"), name="DNA Meth.\n(beta value)",
                       background.panel="white", legend=TRUE, cex.title=0.8,
                       cex.axis=0.8, cex.legend=0.8)
# CpG island track
islandTrack <- AnnotationTrack(range=islandData, genome=gen, name="CpG Is.", 
                               chromosome=chrom,fill="darkgreen")

# DNaseI hypersensitive site data track
dnaseTrack <- DataTrack(range=dnaseData, genome=gen, name="DNAseI", 
                        type="gradient", chromosome=chrom)

# DMR position data track
dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", #设置轨迹列表，并指示不同曲目的相对大小。最后，使用plotTracks函数绘制图
                            chromosome=chrom,fill="darkred")
tracks <- list(iTrack, gTrack, methTrack, dmrTrack, islandTrack, dnaseTrack,
               rTrack)
sizes <- c(2,2,5,2,2,2,3) # set up the relative sizes of the tracks
plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
           add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))



#基因本体测试
{
# Get the significant CpG sites at less than 5% FDR
sigCpGs <- DMPs$Name[DMPs$adj.P.Val<0.05]
sigCpGs[1:10]
length(sigCpGs)
# Get all the CpG sites used in the analysis to form the background
all <- DMPs$Name
length(all)

par(mfrow=c(1,1))
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=TRUE)
topGSA(gst, number=10)
library(tidyverse)
topsa <- subset(gst,FDR<=0.05)
# topsa <- gst %>% group_by(ONTOLOGY) %>% arrange(FDR) %>% slice(1:10)
# topsa$ONTOLOGY <- factor(topsa$ONTOLOGY,levels = c("BP","MF","CC"))
topsa <- arrange(topsa,FDR)
ggplot(data=topsa, aes(x=TERM, y=N, fill=FDR)) +
  geom_bar(stat="identity", width=0.8)+
  labs(x="",y="",title = "")+
  scale_fill_continuous(low="red",high="blue",guide=guide_colorbar(reverse=TRUE))+
  coord_flip()+
  theme_bw(base_size = 15)+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  facet_grid(ONTOLOGY~.,scale="free")

# load Broad human curated (C2) gene sets
load(paste(dataDirectory,"human_c2_v5.rdata",sep="/"))
gsa <- gsameth(sig.cpg=sigCpGs, all.cpg=all, collection=Hs.c2)
topGSA(gsa, number=10)
}





#Differential variability
load(file.path(dataDirectory,"ageData.RData"))
# calculate detection p-values
age.detP <- detectionP(age.rgSet)
# pre-process the data after excluding poor quality samples
age.mSetSq <- preprocessQuantile(age.rgSet)
# add sex information to targets information
age.targets$Sex <- getSex(age.mSetSq)$predictedSex

# ensure probes are in the same order in the mSetSq and detP objects
age.detP <- age.detP[match(featureNames(age.mSetSq),rownames(age.detP)),]
# remove poor quality probes
keep <- rowSums(age.detP < 0.01) == ncol(age.detP) 
age.mSetSqFlt <- age.mSetSq[keep,]

# remove probes with SNPs at CpG or single base extension (SBE) site
age.mSetSqFlt <- dropLociWithSnps(age.mSetSqFlt, snps = c("CpG", "SBE"))

# remove cross-reactive probes
keep <- !(featureNames(age.mSetSqFlt) %in% xReactiveProbes$TargetID)
age.mSetSqFlt <- age.mSetSqFlt[keep,] 
# tag sex chromosome probes for removal
keep <- !(featureNames(age.mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% 
                                                          c("chrX","chrY")])

age.pal <- brewer.pal(8,"Set1")
par(mfrow=c(1,2))
plotMDS(getM(age.mSetSqFlt), top=1000, gene.selection="common", 
        col=age.pal[factor(age.targets$Sample_Group)], labels=age.targets$Sex, 
        main="With Sex CHR Probes")
legend("topleft", legend=levels(factor(age.targets$Sample_Group)), 
       text.col=age.pal)

plotMDS(getM(age.mSetSqFlt[keep,]), top=1000, gene.selection="common", 
        col=age.pal[factor(age.targets$Sample_Group)], labels=age.targets$Sex, 
        main="Without Sex CHR Probes")
legend("top", legend=levels(factor(age.targets$Sample_Group)),
       text.col=age.pal)
# remove sex chromosome probes from data
age.mSetSqFlt <- age.mSetSqFlt[keep,]
# get M-values for analysis
age.mVals <- getM(age.mSetSqFlt)

design <- model.matrix(~factor(age.targets$Sample_Group)) 
# Fit the model for differential variability
# specifying the intercept and age as the grouping factor
fitvar <- varFit(age.mVals, design = design, coef = c(1,2))

# Summary of differential variability
summary(decideTests(fitvar))
topDV <- topVar(fitvar, coef=2)
# Top 10 differentially variable CpGs between old vs. newborns
topDV
# get beta values for ageing data
age.bVals <- getBeta(age.mSetSqFlt)
par(mfrow=c(2,2))
sapply(rownames(topDV)[1:4], function(cpg){
  plotCpg(age.bVals, cpg=cpg, pheno=age.targets$Sample_Group, 
          ylab = "Beta values")
})

#细胞类型组成
# load sorted blood cell data package
library(FlowSorted.Blood.450k)
# ensure that the "Slide" column of the rgSet pheno data is numeric
# to avoid "estimateCellCounts" error
pData(age.rgSet)$Slide <- as.numeric(pData(age.rgSet)$Slide)
# estimate cell counts
cellCounts <- estimateCellCounts(age.rgSet)
# plot cell type proportions by age
par(mfrow=c(1,1))
a = cellCounts[age.targets$Sample_Group == "NewBorns",]
b = cellCounts[age.targets$Sample_Group == "OLD",]
boxplot(a, at=0:5*3 + 1, xlim=c(0, 18), ylim=range(a, b), xaxt="n", 
        col=age.pal[1], main="", ylab="Cell type proportion")
boxplot(b, at=0:5*3 + 2, xaxt="n", add=TRUE, col=age.pal[2])
axis(1, at=0:5*3 + 1.5, labels=colnames(a), tick=TRUE)
legend("topleft", legend=c("NewBorns","OLD"), fill=age.pal)

