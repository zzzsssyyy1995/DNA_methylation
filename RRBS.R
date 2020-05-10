library(methylKit)
# my.methRaw=processBismarkAln(location = "./BMK4R-1_1_bismark_bt2_pe.sorted.bam",
#                              sample.id="test1",
#                              assembly="hg38",
#                              save.context=c("CpG","CHG","CHH"),
#                              read.context = "none",
#                              save.folder=getwd())
file.list=list(system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
               system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
               system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
               system.file("extdata", "control2.myCpG.txt", package = "methylKit"))
myobj=methRead(file.list,
               sample.id=list("test1","test2","ctrl1","ctrl2"),
               assembly="hg18",
               treatment=c(1,1,0,0),
               context="CpG",
               mincov = 10) 
object.size(myobj)
#理想情况下，您应该首先过滤具有极高覆盖率的碱基，以使用filterByCoverage()函数解决PCR偏差，
#然后运行normalizeCoverage()函数以标准化样本之间的覆盖率。
#这两个功能将有助于减少统计测试中由于某些样本读数的系统过采样而可能产生的偏差。
getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)#样品统计信息
getCoverageStats(myobj[[2]],plot=TRUE,both.strands=FALSE)#样品统计信息
myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)
object.size(myobj)
myobj <- normalizeCoverage(myobj)
object.size(myobj)


#样品聚类信息
tiles=tileMethylCounts(myobj,win.size=500,step.size=500)
meth=unite(tiles, destrand=FALSE)#合并样本DMR
# meth=unite(myobj, destrand=FALSE)#合并样本DML
head(meth)
getCorrelation(meth,plot=TRUE)
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
# hc = clusterSamples(meth, dist="correlation", method="ward", plot=F)
PCASamples(meth, screeplot=TRUE)
PCASamples(meth)
##批次效应
# sampleAnnotation=data.frame(batch_id=c("a","a","b","b"),
#                             age=c(19,34,23,40))
# as=assocComp(mBase=meth,sampleAnnotation)
# as
# newObj=removeComp(meth,comp=1)
# (mat=percMethylation(meth))
# mat[mat==100]=80
# newobj=reconstruct(mat,meth)
##平铺窗口分析
# myobj_lowCov = methRead(file.list,
#                         sample.id=list("test1","test2","ctrl1","ctrl2"),
#                         assembly="hg18",
#                         treatment=c(1,1,0,0),
#                         context="CpG",
#                         mincov = 3)
# tiles = tileMethylCounts(myobj_lowCov,win.size=1000,step.size=1000,cov.bases = 10)
# head(tiles[[1]],3)

##查找差异甲基化碱基或区域
myDiff=calculateDiffMeth(meth)
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)
diffMethPerChr(myDiff,plot=T,qvalue.cutoff=0.01, meth.cutoff=25)#可视化每个染色体的低/高甲基化碱基/区域的分布
##过度色散校正
# sim.methylBase1<-dataSim(replicates=6,
#                          sites=1000,
#                          treatment=c(rep(1,3),rep(0,3)),
#                          sample.ids=c(paste0("test",1:3),paste0("ctrl",1:3)))
# my.diffMeth<-calculateDiffMeth(sim.methylBase1[1,],
#                                overdispersion="MN",
#                                test="Chisq",mc.cores=1)
##对协变量做出解释
# covariates=data.frame(age=c(30,80,34,30,80,40))
# sim.methylBase<-dataSim(replicates=6,
#                         sites=1000,
#                         treatment=c(rep(1,3),rep(0,3)),
#                         covariates=covariates,
#                         sample.ids=c(paste0("test",1:3),paste0("ctrl",1:3)))
# my.diffMeth3<-calculateDiffMeth(sim.methylBase,
#                                 covariates=covariates,
#                                 overdispersion="MN",
#                                 test="Chisq",mc.cores=1)
##甲基化位点的基因组分布和cpg岛分布(注释差异甲基化的碱基或区域)
library(genomation)
myDiff_GRanges <- as(myDiff25p,"GRanges")
gene.obj=readTranscriptFeatures("c:/Users/ZSY/Desktop/DNA_methylation/hg38.bed.txt")#读取bed文件，有关启动子/外显子/内含子的注释
diffAnn=annotateWithGeneParts(myDiff_GRanges,gene.obj)# 差异甲基化位点的基因组分布,用启动子、外显子、内含子和基因间区注释给定的对象

cpg.obj=readFeatureFlank("c:/Users/ZSY/Desktop/DNA_methylation/hg38.bed.txt",feature.flank.name=c("CpGi","shores"))#读取CpG岛（CpGi）和CpG海岸注释
diffCpGann=annotateWithFeatureFlank(myDiff_GRanges,cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",
                                    flank.name="shores")#差异甲基化位点的CpG岛分布,函数的作用是:用启动子、外显子、内含子和基因间值注释给定的GRanges对象
##功能注释,在一组定义的区域（例如启动子或CpG岛）上总结甲基化信息
# promoters=regionCounts(myobj,gene.obj$promoters)
# exons=regionCounts(myobj,gene.obj$exons)
# introns=regionCounts(myobj,gene.obj$introns)
# TSSes=regionCounts(myobj,gene.obj$TSSes)

head(getMembers(diffAnn))#获得感兴趣的区域是否与外显子/内含子/启动子重叠
head(getMembers(diffCpGann))#获得感兴趣的区域是否与CpG岛重叠
head(getAssociationWithTSS(diffAnn))#获得距TSS的距离最近的基因名称

getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)#与内含子/外显子/启动子重叠的差异甲基化区域的百分比/数目
plotTargetAnnotation(diffAnn,precedence=TRUE,main="differential methylation annotation")#上面的数据画图
getFeatsWithTargetsStats(diffAnn,percentage=TRUE)#与差异甲基化碱基重叠的内含子/外显子/启动子占所有的内含子/外显子/启动子的百分比
plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),main="differential methylation annotation")#显示CpG岛，CpG岛海岸和其他地区差异甲基化碱基的百分比






library(TxDb.Hsapiens.UCSC.hg38.knownGene)
DMRInfo.ann <- annotateDMRInfo(myDiff_GRanges, 'TxDb.Hsapiens.UCSC.hg38.knownGene')

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)

annoData <- genes(EnsDb.Mmusculus.v79)
annoData;length(annoData)
ranges(annoData)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
txdb_dump <- as.list(txdb)
head(txdb_dump$transcripts)
head(txdb_dump$genes)
