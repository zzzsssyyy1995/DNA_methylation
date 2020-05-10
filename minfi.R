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


baseDir <- system.file("extdata", package="minfiData")
targets <- read.metharray.sheet(baseDir)
RGSet <- read.metharray.exp(targets = targets)
# phenoData <- pData(RGSet)
# phenoData[,1:6]
# manifest <- getManifest(RGSet)
# manifest
# head(getProbeInfo(manifest))

# MSet <- preprocessRaw(RGSet) 
# MSet
# RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
# RSet
# GRset <- mapToGenome(RSet)
# GRset
# beta <- getBeta(GRset)
# M <- getM(GRset)
# CN <- getCN(GRset)
# gr <- granges(GRset)
# head(gr, n= 3)
# annotation <- getAnnotation(GRset)
# names(annotation)


# qc <- getQC(MSet)
# plotQC(qc)
# densityPlot(MSet, sampGroups = phenoData$Sample_Group)
# densityBeanPlot(MSet, sampGroups = phenoData$Sample_Group)
# controlStripPlot(RGSet, controls="BISULFITE CONVERSION II")
# qcReport(RGSet, pdf= "qcReport.pdf")
# predictedSex <- getSex(GRset, cutoff = -2)$predictedSex
# head(predictedSex)


dim(RGSet)
GRset <- preprocessQuantile(RGSet, fixOutliers = TRUE,
                                     removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                     quantileNormalize = TRUE, stratified = TRUE, 
                                     mergeManifest = FALSE, sex = NULL)
dim(GRset)
GRset <- addSnpInfo(GRset)
GRset <- dropLociWithSnps(GRset, snps=c("SBE","CpG"), maf=0)
dim(GRset)


beta <- getBeta(GRset)
age  <- pData(GRset)$age
dmp <- dmpFinder(beta, pheno = age  , type = "continuous")



pheno <- pData(GRset)$status
designMatrix <- model.matrix(~ pheno)
dmrs <- bumphunter(GRset, design = designMatrix, cutoff = 0.2, B=0, type="Beta")
# library(doParallel)
# registerDoParallel(cores = 3)
dmrs <- bumphunter(GRset, design = designMatrix, cutoff = 0.25, B=100, type="Beta")
names(dmrs)
head(dmrs$table, n=3)

