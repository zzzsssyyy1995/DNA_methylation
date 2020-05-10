library(DSS)
require(bsseq)
path = file.path(system.file(package="DSS"), "extdata")
dat1.1 = read.table(file.path(path, "cond1_1.txt"), header=TRUE)
dat1.2 = read.table(file.path(path, "cond1_2.txt"), header=TRUE)
dat2.1 = read.table(file.path(path, "cond2_1.txt"), header=TRUE)
dat2.2 = read.table(file.path(path, "cond2_2.txt"), header=TRUE)
BSobj = makeBSseqData( list(dat1.1, dat1.2, dat2.1, dat2.2),
                       c("C1","C2", "N1", "N2") )[1:1000,]
BSobj


dmlTest = DMLtest(BSobj, group1=c("C1", "C2"), group2=c("N1", "N2"))
head(dmlTest)
# dmlTest.sm = DMLtest(BSobj, group1=c("C1", "C2"), group2=c("N1", "N2"),
#                      smoothing=TRUE)
dmls = callDML(dmlTest, p.threshold=0.001)
head(dmls)
# dmls2 = callDML(dmlTest, delta=0.1, p.threshold=0.001)
# head(dmls2)
dmrs = callDMR(dmlTest, p.threshold=0.01)
head(dmrs)
# dmrs2 = callDMR(dmlTest, delta=0.1, p.threshold=0.05)
# head(dmrs2)
showOneDMR(dmrs[1,], BSobj)


# Strain = rep(c("A", "B", "C"), 4)
# Sex = rep(c("M", "F"), each=6)
# design = data.frame(Strain,Sex)
# design
# X = model.matrix(~Strain+ Sex, design)
# X
# L = cbind(c(0,1,0,0),c(0,0,1,0))
# L
# matrix(c(0,1,-1,0), ncol=1)


data(RRBS)
RRBS
design
DMLfit = DMLfit.multiFactor(RRBS, design=design, formula=~case+cell+case:cell)
colnames(DMLfit$X)
DMLtest.cell = DMLtest.multiFactor(DMLfit, coef="cellrN")
# DMLtest.cell = DMLtest.multiFactor(DMLfit, coef=3)

# ix=sort(DMLtest.cell[,"pvals"], index.return=TRUE)$ix
# head(DMLtest.cell[ix,])

# DMLtest.cell <- DMLtest.cell[order(DMLtest.cell$pvals,decreasing = F),]
# head(DMLtest.cell)
callDMR(DMLtest.cell, p.threshold=0.05)


DMLfit = DMLfit.multiFactor(RRBS, design, ~case+cell)
colnames(DMLfit$X)
# test1 = DMLtest.multiFactor(DMLfit, coef=2)
test2 = DMLtest.multiFactor(DMLfit, coef="caseSLE")
test3 = DMLtest.multiFactor(DMLfit, term="case")

(Contrast = matrix(c(0,1,0), ncol=1))
test4 = DMLtest.multiFactor(DMLfit, Contrast=Contrast)

cor(cbind(test1$pval, test2$pval, test3$pval, test4$pval))
(Treatment = factor(rep(c("Control","Treated"), 3)))
(pair = factor( rep(1:3, each=2)))
design = data.frame(Treatment, pair)
design
DMLfit = DMLfit.multiFactor(BSobj, design, formula = ~ Treatment + pair)
dmlTest = DMLtest.multiFactor(DMLfit, term="Treatment")
