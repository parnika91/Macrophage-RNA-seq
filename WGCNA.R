library(WGCNA)
library(edgeR)
library(org.Mm.eg.db)
library(pheatmap)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
setwd("~/Documents/Data/Kai/")

SampleCond <- read.csv("~/Documents/Data/Kai/SampleCond.txt", sep="", stringsAsFactors=FALSE)
countfile <- read.delim("~/Documents/Data/Kai/Macrophage_RNAseqReadCount.txt", stringsAsFactors=FALSE, header = T) %>%
  tibble::column_to_rownames("ID")
countfile <- countfile[-grep(rownames(countfile), pattern = "PB"),]
colnames(countfile) <- sapply(colnames(countfile), function(x) strsplit(x, split = "X")[[1]][2])

SampleCond$Condition <- gsub(" ", "_", SampleCond$Condition)
# select the samples that I want
selection <- SampleCond[c(which(SampleCond$Condition==c("untr_4h")), 
                          which(SampleCond$Condition==c("SPZhi_4h")), 
                          which(SampleCond$Condition==c("LPS_4h"))),"Sample"]
countfile <- countfile[,selection]
traitData <- read.csv("~/Documents/Data/Kai/Traits_untr_LPS_SPZhi.csv")

# plot boxplot before and after normalisation
countfile %>% 
  gather(Sample, Count) %>% 
  ggplot(aes(Sample, log10(Count))) + 
  geom_boxplot()

group <- selection
  
# remove space so I can use the names later to makeContrasts
# group <- gsub(" ", "_", group)

# make the count file a DGEList file
y <- DGEList(counts=countfile, group=group)

# remove low expression genes so that FDR depends on fewer genes
keep <- filterByExpr(y) 
y <- y[keep, , keep.lib.sizes=FALSE]


  # before filtering, nrow(y) = 26961
# after filtering, nrow(y) = 13019

# normalise gene counts for lib size and RNA composition effect
y <- calcNormFactors(y, method = "upperquartile")
cou <- cpm(y)
rm(y)
as.data.frame(cou) %>% 
  gather(Sample, Count) %>% 
  ggplot(aes(Sample, log10(Count))) + 
  geom_boxplot()
# getwd()
# workingDir = ".";
# setwd(workingDir)

# options(stringsAsFactors = FALSE)
# femData = read.csv("/home/parnika/Documents/Data/Kai/gene_names_with_norm_count_matrix.csv")
# dim(femData)
# names(femData)

# rownames(femData) <- femData[,1]
# #View(femData)
# femData <- femData[,-1]
# #View(femData)
# colnames(femData) <- sapply(colnames(femData), function(x) strsplit(x, split = "X")[[1]][2])

# datExpr0 = as.data.frame(t(femData[, -c(1:8)]))
# names(datExpr0) = femData$substanceBXH
# rownames(datExpr0) = names(femData)[-c(1:8)]

gsg = goodSamplesGenes(cou, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
datExpr0 <- t(cou)
sampleTree = hclust(dist(datExpr0), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
  
# Plot a line to show the cut
abline(h = 1e06, col = "red")
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 1e06, minSize = 15)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# traitData = SampleCond[c(which(SampleCond$Condition==c("untr_4h")), 
#                          which(SampleCond$Condition==c("SPZhi_4h")), 
#                          which(SampleCond$Condition==c("LPS_4h"))),]
# traitData = read.csv("~/Documents/Data/Kai/Traits_untr_LPS_SPZhi.csv")
# dim(traitData); names(traitData)
# # remove columns that hold information we do not need.
# allTraits = traitData[, -c(31, 16)]
# allTraits = allTraits[, c(2, 11:35)]
# dim(allTraits)
# names(allTraits)
# # Form a data frame analogous to expression data that will hold the clinical traits.
# Samples = rownames(datExpr)
# traitRows = match(femaleSamples, allTraits$Mice)
# datTraits = allTraits[traitRows, -1]
# # rownames(datTraits) = allTraits[traitRows, 1]
# datTraits <- traitData
# collectGarbage()
# 
# # Re-cluster samples
# sampleTree2 = hclust(dist(datExpr), method = "average")
# # Convert traits to a color representation: white means low, red means high, grey means missing entry
# traitColors = numbers2colors(datTraits, signed = FALSE)
#   # Plot the sample dendrogram and the colors underneath.
# plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")
# save(datExpr, datTraits, file = "FemaleLiver-01-dataInput.RData")


enableWGCNAThreads()
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", 
     type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net = blockwiseModules(datExpr, checkMissingData = TRUE, power = 8,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,saveTOMFileBase = "femaleMouseTOM",verbose = 3)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree,file = "FemaleLiver-02-networkConstruction-auto.RData")

plotEigengeneNetworks(
  MEs,
  moduleLabels,
  letterSubPlots = FALSE, Letters = NULL,
  excludeGrey = TRUE, greyLabel = "grey",
  plotDendrograms = TRUE, plotHeatmaps = TRUE,
  setMargins = TRUE, marDendro = NULL, marHeatmap = NULL,
  colorLabels = TRUE, signed = TRUE,
  heatmapColors = NULL,
  plotAdjacency = TRUE,
  printAdjacency = FALSE, cex.adjacency = 0.9,
  coloredBarplot = TRUE, barplotMeans = TRUE, barplotErrors = FALSE,
  plotPreservation = "standard",
  zlimPreservation = c(0, 1),
  printPreservation = FALSE, cex.preservation = 0.9,
)
  


############################################################

# here we define the adjacency matrix using soft thresholding with beta=6
ADJ1=abs(cor(datExpr,use="p"))^6
# When you have relatively few genes (<5000) use the following code
# k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
k=softConnectivity(datE=datExpr,power=6)
# Plot a histogram of k and a scale free topology 
plotsizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

# Adjacency can be used to define a separate measure of similarity, 
# the Topological Overlap Matrix(TOM) 
dissTOM=TOMdist(ADJ1)
collectGarbage()
  
# Turn adjacency into a measure of dissimilarity
dissADJ=1-ADJ1
dissADJ=na.omit(dissADJ)

# pam4=pam(as.dist(dissADJ), 4)
# pam5=pam(as.dist(dissADJ), 5)
# pam6=pam(as.dist(dissADJ), 6)
# # Cross-tabulte the detected and the true (simulated) module membership:
# table(pam4$clustering, truemodule)
# table(pam5$clustering, truemodule)
# table(pam6$clustering, truemodule)  

pamTOM4=pam(as.dist(dissTOM), 4)
pamTOM5=pam(as.dist(dissTOM), 5)
pamTOM6=pam(as.dist(dissTOM), 6)
# Cross-tabulte the detected and the true (simulated) module membership:
table(pamTOM4$clustering, truemodule)
table(pamTOM5$clustering, truemodule)
table(pamTOM6$clustering, truemodule)

# hierADJ=hclust(as.dist(dissADJ), method="average" )
# # Plot the resulting clustering tree together with the true color assignment
# sizeGrWindow(10,5);
# plotDendroAndColors(hierADJ, dendroLabels = FALSE, hang = 0.03,
# main = "Gene hierarchical clustering dendrogram and simulated module colors" )

# Plot results of all module detection methods together:
# Calculate the dendrogram7
hierTOM = hclust(as.dist(dissTOM),method="average");
# The reader should vary the height cut-off parameter h1
# (related to the y-axis of dendrogram) in the following
colorStaticTOM = as.character(WGCNA::cutreeStaticColor(hierTOM, cutHeight=.99, minSize=20))
colorDynamicTOM = WGCNA::labels2colors(dynamicTreeCut::cutreeDynamic(hierTOM,method="tree"))
colorDynamicHybridTOM = WGCNA::labels2colors(dynamicTreeCut::cutreeDynamic(hierTOM, distM= dissTOM , 
                                                    cutHeight = 0.998,deepSplit=2, 
                                                    pamRespectsDendro = FALSE))
# Now we plot the results
sizeGrWindow(10,5)
WGCNA::plotDendroAndColors(hierTOM,colors=data.frame(colorStaticTOM,colorDynamicTOM, colorDynamicHybridTOM),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "Gene dendrogram and module colors, TOM dissimilarity")

colorh1 = colorDynamicHybridTOM
# remove the dissimilarities, adjacency matrices etc to free up space
rm(ADJ1); 
rm(dissADJ);
WGCNA::collectGarbage()
save.image("Simulated-NetworkConstruction.RData")

# To  get  a  sense  of  how  related  the  modules  are  one  can  summarize  
# each  module  by  its  eigengene  (first  principalcomponent).
datME=moduleEigengenes(datExpr,colorh1)$eigengenes
#datME=datME[,-which(colnames(datME)=="MEgrey")]
signif(cor(datME, use="p"), 2)


# We  define  a  dissimilarity  measure  between  the  module  eigengenes  
# that  keeps  track  of  the  sign  of  the  correlation between the module eigengenes,
# and use it to cluster the eigengene:
dissimME=(1-t(cor(datME, method="p")))/2
# dissimME=dissimME[-which(rownames(dissimME)=="MEgrey"),]
# dissimME=dissimME[,-which(colnames(dissimME)=="MEgrey")]

hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")

# We create a pairwise scatter plots of the samples (arrays) along the module eigengenes:
samples = rownames(datME)
y = traitData[sapply(rownames(datME), function(x) which(traitData$Sample==x)),]
sizeGrWindow(8,9)
plotMEpairs(datME, y = y[,2])

#We now create a heatmap plots of module expressions.  
# This may help identify modules that are “held together” by 
# spurious correlations caused by outlying arrays.
sizeGrWindow(8,9)
par(mfrow=c(3,1), mar=c(1, 2, 4, 1))
which.module="turquoise";
plotMat(t(scale(datExpr[,colorh1==which.module ])),
        nrgcols=30,rlabels=T,clabels=T,rcols=which.module,title=which.module)

which.module="blue";
plotMat(t(scale(datExpr[,colorh1==which.module ])),
        nrgcols=30,rlabels=T,clabels=T,rcols=which.module,title=which.module)

which.module="brown";
plotMat(t(scale(datExpr[,colorh1==which.module ])),
        nrgcols=30,rlabels=T,clabels=T,rcols=which.module,title=which.module)

#We  now  create  a  plot  that  explains  the  relationships  between  
# modules  (heatmap)  and  the  corresponding  module eigengene (barplot):
sizeGrWindow(8,7);
which.module="midnightblue"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,ylab="eigengene expression",xlab="array sample")

cor.test(y[,2], datME$MEbrown) # cor, pval of only one module
p.values = corPvalueStudent(cor(y[,2],datME, use="p"), nSamples = length(y[,2])) # cor, pval of all modules

# Measure of module significance as average gene significance
GS1=as.numeric(cor(y[,2],datExpr, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance
ModuleSignificance=tapply(GeneSignificance, colorh1, mean, na.rm=T)

sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,colorh1)

#  Intramodular connectivity
ADJ1=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ1, colorh1)

colorlevels=unique(colorh1)
sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels)))
{
  whichmodule=colorlevels[[i]];
  restrict1 = (colorh1==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],GeneSignificance[restrict1], 
                     col=colorh1[restrict1],main=whichmodule,xlab = "Connectivity", 
                     ylab = "Gene Significance", abline = TRUE)
}

# module membership values kME
datKME=signedKME(datExpr, datME, outputColumnName="MM.")
# Display the first few rows of the data frame
head(datKME)

sizeGrWindow(8,6)
par(mfrow=c(2,2))
# We choose 4 modules to plot: turquoise, blue, brown, green.
# For simplicity we write the code out explicitly for each module.
which.color="turquoise";
restrictGenes=colorh1==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color,xlab="Intramodular Connectivity",ylab="(Module Membership)^6")

which.color="blue";
restrictGenes=colorh1==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color,xlab="Intramodular Connectivity",ylab="(Module Membership)^6")

which.color="brown";
restrictGenes=colorh1==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color,xlab="Intramodular Connectivity",ylab="(Module Membership)^6")

which.color="red";
restrictGenes=colorh1==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color,xlab="Intramodular Connectivity",ylab="(Module Membership)^6")

#  Multi-dimensional scaling plots
# cmd1=cmdscale(as.dist(dissTOM),2)
#     sizeGrWindow(7, 6)
# par(mfrow=c(1,1))
# plot(cmd1, col=as.character(colorh1), main="MDS plot",
#      xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")

#We plot the gene significance against intramodular connectivity: 
#  Finding genes with high gene significance and high intramodular connectivity ininteresting modules
FilterGenes= abs(GS1)> .2 & abs(datKME$MM.brown)>.8
table(FilterGenes)

dimnames(data.frame(datExpr))[[2]][FilterGenes]


#  Relationship between the module membership measures (e.g.  MM.turquoise) andintramodular connectivity
sizeGrWindow(8,6)
par(mfrow=c(2,2))
# We choose 4 modules to plot: turquoise, blue, brown, green.
# For simplicity we write the code out explicitly for each module.
which.color="turquoise";
restrictGenes=colorh1==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color,xlab="Intramodular Connectivity",ylab="(Module Membership)^6")

# Gene screening method based on a detailed definition module membership
NS1=networkScreening(y=y[,2], datME=datME, datExpr=datExpr,oddPower=3, blockSize=1000, 
                     minimumSampleSize=4,addMEy=TRUE, removeDiag=FALSE, weightESy=0.5)
# network screening analysis
mean(NoiseGeneIndicator[rank(NS1$p.Weighted,ties.method="first")<=100])
# standard analysis based on the correlation p-values (or Student T test)
mean(NoiseGeneIndicator[rank(NS1$p.Standard,ties.method="first")<=100])

# Gene Screening based on a gene significance measure
# Perform network screening
NS1GS=networkScreeningGS(datExpr=datExpr, datME = datME, GS=GS1)
# Organize its results for easier plotting
GSprediction1=data.frame(GS1,NS1GS$GS.Weighted)
GS.Weighted=NS1GS$GS.Weighted
# Plot a comparison between standard gene significance and network-weighted gene significance
sizeGrWindow(8, 6)
par(mfrow=c(1,1))
verboseScatterplot(GS1, GS.Weighted,main="Weighted gene significance vs. the standard GS\n",col=colorh1)
abline(0,1)


# Visulaise TOM

power=6
color1=colorh1
restGenes= (color1 != "grey")
diss1=1-TOMsimilarityFromExpr( datExpr[, restGenes], power = 6 )
hier1=hclust(as.dist(diss1), method="average" )
diag(diss1) = NA;
sizeGrWindow(7,7)
TOMplot(diss1^4, hier1, as.character(color1[restGenes]),main = "TOM heatmap plot, module genes" )