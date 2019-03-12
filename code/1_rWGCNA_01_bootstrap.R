# Run rWGCNA - first script bootstrap TOMs
# Gokul Ramaswami 8-20-2016

# set up R and load libraries
rm(list=ls());
options(stringsAsFactors = FALSE);
library(WGCNA);
enableWGCNAThreads(nThreads=12);
library(flashClust);
library(gplots)

# Output directory for figures
outputDir <- "../figures/";
# Output directory for data
outputDirData <- "../processed_data/";

# Load Metadata
metaData <- read.table("../Metadata.txt",sep=" ",header=TRUE)
rownames(metaData) = paste(metaData[,"Braincode"],metaData[,"Region"],sep=".")

# Define NCX + other area regions
regions <-       c("NCX","NCX","NCX","AMY","NCX","NCX","NCX","STR","NCX","HIP","NCX","NCX","NCX","NCX","CBC","MD")
names(regions) = c("S1C","V1C","A1C","AMY","OFC","VFC","MFC","STR","ITC","HIP","STC","M1C","DFC","IPC","CBC","MD")

# Define developmental periods
periods <- c(2,2,3,4,5,5,6,6,6,7,7,8,9,9,10,10,10,10,11,11,12,12,12,12,13,13,13,13,13,14)
names(periods) = c("8PCW","9PCW","12PCW","13PCW","16PCW","17PCW","19PCW","21PCW","22PCW","25PCW","26PCW","4M","6M","10M","12M","2Y","3Y","4Y", "8Y","11Y","13Y","15Y","18Y","19Y","21Y","23Y","30Y","36Y","37Y","40Y") 

# convert age into period in metaData
metaData[,"period"] <- apply(metaData,1,function(x) periods[names(periods)==x["Age"]])

# Load data matrix
dataMatrix <- read.table(gzfile("../data/dataMatrix.txt.gz"),sep="\t",header=TRUE,row.name=1)
# Filter metadata to only rows present in data matrix
metaData <- metaData[colnames(dataMatrix),]
metaData[,"GrossRegion"] <- apply(metaData,1,function(x) regions[names(regions)==x["Region"]])

# Get rid of rows with maximum value of zero
dataMatrix_nozero <- dataMatrix[apply(dataMatrix,1,function(x) max(x[which(!is.na(x))])>0),]

# transpose
dataMatrix_nozero <- t(dataMatrix_nozero)

# use WGCNA function to filter data by missingness entries for both genes and samples
if (!file.exists(paste0(outputDirData,"GSG.RData"))) {
	gsg = goodSamplesGenes(dataMatrix_nozero, minNGenes=1, verbose = 3);
	save(gsg,file=paste0(outputDirData,"GSG.RData"))
} else {
	load(paste0(outputDirData,"GSG.RData"))
}
dataMatrix_nozero_removeOutlier = dataMatrix_nozero[gsg$goodSamples, gsg$goodGenes]

# subset metaData
metaData <- metaData[rownames(dataMatrix_nozero_removeOutlier),]

# look for outliers using connectivity
if (!file.exists(paste0(outputDirData,"zConnectivity.RData"))) {
	sdout <- 2; normadj <- (0.5+0.5*bicor(t(dataMatrix_nozero_removeOutlier), use='pairwise.complete.obs'))^2
	netsummary <- fundamentalNetworkConcepts(normadj); 
	K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
	C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
	outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
	print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(t(dataMatrix_nozero_removeOutlier))[outliers]); print(table(outliers))
	pdf(paste0(outputDir,"Outliers_ZConnect.pdf"))
	plot(Z.K, col = as.numeric(as.factor(metaData[,"GrossRegion"])), pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
	abline(h=-2, lty=2)
	dev.off()
	save(Z.K, outliers, file=paste0(outputDirData,"zConnectivity.RData"))
} else {
	load(paste0(outputDirData,"zConnectivity.RData"))
}
dataMatrix_nozero_removeOutlier_keep <- dataMatrix_nozero_removeOutlier[!outliers,]

# Filter metadata to only rows present in data matrix
metaData = metaData[rownames(dataMatrix_nozero_removeOutlier_keep),]

# Bootstrapping 100 TOM matrices - sample with replacement
# set random seed for sampling
for (i in 1:100) {
	ind = i
	set.seed(ind)
	cat('Random Seed index is :- ',ind,'\n')

	# permute samples within saved groups (prenatal and postnatal in this example)
	# This is relevant if you have a study with 2 groups (e.g. case vs control) - need to balance each sampling by same number of cases and controls
	
	ind_pre=which(metaData$period < 8)
	ind_post=which(metaData$period >= 8)

	Pre_Ind=sample(ind_pre,length(ind_pre),replace=T)
	Post_Ind=sample(ind_post,length(ind_post),replace=T)

	meta_Pre=metaData[Pre_Ind,]
	meta_Post=metaData[Post_Ind,]
	
	meta1=as.data.frame(rbind(meta_Pre,meta_Post))
	# If your study doesn't have saved groups, you can just sample with replacement e.g.:
	# meta1=as.data.frame(metaData[sample(seq(nrow(metaData)),nrow(metaData),replace=T),])
	
	data1=dataMatrix_nozero_removeOutlier_keep[match(paste(meta1[,"Braincode"], meta1[,"Region"],sep="."),rownames(dataMatrix_nozero_removeOutlier_keep)),] # Careful with the sample names - especially regarding duplicates when same sample is kept multiple times in a bootstrap

	# calculate TOM
	datExpr=as.data.frame(data1)
		softPower=9 ## These parameters were previously determined from the initial standard WGCNA run
		ds=4  
    	mms=50
    	dthresh=0.1
	cat('Network Analysis Starting ... \n')
	adjacency = adjacency(datExpr, power = softPower, type = "signed",corFnc="bicor");
	cat('Adjacency Calculation Done ... \n')
 
	TOM = TOMsimilarity(adjacency);
	cat('TOM Calculation Done ... \n')
	dissTOM = 1-TOM
	geneTree = flashClust(as.dist(dissTOM), method = "average");
	filename1=paste0(outputDirData,"resampling_TOMs_",ind,".RData")
	save(TOM,adjacency,file=filename1)  						
	cat('TOM calculation complete for run ',ind,'\n')  	
   						
	tree = cutreeHybrid(dendro = geneTree, pamStage=F, minClusterSize =mms, cutHeight = 0.9999, deepSplit = ds, distM = 	as.matrix(dissTOM))
	merged <- mergeCloseModules(exprData = datExpr,colors = tree$labels, cutHeight = dthresh)   				
	filename2=paste0(outputDirData,"resampling_data_",ind,".RData")
	save(geneTree,merged,datExpr,meta1,ind,file=filename2)				
	cat('Completed for run',ind,'\n')
}
