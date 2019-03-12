# Make consensus TOM for rWGCNA - second script
# Gokul Ramaswami 8-20-2016

# set up R and load libraries
rm(list=ls());
options(stringsAsFactors = FALSE);
library(WGCNA);
library(flashClust);
enableWGCNAThreads(nThreads=12);
set.seed(8675309) # Random seed

# Output directory for figures
outputDir <- "../figures/";
# Output directory for data
outputDirData <- "../processed_data/";

## Set the general parameters
nSets <- 100 # Number of bootstrapped TOMs
prob <- 0.95 ## Define the percentile (called prob) for defining the quantile
nSubset <- 100000 ## number of randomly samples entries of the TOM matrix

## Get the original TOM
load(paste0(outputDirData,"TOM_SP9.RData")) # Load in: TOM from standard pipeline
print("Original TOM loaded")
TOM.mat <- as.matrix(TOM)
nGenes <- nrow(TOM.mat)

## Choose the sampled TOM entries
subsetEntries = sample(nGenes*(nGenes-1)/2, size = nSubset)
TOMsubset = list()

## Load the original TOM
TOMsubset.main <- vectorizeMatrix(TOM.mat)[subsetEntries] ## Select the sampled TOM entries
quantile.TOM.main <- quantile(TOMsubset.main,probs=prob,type = 8) ## Calculate the quantile corresponding to prob

quantile.TOM = rep(1, nSets) ## vector of quantiles of the individual TOM matrices
beta.prob = rep(1, nSets) ## Scaling powers to equalize reference TOM values

print("Scaling for Original TOM Complete !")
## This is a tricky part - we can't load all TOMs simultaneously, so we need to load them piecewise and take the minimum
bsize <- 1000

for (set in 1:nSets) { ## Loop over sets
  if (!file.exists(paste0(outputDirData,"piecewiseTOMs/TOM_",set,"_piece_",ceiling(nGenes/bsize),".RData"))) {
    print(paste0("On TOM ",set,"..."))
    load(paste0(outputDirData,"resampling_TOMs_",set,".RData"))
    tmpTOM <- as.matrix(TOM)
    TOMsubset[[set]] <- vectorizeMatrix(tmpTOM)[subsetEntries] ## Select the sampled TOM entries
    quantile.TOM[set] <- quantile(TOMsubset[[set]],probs=prob,type = 8) ## Calculate the quantile corresponding to prob
    beta.prob[set] = log(quantile.TOM.main)/log(quantile.TOM[set]) ## calculate the power of the adjacency function
    print(paste("Scaling power for this TOM is ",signif(beta.prob[set],3),sep=""))
    tmpTOM <- tmpTOM^(beta.prob[set]) ## use the power adjacency function for calibrating the TOM

    for (i in 1:ceiling(nGenes/bsize)) {
      if (i < ceiling(nGenes/bsize)) {
        subTOM <- tmpTOM[,c((i-1)*bsize+1):(i*bsize)]
        save(subTOM,file=paste0(outputDirData,"piecewiseTOMs/TOM_",set,"_piece_",i,".RData"))
      } else {
        subTOM <- tmpTOM[,c((i-1)*bsize+1):nGenes]
        save(subTOM,file=paste0(outputDirData,"piecewiseTOMs/TOM_",set,"_piece_",i,".RData"))
      }
    }
  } else {
    print(set)
    print("Already processed!")
  }
}

if (!file.exists(paste0(outputDirData,"ConsensusTOMscalingQuantiles.RData"))) {
  save(beta.prob,file=paste0(outputDirData,"ConsensusTOMscalingQuantiles.RData"))
} else {
  load(paste0(outputDirData,"ConsensusTOMscalingQuantiles.RData"))
}

rm(TOM.mat,tmpTOM,subTOM)


for (p in 1:ceiling(nGenes/bsize)) {
  if(!file.exists(paste0(outputDirData,"ConsensusTOM_nSets100_piece_",p,".RData"))) {
  consSubTOM <- vector("list",nSets)
  for (set in 1:nSets) {
    print(paste("loading piece",p,"of TOM",set,sep=" "))
    load(paste0(outputDirData,"piecewiseTOMs/TOM_",set,"_piece_",p,".RData"))
    consSubTOM[[set]]$subTOM <- subTOM
  }
  print("Calculting the median quantile across the subTOM")
  consensusTOM <- matrix(NA,nrow=nrow(subTOM),ncol=ncol(subTOM))
  for (i in 1:nrow(subTOM) ) {
    if (i%%100 == 0) { print(paste("On gene",i,sep=" ")) }
    for (j in 1:ncol(subTOM)) {
      thisvec <- c(consSubTOM[[1]]$subTOM[i,j], consSubTOM[[2]]$subTOM[i,j], consSubTOM[[3]]$subTOM[i,j], consSubTOM[[4]]$subTOM[i,j], consSubTOM[[5]]$subTOM[i,j], consSubTOM[[6]]$subTOM[i,j], consSubTOM[[7]]$subTOM[i,j], consSubTOM[[8]]$subTOM[i,j], consSubTOM[[9]]$subTOM[i,j], consSubTOM[[10]]$subTOM[i,j], consSubTOM[[11]]$subTOM[i,j], consSubTOM[[12]]$subTOM[i,j], consSubTOM[[13]]$subTOM[i,j], consSubTOM[[14]]$subTOM[i,j], consSubTOM[[15]]$subTOM[i,j], consSubTOM[[16]]$subTOM[i,j], consSubTOM[[17]]$subTOM[i,j], consSubTOM[[18]]$subTOM[i,j], consSubTOM[[19]]$subTOM[i,j], consSubTOM[[20]]$subTOM[i,j], consSubTOM[[21]]$subTOM[i,j], consSubTOM[[22]]$subTOM[i,j], consSubTOM[[23]]$subTOM[i,j], consSubTOM[[24]]$subTOM[i,j], consSubTOM[[25]]$subTOM[i,j], consSubTOM[[26]]$subTOM[i,j], consSubTOM[[27]]$subTOM[i,j], consSubTOM[[28]]$subTOM[i,j], consSubTOM[[29]]$subTOM[i,j], consSubTOM[[30]]$subTOM[i,j], consSubTOM[[31]]$subTOM[i,j], consSubTOM[[32]]$subTOM[i,j], consSubTOM[[33]]$subTOM[i,j], consSubTOM[[34]]$subTOM[i,j], consSubTOM[[35]]$subTOM[i,j], consSubTOM[[36]]$subTOM[i,j], consSubTOM[[37]]$subTOM[i,j], consSubTOM[[38]]$subTOM[i,j], consSubTOM[[39]]$subTOM[i,j], consSubTOM[[40]]$subTOM[i,j], consSubTOM[[41]]$subTOM[i,j], consSubTOM[[42]]$subTOM[i,j], consSubTOM[[43]]$subTOM[i,j], consSubTOM[[44]]$subTOM[i,j], consSubTOM[[45]]$subTOM[i,j], consSubTOM[[46]]$subTOM[i,j], consSubTOM[[47]]$subTOM[i,j], consSubTOM[[48]]$subTOM[i,j], consSubTOM[[49]]$subTOM[i,j], consSubTOM[[50]]$subTOM[i,j], consSubTOM[[51]]$subTOM[i,j], consSubTOM[[52]]$subTOM[i,j], consSubTOM[[53]]$subTOM[i,j], consSubTOM[[54]]$subTOM[i,j], consSubTOM[[55]]$subTOM[i,j], consSubTOM[[56]]$subTOM[i,j], consSubTOM[[57]]$subTOM[i,j], consSubTOM[[58]]$subTOM[i,j], consSubTOM[[59]]$subTOM[i,j], consSubTOM[[60]]$subTOM[i,j], consSubTOM[[61]]$subTOM[i,j], consSubTOM[[62]]$subTOM[i,j], consSubTOM[[63]]$subTOM[i,j], consSubTOM[[64]]$subTOM[i,j], consSubTOM[[65]]$subTOM[i,j], consSubTOM[[66]]$subTOM[i,j], consSubTOM[[67]]$subTOM[i,j], consSubTOM[[68]]$subTOM[i,j], consSubTOM[[69]]$subTOM[i,j], consSubTOM[[70]]$subTOM[i,j], consSubTOM[[71]]$subTOM[i,j], consSubTOM[[72]]$subTOM[i,j], consSubTOM[[73]]$subTOM[i,j], consSubTOM[[74]]$subTOM[i,j], consSubTOM[[75]]$subTOM[i,j], consSubTOM[[76]]$subTOM[i,j], consSubTOM[[77]]$subTOM[i,j], consSubTOM[[78]]$subTOM[i,j], consSubTOM[[79]]$subTOM[i,j], consSubTOM[[80]]$subTOM[i,j], consSubTOM[[81]]$subTOM[i,j], consSubTOM[[82]]$subTOM[i,j], consSubTOM[[83]]$subTOM[i,j], consSubTOM[[84]]$subTOM[i,j], consSubTOM[[85]]$subTOM[i,j], consSubTOM[[86]]$subTOM[i,j], consSubTOM[[87]]$subTOM[i,j], consSubTOM[[88]]$subTOM[i,j], consSubTOM[[89]]$subTOM[i,j], consSubTOM[[90]]$subTOM[i,j], consSubTOM[[91]]$subTOM[i,j], consSubTOM[[92]]$subTOM[i,j], consSubTOM[[93]]$subTOM[i,j], consSubTOM[[94]]$subTOM[i,j], consSubTOM[[95]]$subTOM[i,j], consSubTOM[[96]]$subTOM[i,j], consSubTOM[[97]]$subTOM[i,j], consSubTOM[[98]]$subTOM[i,j], consSubTOM[[99]]$subTOM[i,j], consSubTOM[[100]]$subTOM[i,j])
      consensusTOM[i,j] <- median(thisvec)
    }
  }
  save(consensusTOM,file=paste0(outputDirData,"ConsensusTOM_nSets100_piece_",p,".RData"))
  }
}

cat('Done .... \n')

##create consensus TOM
if (!file.exists(paste0(outputDirData,"consensusTOM_final.RData"))) { 
consensusTOM_final=matrix(NA,nrow=nGenes,ncol=1)

for (i in 1:ceiling(nGenes/bsize))
	{
		filename=paste0(outputDirData,"ConsensusTOM_nSets100_piece_",i,".RData",sep="")
		load(filename);
		cat(dim(consensusTOM),'\n')
		consensusTOM_final=cbind(consensusTOM_final,consensusTOM)
		cat(dim(consensusTOM_final),'\n')

	}

	consensusTOM_final= consensusTOM_final[,-1]

dim(consensusTOM_final)  ###should be a square matrix

save(consensusTOM_final,file=paste0(outputDirData,"consensusTOM_final.RData"))
} else {
	load(paste0(outputDirData,"consensusTOM_final.RData"))	
}

####Analysing consensus TOM

if (!file.exists(paste0(outputDirData,"multiData_Resample.RData"))) { 
nSets=100 ### from resampling
setLabels=paste("Data_",c(1:100),sep="")
shortLabels=paste("Data_",c(1:100),sep="")

multiExpr=vector(mode="list",length=nSets)
for (i in 1: nSets){
	load(paste0(outputDirData,"resampling_data_",i,".RData"))
	multiExpr[[i]] = list(data=as.data.frame(datExpr))
	names(multiExpr[[i]]$data)=colnames(datExpr)
	rownames(multiExpr[[i]]$data)=seq(1,nrow(datExpr),1)
	cat('Done for data ....',i,'\n')


}


multiMeta=vector(mode="list",length=nSets)
for (i in 1: nSets){
	load(paste0(outputDirData,"resampling_data_",i,".RData"))
	multiMeta[[i]] = list(data=as.data.frame(meta1))
	names(multiMeta[[i]]$data)=colnames(meta1)
	rownames(multiMeta[[i]]$data)=rownames(meta1)
	cat('Done for data ....',i,'\n')


}


checkSets(multiExpr) # check data size
checkSets(multiMeta) # check data size

rm(datExpr,geneTree,i,ind,merged,meta1)

save(multiExpr,multiMeta,file=paste0(outputDirData,"multiData_Resample.RData"))
}

# Load meta data colors for dendrogram
load(paste0(outputDirData,"DendroTraits_SP9.RData")) # Load in: geneSigsColor

################

load(paste0(outputDirData,"multiData_Resample.RData"))

if (!file.exists(paste0(outputDirData,"consensusTreeCut.RData"))) {
consTOM=consensusTOM_final;

dissTOM = 1-consTOM
consTree = flashClust(as.dist(dissTOM), method = "average");

mColorh <- mLabelh <- colorLabels <- NULL
  for (minModSize in c(20,50,100)) {
    for (dthresh in c(0.1,0.2,0.25)) {
      for (ds in c(2,4)) {
        print("Trying parameters:")
        print(c(minModSize,dthresh,ds))
        tree = cutreeHybrid(dendro = consTree, pamStage=FALSE,
          minClusterSize = minModSize, cutHeight = 0.9999,
          deepSplit = ds, distM = as.matrix(dissTOM))

        merged <- mergeCloseModules(exprData = multiExpr,colors = tree$labels,
                                    cutHeight = dthresh)
        mColorh <- cbind(mColorh,labels2colors(merged$colors))
        mLabelh <- c(mLabelh,paste("DS=",ds," mms=\n",minModSize," dcor=",dthresh))
      }
    }
  }

rm(consTOM, dissTOM)
save(consTree,mColorh,mLabelh,merged,file=paste0(outputDirData,"consensusTreeCut.RData"))
} else {
	load(paste0(outputDirData,"consensusTreeCut.RData"))	
}

###################################


mColorh1=cbind(mColorh,geneSigsColor])

mLabelh1=c(mLabelh,rownames(geneSigsColor))

pdf(paste0(outputDir,"Signed_New_Dendro_Cons.pdf"),height=25,width=20)
plotDendroAndColors(consTree,mColorh1,groupLabels=mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main="Dendrogram With Different Module Cutting Parameters")
dev.off()

if (!file.exists(paste0(outputDirData,"FinalTreeCut.RData"))) {
   ########Final Cut

		mms=50
		ds =4
		dthresh=0.2
        tree = cutreeHybrid(dendro = consTree, pamStage=F, minClusterSize =mms, cutHeight = 0.9999, deepSplit = ds, distM = as.matrix(1-consensusTOM_final))

        merged <- mergeCloseModules(exprData = multiExpr,colors = tree$labels, cutHeight = dthresh)
        mColorh.cons <- labels2colors(merged$colors)
        mLabelh.cons <- "Merged Colors"

	moduleColors.cons = labels2colors(merged$colors);
	mColorh.cons_Label=as.data.frame(cbind(moduleColors.cons,merged$colors))
	rownames(mColorh.cons_Label)=colnames(multiExpr[[1]]$data)

	mColorh2=cbind(mColorh.cons,geneSigsColor)
mLabelh2=c(mLabelh.cons,rownames(geneSigsColor))
save(mColorh2,mLabelh2,mColorh.cons_Label,mms,ds,dthresh,file=paste0(outputDirData,"FinalTreeCut.RData"))
} else {
	load(paste0(outputDirData,"FinalTreeCut.RData"))
}


 pdf(paste0(outputDir,"Final_ModuleDendro_Cons.pdf"),height=10,width=16)
  plotDendroAndColors(consTree, mColorh2, groupLabels = mLabelh2,addGuide=TRUE,dendroLabels=FALSE,main= paste("Consensus Signed bicor network with power = 9, mms=",mms,"ds=",ds,"dthresh=",dthresh));
  dev.off()

##Use the module colors from this consensus network for further analysis and eigengene calculation
# mColorh.cons_Label
