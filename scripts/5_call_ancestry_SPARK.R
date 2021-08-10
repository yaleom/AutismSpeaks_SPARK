####==========================================================================================####
####========================Step5: call ancestry analysis                   ============================####
####==========================================================================================####

library(parallel)

grmdir="/afm01/UQ/Q3457/SPARK/hrc_imputed/QCd/grm/" 
projc="SPARK" ##cohort to be assigned
r2="0.3" ##r2 threshold used to prune SNPs

refc="1000G"
setwd(grmdir)
## Load PCS
eig <- sqrt( read.table(paste0(refc, "_pruned", r2, ".pca20.eigenval"))[,1] )
ref <- read.table(paste0(refc, "_pruned", r2, ".pca20.eigenvec"), stringsAsFactors=FALSE)
toAssign <- read.table(paste0(refc, "_", projc, "_proj_pruned", r2, ".proj.eigenvec"),stringsAsFactors=FALSE) 
ids <- toAssign[1:2]
##Note the FID is not unique!!!!!!!!

## get parameters
spop  <- read.table("/afm01/UQ/Q3457/AGP/hrc_imputed/QCd/grm/1kg_pop2504.txt",header=TRUE,stringsAsFactor =FALSE) ##note the FID is unique
spop  <- spop[which(spop[,"FID"]%in%ref[,2]),]; rownames(spop) <- spop[,"FID"]
uspop <- unique(spop[,"SuperPOP_CODE"])

nPCS  <- ncol(ref)-2
for(k in 1:nPCS){
  ref[,2+k] <- ref[,2+k] * eig[k]
  toAssign[,2+k] <- toAssign[,2+k] * eig[k]
}

rownames(ref) <- ref[,2]
rownames(toAssign) <- toAssign[,2] ##FID

ref <- as.matrix(ref[,-(1:2)])
toAssign <- as.matrix(toAssign[,-(1:2)])

Pops        <- c("EUR","EAS","SAS","AFR", "AMR")
Cols        <- c("pink","lightblue","lightgreen","grey", "red")
names(Cols) <- Pops


npcToUse           <- 3
meanSPop           <- do.call("cbind",lapply(uspop,function(u) colMeans(ref[which(spop[rownames(ref),"SuperPOP_CODE"]==u),1:npcToUse])))
varSPop            <- do.call("cbind",lapply(uspop,function(u) apply(ref[which(spop[rownames(ref),"SuperPOP_CODE"]==u),1:npcToUse],2,var)))
covSPop            <- lapply(uspop,function(u) solve( cov(ref[which(spop[rownames(ref),"SuperPOP_CODE"]==u),1:npcToUse]) ) )

colnames(meanSPop) <- uspop
colnames(varSPop)  <- uspop
names(covSPop)     <- uspop
detSPop            <- sapply(uspop,function(u) determinant( covSPop[[u]] ,logarithm = TRUE)$modulus)
names(detSPop)     <- uspop

priorProp          <- table(spop[rownames(ref),"SuperPOP_CODE"])
priorProp          <- priorProp/sum(priorProp)

postProb <- function(Coordinates,meanList,covList,detList){
  K <- length(detList)
  z <- Coordinates-meanList
  S <- 0.5*sapply(1:K,function(k) detList[k]-crossprod(z[,k],covList[[k]])%*%z[,k])
  L <- max(S)
  P <- priorProp*exp(S-L)
  P <- P/sum(P)
  names(P) <- names(detList)
  return(P)
}

ntoAssign         <- nrow(toAssign)
predSPoptoAssign  <- do.call("rbind",mclapply(1:ntoAssign,function(i) postProb(Coordinates=toAssign[i,1:npcToUse],meanList=meanSPop,covList=covSPop,detList=detSPop),mc.cores=2))
colnames(predSPoptoAssign) <- paste0("Prob_",colnames(predSPoptoAssign))
rownames(predSPoptoAssign) <- rownames(toAssign)
colSums(predSPoptoAssign)

predSPoptoAssign                   <- as.data.frame(predSPoptoAssign)
predSPoptoAssign$POP_Gauss         <- Pops[apply(predSPoptoAssign[,paste0("Prob_",Pops)],1,which.max)]
predSPoptoAssign$maxPostProb_Gauss <- apply(predSPoptoAssign[,paste0("Prob_",Pops)],1,max)

ancestry_toAssign <- cbind(IID=rownames(predSPoptoAssign),predSPoptoAssign)
ancestry_toAssign$FID <- ids[match(ancestry_toAssign$IID, ids$V2), "V1"]

thrsh <- 0.9
eurid <- ancestry_toAssign[which(ancestry_toAssign$POP_Gauss=="EUR" & predSPoptoAssign$maxPostProb_Gauss>thrsh),c("FID", "IID")]
afrid <- ancestry_toAssign[which(ancestry_toAssign$POP_Gauss=="AFR" & predSPoptoAssign$maxPostProb_Gauss>thrsh),c("FID", "IID")]
easid <- ancestry_toAssign[which(ancestry_toAssign$POP_Gauss=="EAS" & predSPoptoAssign$maxPostProb_Gauss>thrsh),c("FID", "IID")]
sasid <- ancestry_toAssign[which(ancestry_toAssign$POP_Gauss=="SAS" & predSPoptoAssign$maxPostProb_Gauss>thrsh),c("FID", "IID")]
##newadded
amrid <- ancestry_toAssign[which(ancestry_toAssign$POP_Gauss=="AMR" & predSPoptoAssign$maxPostProb_Gauss>thrsh),c("FID", "IID")]


write.table(eurid,"EUR.id",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
write.table(afrid,"AFR.id",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
write.table(easid,"EAS.id",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
write.table(sasid,"SAS.id",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

##new added
write.table(amrid,"AMR.id",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

write.table(ancestry_toAssign, file = "ancestry_assigned.txt",quote=FALSE,row.names=FALSE,col.names=T,sep="\t")

toAssign$Assigned <- ancestry_toAssign[match(toAssign$V2, ancestry_toAssign$IID), "POP_Gauss"]
names(toAssign) <- c("FID", "IID", paste0("PC", 1:20), "Assigned")
toAssign$Group <- "Other"

# phe1 <- read.table("/gpfs1/scratch/group30days/cnsg_park/uqywan/auti/SSC/phenotypes/RRBs/Probands_Key_Variables_matchedID_cleaned.txt", stringsAsFactors = F, header = T) ##2856
# phe2 <- read.table("/gpfs1/scratch/group30days/cnsg_park/uqywan/auti/SSC/phenotypes/RRBs/Unaffected_Siblings_matchedID_cleaned.txt", stringsAsFactors = F, header = T) ##2365
# phe3 <- read.table("/gpfs1/scratch/group30days/cnsg_park/uqywan/auti/SSC/phenotypes/RRBs/Other_Siblings_matchedID_cleaned.txt", stringsAsFactors = F, header = T) ##337
# toAssign[toAssign$IID %in% phe1$IID, "Group"] <- "Probands"
# toAssign[toAssign$IID %in% phe2$IID, "Group"] <- "Unaffected_Siblings"
# toAssign[toAssign$IID %in% phe3$IID, "Group"] <- "Other_Siblings"
# 
# ##probands only
# toAssign1 <- toAssign[toAssign$Group == "Probands",]
# toAssign1$Self_Reported <- phe1[match(toAssign1$IID, phe1$IID), "race"]

write.table(toAssign, file = "SPARK_anc.txt", col.names = T, row.names = F, quote = F)


