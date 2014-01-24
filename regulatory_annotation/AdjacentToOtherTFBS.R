args <- commandArgs(T)
tfname <- args[1]
quant <- as.numeric(args[2])
window <- as.numeric(args[3])
alltfs <- c("E-box","LEF1","RUSH")
poldir <- "/scratch/wynn/FilesForCluster/Pigmentation/Polarizing"
regdir <- "/scratch/wynn/FilesForCluster/TFBSOut"
#results <- read.table(paste(dir,"/",tfname,".results.txt",sep=""),header=F,as.is=T)
#names(results) <- c("TFStart","RefBl","AltBl","RefBr","AltBr","SNPPos")
polresults <- read.table(paste(poldir,"/DerivedInfoFor",tfname,".txt",sep=""),header=T,as.is=T)
#merged <- merge(results,polresults,by="SNPPos")
efder <- polresults[polresults$OneDer=="Ef",]
othertfs <- alltfs[-which(alltfs==tfname)]
TF1 <- TF2 <- c(rep(NA,nrow(efder)))
TF1Pos <- TF2Pos <- c(rep(NA,nrow(efder)))
tf1sites <- read.table(paste(regdir,"/",othertfs[1],".pwmscores.txt",sep=""),header=F,as.is=T)
tf2sites <- read.table(paste(regdir,"/",othertfs[2],".pwmscores.txt",sep=""),header=F,as.is=T)
names(tf1sites) <- names(tf2sites) <- c("TFStart","RefBl","AltBl","RefBr","AltBr")
cutoff1 <- quantile(as.numeric(tf1sites$RefBr),1-quant)
cutoff2 <- quantile(as.numeric(tf2sites$RefBr),1-quant)
tf1sites <- tf1sites[as.numeric(tf1sites$RefBr) >= cutoff1,]
tf2sites <- tf1sites[as.numeric(tf2sites$RefBr) >= cutoff2,]
for (s in 1:nrow(efder)) {
	subtf1 <- tf1sites[intersect(which(as.numeric(tf1sites$TFStart)>=as.numeric(efder$SNPPos[s])-window),which(as.numeric(tf1sites$TFStart)<=as.numeric(efder$SNPPos[s])+window)),]
	subtf2 <- tf2sites[intersect(which(as.numeric(tf2sites$TFStart)>=as.numeric(efder$SNPPos[s])-window),which(as.numeric(tf2sites$TFStart)<=as.numeric(efder$SNPPos[s])+window)),]
	if (nrow(subtf1) > 0) {
		TF1[s] <- othertfs[1]
		TF1Pos[s] <- paste(as.vector(subtf1$TFStart),collapse=",")
                #if (nrow(subtf1) > 1) {
                #    print(paste("Should have multiple sites for ",othertfs[1],"adjacent to ",tfname,sep=""))
                #}
	}
	if (nrow(subtf2) > 0) {
		TF2[s] <- othertfs[2]
		TF2Pos[s] <- paste(as.vector(subtf2$TFStart),collapse=",")
                #if (nrow(subtf2) > 1) {
                #    print(paste("Should have multiple sites for ",othertfs[2],"adjacent to ",tfname,sep=""))
                #}
	}
}
cands <- cbind(efder, TF1, TF1Pos, TF2, TF2Pos)
cands <- cands[union(which(is.na(cands$TF1)==FALSE),which(is.na(cands$TF2)==FALSE)),]
write.table(cands,file=paste(regdir,"/CandidateSitesFor",tfname,"EfDerOthersAdjacent.txt",sep=""),row.names=F,quote=F)
