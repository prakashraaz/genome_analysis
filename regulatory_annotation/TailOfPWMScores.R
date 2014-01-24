quant <- as.numeric(commandArgs(T)[1])
dir <- "/scratch/wynn/FilesForCluster/TFBSOut"
tfs <- c("E-box","LEF1","RUSH")
cutoffs <- c(rep(10,length(tfs)))
for (t in 1:length(tfs)) {
	allsites <- read.table(paste(dir,"/",tfs[t],".pwmscores.txt",sep=""),header=F,as.is=T)
	names(allsites) <- c("TFStart","RefBl","AltBl","RefBr","AltBr")
	cutoffs[t] <- quantile(as.numeric(allsites$RefBr),1-quant)
}
outdf <- data.frame(tfs, cutoffs)
write.table(outdf,file=paste(dir,"/PWMScoreCutoffsFor",quant,"TailOfRefBr.txt",sep=""),row.names=F,quote=F)