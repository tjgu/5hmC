library(parallel)
library(DiffBind)

dir = "/blue/icbrbi/tgu/projects/Jinying_Zhao/sec/tests_v3/data/"
dirres = "/orange/zhao/share/Gu5hmC/bams/q01_v6/"
samples = read.csv(file.path(dir, "sampleSheet_noControl_all_le10M_p5f2_fdr01.txt"), sep="\t")
novo = dba(sampleSheet=samples, minOverlap=200)
dba.save(novo, file="all_noControl_le5M_1060_200.dba", dir=dirres)
info = dba.show(novo)
write.table(info, paste(dirres, "all_noControl_peaksInfo_le5M_1060_200.txt", sep=""), col.names=T, sep="\t", quote=F, row.names=F)
olap.rate <- dba.overlap(novo,mode=DBA_OLAP_RATE)

peakNum = as.numeric(as.character(info[,4]))
tiff(paste(dirres,"all_noControl_peakNumSum_le5M_1060_200.tiff", sep=""), width = 2500, height = 2000, res=300)
hist(peakNum, xlab="peakNumber", main="all samples", breaks=30)
dev.off()

tiff(paste(dirres,"all_noControl_overLapRate_le5M_1060_200.tiff", sep=""), width = 2500, height = 2000, res=300)
plot(olap.rate, type='b', ylab='# peaks', xlab='Overlap at least this many peaksets', main="all samples")
dev.off()
write.table(t(summary(peakNum)), paste(dirres, "all_noControl_peakNumSum_le5M_1060_200.txt", sep=""), sep="\t", quote=F, row.names=F)
write.table(olap.rate, paste(dirres, "all_noControl_overLapRate_le5M_1060_200.txt", sep=""), sep="\t", quote=F, col.names=F)

consensus.peaks <- dba.peakset(novo, bRetrieve=TRUE)
write.table(consensus.peaks, paste(dirres, "all_noControl_consensusPeaks_le5M_1060_200.txt", sep=""), quote=F, sep="\t")

novoLen = summary(novo$merged[,3] - novo$merged[,2] + 1)
consensus.peaks.df = as.data.frame(consensus.peaks)
novoLenF50 = summary(consensus.peaks.df$end - consensus.peaks.df$start + 1)
novoLenAll = rbind(novoLen, novoLenF50)
write.table(novoLenAll, paste(dirres, "all_noControl_consensusPeaks_lengthSummary_le5M_1060_200.txt", sep=""), quote=F, sep="\t")

novo2 <- dba.count(novo, peaks=consensus.peaks, summits=FALSE, bSubControl=FALSE, score=DBA_SCORE_RPKM, bParallel=FALSE, bUseSummarizeOverlaps=TRUE, mapQCth=20)
info = dba.show(novo2)
dba.save(novo2, file="all_noControl_le5M_1060_200.dba", dir=dirres)

readsy = NULL
for(i in 1:length(novo2$peaks)){
        readso = novo2$peaks[[i]]$Reads
        readsy = cbind(readsy, readso)
}
colnames(readsy) = novo2$samples$SampleID
readsy = cbind(novo2$peaks[[1]]$Chr, novo2$peaks[[1]]$Start, novo2$peaks[[1]]$End, readsy)
colnames(readsy)[1:3] = c("Chr", "Start", "End")
write.table(readsy, paste(dirres, "all_noControl_consensusPeaks_reCounts_reads_le5M_1060_200.txt", sep=""), quote=F, sep="\t")
gc()

libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP, PeakReads=round(info$Reads*info$FRiP))
rownames(libsizes) <- info$ID
write.table(libsizes, paste(dirres, "all_noControl_consensusPeaks_readsRecount_info_le5M_1060_200.txt", sep=""), quote=F, sep="\t", row.names=T)

reCountPeaks <- dba.peakset(novo2, bRetrieve=TRUE)
write.table(reCountPeaks, paste(dirres, "all_noControl_consensusPeaks_reCounts_le5M_1060_200.txt", sep=""), quote=F, sep="\t")

dba.save(novo2, file="all_noControl_le5M_1060_200.dba", dir=dirres)
