####
## Script to preform permutation test on the specified matrices of gene CN counts
####

library(devtools)
library(Hmisc)

#columns: BK_SAMN01057691	BK_SRR7813602	BB_SRR7758718	BB_SAMN02256313	BB_SAMN02256315	BB_SAMN02256316	BB_SAMN02256317	BB_SAMN02256318	BB_SAMN02256319	BB_SAMN02256320	BB_SAMN02256321	BB_SAMN02256322	PB_SAMN02261811	PB_SAMN02261819	PB_SAMN02261821	PB_SAMN02261826	PB_SAMN02261840	PB_SAMN02261845	PB_SAMN02261851	PB_SAMN02261853	PB_SAMN02261854	PB_SAMN02261856	PB_SAMN02261858	PB_SAMN02261865	PB_SAMN02261868	PB_SAMN02261870	PB_SAMN02261871	PB_SAMN02261878	PB_SAMN02261880

cnvFC = read.csv(file="freeC_rounded_cnv_matrix.tsv", sep = '\t', row.names = 1, header=TRUE)
cnvST = read.csv(file="samtools_rounded_cnv_matrix.tsv", sep = '\t', row.names = 1, header=TRUE)


###remove BB2018 from samples

BB2018 = as.list("BB_SRR7758718")

cnvFC = cnvFC[,!(names(cnvFC) %in% BB2018)]
cnvST = cnvST[,!(names(cnvST) %in% BB2018)]

#define sample lists

samp=c("BK_SAMN*91","BK_SRR*02","BB_SAMN*13","BB_SAMN*15","BB_SAMN*16","BB_SAMN*17","BB_SAMN*18","BB_SAMN*19","BB_SAMN*20","BB_SAMN*21","BB_SAMN*22","PB_SAMN*11","PB_SAMN*19","PB_SAMN*21","PB_SAMN*26","PB_SAMN*40","PB_SAMN*45","PB_SAMN*51","PB_SAMN*53","PB_SAMN*54","PB_SAMN*56","PB_SAMN*58","PB_SAMN*65","PB_SAMN*68","PB_SAMN*70","PB_SAMN*71","PB_SAMN*78","PB_SAMN*80")

spec=c("BkB","BkB","BrB","BrB","BrB","BrB","BrB","BrB","BrB","BrB","BrB","PB","PB","PB","PB","PB","PB","PB","PB","PB","PB","PB","PB","PB","PB","PB","PB","PB")

#spec=c("Black Bear","Black Bear","Brown Bear","Brown Bear","Brown Bear","Brown Bear","Brown Bear","Brown Bear","Brown Bear","Brown Bear","Brown Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear")


##General format for Vst calculation for each row
####
##Vp*n==[var(BB)*pBB+var(PB)*pPB]
##[var(BBPB) - (Vp*n/pPBBB)]/var(BBPB)
####
##[var(BBPB) - ([var(BB)*pBB+var(PB)*pPB]/pPBBB)]/var(BBPB)
###
##[var(BBPB) - ([var(BB)*9+var(PB)*17]/26)]/var(BBPB)
####

#calculate per-gene CN variance in BB, PB and BBPB
V1all = t(apply(cnvFC, 1, function(y) c(var(y[3:11])*8/9, var(y[12:28])*16/17, var(y[3:28])*25/26))) ##equal var.p in Excel
V2all = t(apply(cnvST, 1, function(y) c(var(y[3:11])*8/9, var(y[12:28])*16/17, var(y[3:28])*25/26))) ##equal var.p in Excel

#calculate per-gene Vst
VstFCall = (V1all[,3] - ((V1all[,1]*9 + V1all[,2]*17)/26))/V1all[,3]
VstSTall = (V2all[,3] - ((V2all[,1]*9 + V2all[,2]*17)/26))/V2all[,3]

#report correlation of Vst values between CN calling methods

cor(VstFCall,VstSTall, use="complete.obs", method = "pearson")

###
#Remove all genes with variance == 0 to create a reduced matric containing on CN variable genes
###

NonzeroVarFC = names(which(V1all[,3]>0))
NonzeroVarST = names(which(V2all[,3]>0))

cnvFCreduced = cnvFC[(rownames(cnvFC) %in% NonzeroVarFC),]
cnvSTreduced = cnvST[(rownames(cnvST) %in% NonzeroVarST),]
####
#calculate Vst (this is redundant from above but it only takes a few seconds)

V1 = t(apply(cnvFCreduced, 1, function(y) c(var(y[3:11])*8/9, var(y[12:28])*16/17, var(y[3:28])*25/26)))
V2 = t(apply(cnvSTreduced, 1, function(y) c(var(y[3:11])*8/9, var(y[12:28])*16/17, var(y[3:28])*25/26)))

VstFCreduced = (V1[,3] - ((V1[,1]*9 + V1[,2]*17)/26))/V1[,3]
VstSTreduced = (V2[,3] - ((V2[,1]*9 + V2[,2]*17)/26))/V2[,3]

#shuffle reduced CN matrix
cnvFCrdm = cnvFCreduced[, c(1:2, sample(3:ncol(cnvFCreduced)))]
cnvSTrdm = cnvSTreduced[, c(1:2, sample(3:ncol(cnvSTreduced)))]

#calculate Vst on shuffled matrix

V1rdm = t(apply(cnvFCrdm , 1, function(y) c(var(y[3:11])*8/9, var(y[12:28])*16/17, var(y[3:28])*25/26)))
V2rdm = t(apply(cnvSTrdm , 1, function(y) c(var(y[3:11])*8/9, var(y[12:28])*16/17, var(y[3:28])*25/26)))

VstFCrdm = (V1rdm[,3] - ((V1rdm[,1]*9 + V1rdm[,2]*17)/26))/V1rdm[,3]
VstSTrdm = (V2rdm[,3] - ((V2rdm[,1]*9 + V2rdm[,2]*17)/26))/V2rdm[,3]

VstPerFC = VstFCrdm
VstPerST = VstSTrdm

#repeat x1000
for (i in 1:999){

cnvFCrdm = cnvFCreduced[, c(1:2, sample(3:ncol(cnvFCreduced)))]
cnvSTrdm = cnvSTreduced[, c(1:2, sample(3:ncol(cnvSTreduced)))]

V1rdm = t(apply(cnvFCrdm , 1, function(y) c(var(y[3:11])*8/9, var(y[12:28])*16/17, var(y[3:28])*25/26)))
V2rdm = t(apply(cnvSTrdm , 1, function(y) c(var(y[3:11])*8/9, var(y[12:28])*16/17, var(y[3:28])*25/26)))

VstFCrdm = (V1rdm[,3] - ((V1rdm[,1]*9 + V1rdm[,2]*17)/26))/V1rdm[,3]
VstSTrdm = (V2rdm[,3] - ((V2rdm[,1]*9 + V2rdm[,2]*17)/26))/V2rdm[,3]

#add vector to permutation test matrix
VstPerFC = cbind(VstFCrdm, VstPerFC)
VstPerST = cbind(VstSTrdm, VstPerST)
}

############
##Extract cutoffs for each gene from permutations
############
VstFCcutoff.50 = apply(VstPerFC, 1, function(x) quantile(x, probs=.50))
VstFCcutoff.95 = apply(VstPerFC, 1, function(x) quantile(x, probs=.95))
VstFCcutoff.99 = apply(VstPerFC, 1, function(x) quantile(x, probs=.99))


VstSTcutoff.50 = apply(VstPerST, 1, function(x) quantile(x, probs=.50))
VstSTcutoff.95 = apply(VstPerST, 1, function(x) quantile(x, probs=.95))
VstSTcutoff.99 = apply(VstPerST, 1, function(x) quantile(x, probs=.99))

VstCutoffsFC = cbind(VstFCcutoff.50, VstFCcutoff.95, VstFCcutoff.99, VstFCreduced)
VstCutoffsST = cbind(VstSTcutoff.50, VstSTcutoff.95, VstSTcutoff.99, VstSTreduced)