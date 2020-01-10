#####################
#RNA-seq AST analysis for gonococcus
#####################

#Generate PCA of all data
library(DESeq2)
library(ggplot2)

dat <- read.csv("counts.csv", head=T)
rownames(dat) <- dat$Geneid
dat <- dat[,-1]

strain <- factor(c(rep("DAL.S", 9),rep("ALB", 9),rep("ORA", 9),rep("SEA", 9),rep("AR.0181", 9),rep("AR.0179", 9),rep("AR.0167", 9),rep("KCY", 9), rep("DTR", 9),rep("LAX", 9)))
treatment <- factor(rep(c(rep("control", 3), rep("T60", 3), rep("T180", 3)),10))
resistance <- factor(c(rep("sus", 36), rep("res", 54)))
BAPS <- factor(c(rep("8", 9),rep("4", 9),rep("7", 9),rep("6", 9),rep("4", 9),rep("7", 9),rep("7", 9),rep("4", 9), rep("8", 9),rep("6", 9)))

coldata <- data.frame("SampleName"=colnames(dat), strain, treatment, resistance, BAPS)
dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata, design = ~ resistance + treatment + BAPS + resistance:treatment)
dds <- DESeq(dds)

vsd <- varianceStabilizingTransformation(dds)
meanSdPlot(assay(vsd[notAllZero,]))
data <- plotPCA(vsd, intgroup = c( "strain", "treatment"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=BAPS, size=treatment, shape=resistance)) + geom_point() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_bw()

#####################
#Single time-point contrast to control (no drug condition)
#####################
  
dat <- read.csv("limited_counts.csv", head=T)
rownames(dat) <- dat$Geneid
dat <- dat[,-1]
strain <- factor(c(rep("DAL.S", 6),rep("ALB", 6),rep("ORA", 6),rep("SEA", 6),rep("AR.0181", 6),rep("AR.0179", 6),rep("AR.0167", 6),rep("KCY", 6), rep("DTR", 6),rep("LAX", 6))) 
treatment <- factor(rep(c(rep("control", 3), rep("T180", 3)),10))
resistance <- factor(c(rep("sus", 24), rep("res", 36)))

coldata <- data.frame("SampleName"=colnames(dat), strain, treatment, resistance)
dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata, design = ~ resistance + treatment + resistance:treatment)

dds$treatment <- relevel(dds$treatment, "T180")
dds$resistance <- relevel(dds$resistance, "sus")

dds <- DESeq(dds)

###################
# MAIN Phenotype effect in REFERENCE treatment (T180)
###################

resPhenoInT180 <- results(dds, contrast=c("resistance","sus","res"))
table(resPhenoInT180$padj < 0.01)
resPhenoInT180Sig <- subset(resPhenoInT180, padj < 0.01)
write.table(resPhenoInT180Sig, "~/Desktop/resPhenoInT180Sig_No28bl.txt", quote=F, sep="\t")

###################
# MAIN TREATMENT effect in REFERENCE strain (sus)
###################

resTreatmentsus <- results(dds, contrast=c("treatment","T180","control"))
table(resTreatmentsus$padj < 0.01)
resTreatmentsusSig <- subset(resTreatmentsus, padj < 0.01)
write.table(resTreatmentsusSig, "~/Desktop/resTreatmentsusSig_No28bl.txt", quote=F, sep="\t")

###################
# MAIN TREATMENT effect in Non-REFERENCE strain (res) (add the interaction to the main effect)
###################

resTreatmentres <- results(dds, contrast=list(c("treatment_control_vs_T180","resistanceres.treatmentcontrol")))
table(resTreatmentres$padj < 0.01)
resTreatmentresSig <- subset(resTreatmentres, padj < 0.01)
write.table(resTreatmentresSig, "~/Desktop/resTreatmentresSig_No28bl.txt", quote=F, sep="\t")

###################
# INTERACTION effect - is the treatment effect *different* across strains?
###################

resInteract <- results(dds, name="resistanceres.treatmentcontrol")
table(resInteract$padj < 0.01)
resInteractSig <- subset(resInteract, padj < 0.01)
write.table(resInteractSig, "~/Desktop/resInteractSig_No28bl.txt", quote=F, sep="\t")

###################
# MAIN Phenotype effect in Alternate treatment (control)
###################

resPhenoInControl <- results(dds, contrast=list(c("resistance_res_vs_sus","resistanceres.treatmentcontrol")))
table(resPhenoInControl$padj < 0.01)
resPhenoInControlSig <- subset(resPhenoInControl, padj < 0.01)
write.table(resPhenoInControlSig, "~/Desktop/resPhenoInControlSig_No28bl.txt", quote=F, sep="\t")

###################
#Comparing contrast to get expression patters
###################

exp_pats_I <- read.csv("~/Desktop/T180_trans/resInteractSig_No28bl.txt", head=T, sep="\t")
exp_pats_I$I <- c(rep("I",dim(exp_pats_I)[1]))

exp_pats_C <- read.csv("~/Desktop/T180_trans/resPhenoInControlSig_No28bl.txt", head=T, sep="\t")
exp_pats_C$C <- c(rep("C",dim(exp_pats_C)[1]))

exp_pats_A <- read.csv("~/Desktop/T180_trans/resPhenoInT180Sig_No28bl.txt", head=T, sep="\t")
exp_pats_A$A <- c(rep("A",dim(exp_pats_A)[1]))

exp_pats_R <- read.csv("~/Desktop/T180_trans/resTreatmentresSig_No28bl.txt", head=T, sep="\t")
exp_pats_R$R <- c(rep("R",dim(exp_pats_R)[1]))

exp_pats_S <- read.csv("~/Desktop/T180_trans/resTreatmentsusSig_No28bl.txt", head=T, sep="\t")
exp_pats_S$S <- c(rep("S",dim(exp_pats_S)[1]))

tot <- merge(exp_pats_I, exp_pats_C, by="geneID", all.x=T, all.y=T)
tot <- merge(tot, exp_pats_A, by="geneID", all.x=T, all.y=T)
tot <- merge(tot, exp_pats_R, by="geneID", all.x=T, all.y=T)
tot <- merge(tot, exp_pats_S, by="geneID", all.x=T, all.y=T)

write.csv(tot, "~/Desktop/tot_No28bl.csv")

###################################
#Plotting by gene
###################################


#Plot by gene
data3 <- plotCounts(dds, gene=c("CDS:NGO0049:-:hypothetical_protein:cds38:-:533:NC_002946.2"), intgroup=c("strain","treatment"), returnData=TRUE)
ggplot(data3, aes(x= treatment, y=count, color= resistance)) + geom_point(position=position_jitter(width=.1,height=0), size=2) + stat_summary(fun.y="mean", geom="line", aes(group=factor(strain)), size=1) +scale_x_discrete(limits=c("control","T180")) + theme_bw()

