# SNG transcriptomic analysis
library('stringr');library('DESeq2');library('ggplot2');library('readr');
library('grid');library('gridExtra');library('igraph'); library('limma')

# Set directories
MAINDIR  <- "/Users/tkay/Desktop/Work/SNG/data"
FIGDIR   <- "/Users/tkay/Desktop/Work/SNG/figures"

# Import data
setwd(MAINDIR)
cts  <- read.csv("brains_ABC.csv", row.names = 1)
cts_O <- read.csv("brains_O.csv", check.names = F, row.names = 1)
pairwise_A <- read.csv("EdgeData_A.csv")
pairwise_B <- read.csv("EdgeData_B.csv")
pairwise_C <- read.csv("EdgeData_C.csv")
pairwise_O <- read.csv("EdgeData_O.csv")

COLONIES <- c("A", "B", "C", "O")
for (col in COLONIES){
  assign(paste0("node_", col), read.csv(paste0("md_", col, "3.csv"), row.names=NULL))
}

# Separate transcriptomes by colony
cts_A <- cts[,1:90]
cts_B <- cts[,91:194]
cts_C <- cts[,195:299]

# Rename individuals for comparison across data-sets
newnames_A <- c()
for(id in 1:ncol(cts_A)){
  newnames_A <- c(newnames_A, as.numeric(str_split(colnames(cts_A)[id], pattern="_")[[1]][2]))
}
colnames(cts_A) <- newnames_A - 1
newnames_B <- c()
for(id in 1:ncol(cts_B)){
  newnames_B <- c(newnames_B, as.numeric(str_split(colnames(cts_B)[id], pattern="_")[[1]][2]))
}
colnames(cts_B) <- newnames_B - 1
newnames_C <- c()
for(id in 1:ncol(cts_C)){
  newnames_C <- c(newnames_C, as.numeric(str_split(colnames(cts_C)[id], pattern="_")[[1]][2]))
}
colnames(cts_C) <- newnames_C - 1
newnames_O <- c()
for(id in 1:ncol(cts_O)){
  newnames_O <- c(newnames_O, as.numeric(str_split(colnames(cts_O)[id], pattern="_")[[1]][2]))
}
colnames(cts_O) <- newnames_O

# Remove IDs which are not present in the individual meta-data - i.e., labeling errors
cts_A <- cts_A[,colnames(cts_A) %in% node_A$Tag] # removes 3
cts_B <- cts_B[,colnames(cts_B) %in% node_B$Tag] # removes 5
cts_C <- cts_C[,colnames(cts_C) %in% node_C$Tag] # removes 1

# Convert colnames from tagIDs to antIDs
for(i in 1:ncol(cts_A)){
  tag <- colnames(cts_A)[i]
  colnames(cts_A)[i] <- node_A$antID[which(node_A$Tag == tag)]
}
for(i in 1:ncol(cts_B)){
  tag <- colnames(cts_B)[i]
  colnames(cts_B)[i] <- node_B$antID[which(node_B$Tag == tag)]
}
for(i in 1:ncol(cts_C)){
  tag <- colnames(cts_C)[i]
  colnames(cts_C)[i] <- node_C$antID[which(node_C$Tag == tag)]
}

# Trim node data to match count data
node_A_t <- node_A[node_A$antID %in% colnames(cts_A),]
node_B_t <- node_B[node_B$antID %in% colnames(cts_B),]
node_C_t <- node_C[node_C$antID %in% colnames(cts_C),]
node_O_t <- node_O[node_O$antID %in% colnames(cts_O),]
node_O_t <- node_O_t[!is.na(node_O_t$Age),]

# reformat old brains data
cts_O_t <- cts_O[,order(as.numeric(colnames(cts_O)))]
cts_O_t <- cts_O_t[,colnames(cts_O_t) %in% node_O_t$antID]

# Name rows
row.names(node_A_t) <- node_A_t$antID
row.names(node_B_t) <- node_B_t$antID
row.names(node_C_t) <- node_C_t$antID
row.names(node_O_t) <- node_O_t$antID

# Remove one from A where there is no microbiota data
node_A_t <- na.omit(node_A_t)

# Normalise age between 0 and 1 (for DESeq2)
node_A_t$Age <- (node_A_t$Age-min(node_A_t$Age))/(max(node_A_t$Age)-min(node_A_t$Age))
node_B_t$Age <- (node_B_t$Age-min(node_B_t$Age))/(max(node_B_t$Age)-min(node_B_t$Age))
node_C_t$Age <- (node_C_t$Age-min(node_C_t$Age))/(max(node_C_t$Age)-min(node_C_t$Age))
node_O_t$Age <- (node_O_t$Age-min(node_O_t$Age))/(max(node_O_t$Age)-min(node_O_t$Age))

# Make extraction group non-numerical
node_A_t$Extraction <- paste0("g", node_A_t$Extraction)
node_B_t$Extraction <- paste0("g", node_B_t$Extraction)
node_C_t$Extraction <- paste0("g", node_C_t$Extraction)

# Reorder columns in A, where the IDs are not sequential
cts_A <- cts_A[,order(as.numeric(colnames(cts_A)))]
cts_A <- cts_A[,names(cts_A) %in% node_A_t$antID]

# Normalize old colony space use values
node_O_t$spacePCA <- (node_O_t$spacePCA - min(node_O_t$spacePCA)) / (min(node_O_t$spacePCA) - max(node_O_t$spacePCA))

ddsTask_A <- DESeqDataSetFromMatrix(countData = cts_A,
                                   colData = node_A_t,
                                   design = ~taskPCA)
ddsTask_B <- DESeqDataSetFromMatrix(countData = cts_B,
                                   colData = node_B_t,
                                   design = ~taskPCA)
ddsTask_C <- DESeqDataSetFromMatrix(countData = cts_C,
                                   colData = node_C_t,
                                   design = ~taskPCA)
ddsTask_O <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                   colData = node_O_t,
                                   design = ~taskPCA)
ddsSpace_A <- DESeqDataSetFromMatrix(countData = cts_A,
                                    colData = node_A_t,
                                    design = ~spacePCA)
ddsSpace_B <- DESeqDataSetFromMatrix(countData = cts_B,
                                    colData = node_B_t,
                                    design = ~spacePCA)
ddsSpace_C <- DESeqDataSetFromMatrix(countData = cts_C,
                                    colData = node_C_t,
                                    design = ~spacePCA)
ddsSpace_O <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                    colData = node_O_t,
                                    design = ~spacePCA)
ddsMat_A <- DESeqDataSetFromMatrix(countData = cts_A,
                                     colData = node_A_t,
                                     design = ~maturity)
ddsMat_B <- DESeqDataSetFromMatrix(countData = cts_B,
                                     colData = node_B_t,
                                     design = ~maturity)
ddsMat_C <- DESeqDataSetFromMatrix(countData = cts_C,
                                     colData = node_C_t,
                                     design = ~maturity)
ddsMat_O <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                   colData = node_O_t,
                                   design = ~maturity)
ddsAge_A <- DESeqDataSetFromMatrix(countData = cts_A,
                                   colData = node_A_t,
                                   design = ~Age)
ddsAge_B <- DESeqDataSetFromMatrix(countData = cts_B,
                                   colData = node_B_t,
                                   design = ~Age)
ddsAge_C <- DESeqDataSetFromMatrix(countData = cts_C,
                                   colData = node_C_t,
                                   design = ~Age)
ddsAge_O <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                   colData = node_O_t,
                                   design = ~Age)
ddsMic_A <- DESeqDataSetFromMatrix(countData = cts_A,
                                   colData = node_A_t,
                                   design = ~MicPCAm)
ddsMic_B <- DESeqDataSetFromMatrix(countData = cts_B,
                                   colData = node_B_t,
                                   design = ~MicPCAm)
ddsMic_C <- DESeqDataSetFromMatrix(countData = cts_C,
                                   colData = node_C_t,
                                   design = ~MicPCAm)
ddsMic_O <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                   colData = node_O_t,
                                   design = ~MicPCAm)

keep <- rowSums(counts(ddsAge_A)) + rowSums(counts(ddsAge_B)) + rowSums(counts(ddsAge_C)+ rowSums(counts(ddsAge_O)))  >= 100

# Remove low-count genes
ddsTask_A <- ddsTask_A[keep,]
ddsTask_B <- ddsTask_B[keep,]
ddsTask_C <- ddsTask_C[keep,] 
ddsTask_O <- ddsTask_O[keep,]
ddsSpace_A <- ddsSpace_A[keep,]
ddsSpace_B <- ddsSpace_B[keep,]
ddsSpace_C <- ddsSpace_C[keep,]
ddsSpace_O <- ddsSpace_O[keep,]
ddsMat_A <- ddsMat_A[keep,]
ddsMat_B <- ddsMat_B[keep,]
ddsMat_C <- ddsMat_C[keep,]
ddsMat_O <- ddsMat_O[keep,]
ddsAge_A <- ddsAge_A[keep,]
ddsAge_B <- ddsAge_B[keep,]
ddsAge_C <- ddsAge_C[keep,]
ddsAge_O <- ddsAge_O[keep,]
ddsMic_A <- ddsMic_A[keep,]
ddsMic_B <- ddsMic_B[keep,]
ddsMic_C <- ddsMic_C[keep,]
ddsMic_O <- ddsMic_O[keep,]

# Run DEA
ddsTask_A <- DESeq(ddsTask_A)
ddsTask_B <- DESeq(ddsTask_B)
ddsTask_C <- DESeq(ddsTask_C)
ddsTask_O <- DESeq(ddsTask_O)
ddsSpace_A <- DESeq(ddsSpace_A)
ddsSpace_B <- DESeq(ddsSpace_B)
ddsSpace_C <- DESeq(ddsSpace_C)
ddsSpace_O <- DESeq(ddsSpace_O)
ddsMat_A <- DESeq(ddsMat_A)
ddsMat_B <- DESeq(ddsMat_B)
ddsMat_C <- DESeq(ddsMat_C)
ddsMat_O <- DESeq(ddsMat_O)
ddsAge_A <- DESeq(ddsAge_A)
ddsAge_B <- DESeq(ddsAge_B)
ddsAge_C <- DESeq(ddsAge_C)
ddsAge_O <- DESeq(ddsAge_O)
ddsMic_A <- DESeq(ddsMic_A)
ddsMic_B <- DESeq(ddsMic_B)
ddsMic_C <- DESeq(ddsMic_C)
ddsMic_O <- DESeq(ddsMic_O)

# Print results
summary(results(ddsMat_A))
summary(results(ddsMat_B))
summary(results(ddsMat_C))
summary(results(ddsMat_O))

taskPCA_geneIDs <- unique(c(row.names(na.omit(as.data.frame(results(ddsTask_A))[as.data.frame(results(ddsTask_A))$padj < 0.05,])),
                     row.names(na.omit(as.data.frame(results(ddsTask_B))[as.data.frame(results(ddsTask_B))$padj < 0.05,])),
                     row.names(na.omit(as.data.frame(results(ddsTask_C))[as.data.frame(results(ddsTask_C))$padj < 0.05,])),
                     row.names(na.omit(as.data.frame(results(ddsTask_O))[as.data.frame(results(ddsTask_O))$padj < 0.05,]))))
spacePCA_geneIDs <- unique(c(row.names(na.omit(as.data.frame(results(ddsSpace_A))[as.data.frame(results(ddsSpace_A))$padj < 0.05,])),
                            row.names(na.omit(as.data.frame(results(ddsSpace_B))[as.data.frame(results(ddsSpace_B))$padj < 0.05,])),
                            row.names(na.omit(as.data.frame(results(ddsSpace_C))[as.data.frame(results(ddsSpace_C))$padj < 0.05,])),
                            row.names(na.omit(as.data.frame(results(ddsSpace_O))[as.data.frame(results(ddsSpace_O))$padj < 0.05,]))))
maturity_geneIDs <- unique(c(row.names(na.omit(as.data.frame(results(ddsMat_A))[as.data.frame(results(ddsMat_A))$padj < 0.05,])),
                             row.names(na.omit(as.data.frame(results(ddsMat_B))[as.data.frame(results(ddsMat_B))$padj < 0.05,])),
                             row.names(na.omit(as.data.frame(results(ddsMat_C))[as.data.frame(results(ddsMat_C))$padj < 0.05,])),
                             row.names(na.omit(as.data.frame(results(ddsMat_O))[as.data.frame(results(ddsMat_O))$padj < 0.05,]))))
age_geneIDs <- unique(c(row.names(na.omit(as.data.frame(results(ddsAge_A))[as.data.frame(results(ddsAge_A))$padj < 0.05,])),
                             row.names(na.omit(as.data.frame(results(ddsAge_B))[as.data.frame(results(ddsAge_B))$padj < 0.05,])),
                             row.names(na.omit(as.data.frame(results(ddsAge_C))[as.data.frame(results(ddsAge_C))$padj < 0.05,])),
                             row.names(na.omit(as.data.frame(results(ddsAge_O))[as.data.frame(results(ddsAge_O))$padj < 0.05,]))))
mice_geneIDs <- unique(c(row.names(na.omit(as.data.frame(results(ddsMic_A))[as.data.frame(results(ddsMic_A))$padj < 0.05,])),
                        row.names(na.omit(as.data.frame(results(ddsMic_B))[as.data.frame(results(ddsMic_B))$padj < 0.05,])),
                        row.names(na.omit(as.data.frame(results(ddsMic_C))[as.data.frame(results(ddsMic_C))$padj < 0.05,])),
                        row.names(na.omit(as.data.frame(results(ddsMic_O))[as.data.frame(results(ddsMic_O))$padj < 0.05,]))))


genesDE <- unique(c(taskPCA_geneIDs, spacePCA_geneIDs, maturity_geneIDs, age_geneIDs, mice_geneIDs))

cts_At <- cts_A[rownames(cts_A) %in% genesDE, ]
for (c in 1:ncol(cts_At)){
  cts_At[,c] <- cts_At[,c]  / sum(cts_At[,1])
}
cts_Bt <- cts_B[rownames(cts_B) %in% genesDE, ]
for (c in 1:ncol(cts_Bt)){
  cts_Bt[,c] <- cts_Bt[,c]  / sum(cts_Bt[,1])
}
cts_Ct <- cts_C[rownames(cts_C) %in% genesDE, ]
for (c in 1:ncol(cts_Ct)){
  cts_Ct[,c] <- cts_Ct[,c]  / sum(cts_Ct[,1])
}
cts_Ot <- cts_O[rownames(cts_O) %in% genesDE, ]
for (c in 1:ncol(cts_Ot)){
  cts_Ot[,c] <- cts_Ot[,c]  / sum(cts_Ot[,1])
}

# vst transformation
vsd_A <- vst(ddsAge_A, blind = FALSE)
vsd_B <- vst(ddsAge_B, blind = FALSE)
vsd_C <- vst(ddsAge_C, blind = FALSE)
vsd_O <- vst(ddsAge_O, blind = FALSE)

# Create PCA
pcaData_A <- plotPCA(vsd_A, intgroup=c("maturity"), returnData=TRUE)
pcaData_B <- plotPCA(vsd_B, intgroup=c("maturity"), returnData=TRUE)
pcaData_C <- plotPCA(vsd_C, intgroup=c("maturity"), returnData=TRUE)
pcaData_O <- plotPCA(vsd_O, intgroup=c("maturity"), returnData=TRUE)

# How much variance along PCs 1&2
percentVar_A <- round(100*attr(pcaData_A, "percentVar"))
percentVar_B <- round(100*attr(pcaData_B, "percentVar"))
percentVar_C <- round(100*attr(pcaData_C, "percentVar"))
percentVar_O <- round(100*attr(pcaData_O, "percentVar"))

# Plot PCAs
PCA_A <- ggplot(pcaData_A, aes(PC1, PC2, color=maturity)) +
  geom_point(size=50) +
  xlab(paste0(percentVar_A[1], "% variance")) +
  ylab(paste0(percentVar_A[2], "% variance")) +
  coord_fixed() +
  theme(axis.text=element_text(size=0), axis.title=element_text(size=100,face="bold"),
        legend.text=element_text(size=0), legend.title=element_text(size=0))
PCA_B <- ggplot(pcaData_B, aes(PC1, PC2, color=maturity)) +
  geom_point(size=50) +
  xlab(paste0(percentVar_B[1], "% variance")) +
  ylab(paste0(percentVar_B[2], "% variance")) +
  coord_fixed() +
  theme(axis.text=element_text(size=0), axis.title=element_text(size=100,face="bold"),
        legend.text=element_text(size=0), legend.title=element_text(size=0))
PCA_C <- ggplot(pcaData_C, aes(PC1, PC2, color=maturity)) +
  geom_point(size=50) +
  xlab(paste0(percentVar_C[1], "% variance")) +
  ylab(paste0(percentVar_C[2], "% variance")) +
  coord_fixed() +
  theme(axis.text=element_text(size=0), axis.title=element_text(size=100,face="bold"),
        legend.text=element_text(size=0), legend.title=element_text(size=0))
PCA_O <- ggplot(pcaData_O, aes(PC1, PC2, color=maturity)) +
  geom_point(size=50) +
  xlab(paste0(percentVar_C[1], "% variance")) +
  ylab(paste0(percentVar_C[2], "% variance")) +
  coord_fixed() +
  theme(axis.text=element_text(size=0), axis.title=element_text(size=100,face="bold"),
        legend.text=element_text(size=0), legend.title=element_text(size=0))

setwd(FIGDIR)
jpeg('PCA_Brains_ABCO.jpg', width=6000, height=2000, unit='px')
grid.arrange(PCA_A,
             PCA_B,
             PCA_C,
             PCA_O,
             ncol=2)
dev.off()

c(summary(lm(pcaData_A$PC2 ~ node_A_t$maturity))$adj.r.squared,
summary(lm(pcaData_B$PC2 ~ node_B_t$maturity))$adj.r.squared,
summary(lm(pcaData_C$PC2 ~ node_C_t$maturity))$adj.r.squared,
summary(lm(pcaData_O$PC2 ~ node_O_t$maturity))$adj.r.squared)

node_A_t$genePC1 <- pcaData_A$PC1
node_B_t$genePC1 <- pcaData_B$PC1
node_C_t$genePC1 <- pcaData_C$PC1
node_O_t$genePC1 <- pcaData_O$PC1
node_A_t$genePC2 <- pcaData_A$PC2
node_B_t$genePC2 <- pcaData_B$PC2
node_C_t$genePC2 <- pcaData_C$PC2
node_O_t$genePC2 <- pcaData_O$PC2

setwd(MAINDIR)
write.csv(node_A_t, "md_A4.csv", row.names=FALSE)
write.csv(node_B_t, "md_B4.csv", row.names=FALSE)
write.csv(node_C_t, "md_C4.csv", row.names=FALSE)
write.csv(node_O_t, "md_O4.csv", row.names=FALSE)

# Repeat DEAs while controlling for different combinations of variables

# Social and behavioral
dds_mtA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ maturity + taskPCA)
dds_tmA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ taskPCA + maturity)
dds_mtA <- dds_mtA[keep,]
dds_tmA <- dds_tmA[keep,]
dds_mtA <- DESeq(dds_mtA)
dds_tmA <- DESeq(dds_tmA)
summary(results(dds_mtA))
summary(results(dds_tmA))

dds_mtB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ maturity + taskPCA)
dds_tmB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ taskPCA + maturity)
dds_mtB <- dds_mtB[keep,]
dds_tmB <- dds_tmB[keep,]
dds_mtB <- DESeq(dds_mtB)
dds_tmB <- DESeq(dds_tmB)
summary(results(dds_mtB))
summary(results(dds_tmB))

dds_mtC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ maturity + taskPCA)
dds_tmC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ taskPCA + maturity)
dds_mtC <- dds_mtC[keep,]
dds_tmC <- dds_tmC[keep,]
dds_mtC <- DESeq(dds_mtC)
dds_tmC <- DESeq(dds_tmC)
summary(results(dds_mtC))
summary(results(dds_tmC))


dds_mtO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ maturity + taskPCA)
dds_tmO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ taskPCA + maturity)
dds_mtO <- dds_mtO[keep,]
dds_tmO <- dds_tmO[keep,]
dds_mtO <- DESeq(dds_mtO)
dds_tmO <- DESeq(dds_tmO)
summary(results(dds_mtO))
summary(results(dds_tmO))

# Social and age
dds_maA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ maturity + Age)
dds_amA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ Age + maturity)
dds_maA <- dds_maA[keep,]
dds_amA <- dds_amA[keep,]
dds_maA <- DESeq(dds_maA)
dds_amA <- DESeq(dds_amA)
summary(results(dds_maA))
summary(results(dds_amA))

dds_maB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ maturity + Age)
dds_amB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ Age + maturity)
dds_maB <- dds_maB[keep,]
dds_amB <- dds_amB[keep,]
dds_maB <- DESeq(dds_maB)
dds_amB <- DESeq(dds_amB)
summary(results(dds_maB))
summary(results(dds_amB))

dds_maC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ maturity + Age)
dds_amC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ Age + maturity)
dds_maC <- dds_maC[keep,]
dds_amC <- dds_amC[keep,]
dds_maC <- DESeq(dds_maC)
dds_amC <- DESeq(dds_amC)
summary(results(dds_maC))
summary(results(dds_amC))


dds_maO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ maturity + Age)
dds_amO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ Age + maturity)
dds_maO <- dds_maO[keep,]
dds_amO <- dds_amO[keep,]
dds_maO <- DESeq(dds_maO)
dds_amO <- DESeq(dds_amO)
summary(results(dds_maO))
summary(results(dds_amO))

# Maturity & space

dds_mpA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ maturity + spacePCA)
dds_pmA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ spacePCA + maturity)
dds_mpA <- dds_mpA[keep,]
dds_pmA <- dds_pmA[keep,]
dds_mpA <- DESeq(dds_mpA)
dds_pmA <- DESeq(dds_pmA)
summary(results(dds_mpA))
summary(results(dds_pmA))

dds_mpB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ maturity + spacePCA)
dds_pmB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ spacePCA + maturity)
dds_mpB <- dds_mpB[keep,]
dds_pmB <- dds_pmB[keep,]
dds_mpB <- DESeq(dds_mpB)
dds_pmB <- DESeq(dds_pmB)
summary(results(dds_mpB))
summary(results(dds_pmB))

dds_mpC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ maturity + spacePCA)
dds_pmC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ spacePCA + maturity)
dds_mpC <- dds_mpC[keep,]
dds_pmC <- dds_pmC[keep,]
dds_mpC <- DESeq(dds_mpC)
dds_pmC <- DESeq(dds_pmC)
summary(results(dds_mpC))
summary(results(dds_pmC))


dds_mpO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ maturity + spacePCA)
dds_pmO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ spacePCA + maturity)
dds_mpO <- dds_mpO[keep,]
dds_pmO <- dds_pmO[keep,]
dds_mpO <- DESeq(dds_mpO)
dds_pmO <- DESeq(dds_pmO)
summary(results(dds_mpO))
summary(results(dds_pmO))

# Maturity and microbiota

dds_mgA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ maturity + MicPCAm)
dds_gmA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ MicPCAm + maturity)
dds_mgA <- dds_mgA[keep,]
dds_gmA <- dds_gmA[keep,]
dds_mgA <- DESeq(dds_mgA)
dds_gmA <- DESeq(dds_gmA)
summary(results(dds_mgA))
summary(results(dds_gmA))

dds_mgB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ maturity + MicPCAm)
dds_gmB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ MicPCAm + maturity)
dds_mgB <- dds_mgB[keep,]
dds_gmB <- dds_gmB[keep,]
dds_mgB <- DESeq(dds_mgB)
dds_gmB <- DESeq(dds_gmB)
summary(results(dds_mgB))
summary(results(dds_gmB))

dds_mgC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ maturity + MicPCAm)
dds_gmC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ MicPCAm + maturity)
dds_mgC <- dds_mgC[keep,]
dds_gmC <- dds_gmC[keep,]
dds_mgC <- DESeq(dds_mgC)
dds_gmC <- DESeq(dds_gmC)
summary(results(dds_mgC))
summary(results(dds_gmC))


dds_mgO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ maturity + MicPCAm)
dds_gmO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ MicPCAm + maturity)
dds_mgO <- dds_mgO[keep,]
dds_gmO <- dds_gmO[keep,]
dds_mgO <- DESeq(dds_mgO)
dds_gmO <- DESeq(dds_gmO)
summary(results(dds_mgO))
summary(results(dds_gmO))

# Behabior v microbiota

dds_tgA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ taskPCA + MicPCAm)
dds_gtA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ MicPCAm + taskPCA)
dds_tgA <- dds_tgA[keep,]
dds_gtA <- dds_gtA[keep,]
dds_tgA <- DESeq(dds_tgA)
dds_gtA <- DESeq(dds_gtA)
summary(results(dds_tgA))
summary(results(dds_gtA))

dds_tgB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ taskPCA + MicPCAm)
dds_gtB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ MicPCAm + taskPCA)
dds_tgB <- dds_tgB[keep,]
dds_gtB <- dds_gtB[keep,]
dds_tgB <- DESeq(dds_tgB)
dds_gtB <- DESeq(dds_gtB)
summary(results(dds_tgB))
summary(results(dds_gtB))

dds_tgC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ taskPCA + MicPCAm)
dds_gtC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ MicPCAm + taskPCA)
dds_tgC <- dds_tgC[keep,]
dds_gtC <- dds_gtC[keep,]
dds_tgC <- DESeq(dds_tgC)
dds_gtC <- DESeq(dds_gtC)
summary(results(dds_tgC))
summary(results(dds_gtC))


dds_tgO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ taskPCA + MicPCAm)
dds_gtO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ MicPCAm + taskPCA)
dds_tgO <- dds_tgO[keep,]
dds_gtO <- dds_gtO[keep,]
dds_tgO <- DESeq(dds_tgO)
dds_gtO <- DESeq(dds_gtO)
summary(results(dds_tgO))
summary(results(dds_gtO))

# Behavior and age

dds_taA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ taskPCA + Age)
dds_atA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ Age + taskPCA)
dds_taA <- dds_taA[keep,]
dds_atA <- dds_atA[keep,]
dds_taA <- DESeq(dds_taA)
dds_atA <- DESeq(dds_atA)
summary(results(dds_taA))
summary(results(dds_atA))

dds_taB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ taskPCA + Age)
dds_atB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ Age + taskPCA)
dds_taB <- dds_taB[keep,]
dds_atB <- dds_atB[keep,]
dds_taB <- DESeq(dds_taB)
dds_atB <- DESeq(dds_atB)
summary(results(dds_taB))
summary(results(dds_atB))

dds_taC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ taskPCA + Age)
dds_atC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ Age + taskPCA)
dds_taC <- dds_taC[keep,]
dds_atC <- dds_atC[keep,]
dds_taC <- DESeq(dds_taC)
dds_atC <- DESeq(dds_atC)
summary(results(dds_taC))
summary(results(dds_atC))


dds_taO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ taskPCA + Age)
dds_atO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ Age + taskPCA)
dds_taO <- dds_taO[keep,]
dds_atO <- dds_atO[keep,]
dds_taO <- DESeq(dds_taO)
dds_atO <- DESeq(dds_atO)
summary(results(dds_taO))
summary(results(dds_atO))

# Behavior and physical

dds_tpA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ taskPCA + spacePCA)
dds_ptA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ spacePCA + taskPCA)
dds_tpA <- dds_tpA[keep,]
dds_ptA <- dds_ptA[keep,]
dds_tpA <- DESeq(dds_tpA)
dds_ptA <- DESeq(dds_ptA)
summary(results(dds_tpA))
summary(results(dds_ptA))

dds_tpB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ taskPCA + spacePCA)
dds_ptB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ spacePCA + taskPCA)
dds_tpB <- dds_tpB[keep,]
dds_ptB <- dds_ptB[keep,]
dds_tpB <- DESeq(dds_tpB)
dds_ptB <- DESeq(dds_ptB)
summary(results(dds_tpB))
summary(results(dds_ptB))

dds_tpC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ taskPCA + spacePCA)
dds_ptC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ spacePCA + taskPCA)
dds_tpC <- dds_tpC[keep,]
dds_ptC <- dds_ptC[keep,]
dds_tpC <- DESeq(dds_tpC)
dds_ptC <- DESeq(dds_ptC)
summary(results(dds_tpC))
summary(results(dds_ptC))


dds_tpO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ taskPCA + spacePCA)
dds_ptO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ spacePCA + taskPCA)
dds_tpO <- dds_tpO[keep,]
dds_ptO <- dds_ptO[keep,]
dds_tpO <- DESeq(dds_tpO)
dds_ptO <- DESeq(dds_ptO)
summary(results(dds_tpO))
summary(results(dds_ptO))

# Space and age

dds_apA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ Age + spacePCA)
dds_paA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ spacePCA + Age)
dds_apA <- dds_apA[keep,]
dds_paA <- dds_paA[keep,]
dds_apA <- DESeq(dds_apA)
dds_paA <- DESeq(dds_paA)
summary(results(dds_apA))
summary(results(dds_paA))

dds_apB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ Age + spacePCA)
dds_paB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ spacePCA + Age)
dds_apB <- dds_apB[keep,]
dds_paB <- dds_paB[keep,]
dds_apB <- DESeq(dds_apB)
dds_paB <- DESeq(dds_paB)
summary(results(dds_apB))
summary(results(dds_paB))

dds_apC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ Age + spacePCA)
dds_paC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ spacePCA + Age)
dds_apC <- dds_apC[keep,]
dds_paC <- dds_paC[keep,]
dds_apC <- DESeq(dds_apC)
dds_paC <- DESeq(dds_paC)
summary(results(dds_apC))
summary(results(dds_paC))


dds_apO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ Age + spacePCA)
dds_paO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ spacePCA + Age)
dds_apO <- dds_apO[keep,]
dds_paO <- dds_paO[keep,]
dds_apO <- DESeq(dds_apO)
dds_paO <- DESeq(dds_paO)
summary(results(dds_apO))
summary(results(dds_paO))

# Age and microbiota

dds_agA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ Age + MicPCAm)
dds_gaA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ MicPCAm + Age)
dds_agA <- dds_agA[keep,]
dds_gaA <- dds_gaA[keep,]
dds_agA <- DESeq(dds_agA)
dds_gaA <- DESeq(dds_gaA)
summary(results(dds_agA))
summary(results(dds_gaA))

dds_agB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ Age + MicPCAm)
dds_gaB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ MicPCAm + Age)
dds_agB <- dds_agB[keep,]
dds_gaB <- dds_gaB[keep,]
dds_agB <- DESeq(dds_agB)
dds_gaB <- DESeq(dds_gaB)
summary(results(dds_agB))
summary(results(dds_gaB))

dds_agC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ Age + MicPCAm)
dds_gaC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ MicPCAm + Age)
dds_agC <- dds_agC[keep,]
dds_gaC <- dds_gaC[keep,]
dds_agC <- DESeq(dds_agC)
dds_gaC <- DESeq(dds_gaC)
summary(results(dds_agC))
summary(results(dds_gaC))


dds_agO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ Age + MicPCAm)
dds_gaO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ MicPCAm + Age)
dds_agO <- dds_agO[keep,]
dds_gaO <- dds_gaO[keep,]
dds_agO <- DESeq(dds_agO)
dds_gaO <- DESeq(dds_gaO)
summary(results(dds_agO))
summary(results(dds_gaO))

# Physical and microbiota

dds_pgA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ spacePCA + MicPCAm)
dds_gpA <- DESeqDataSetFromMatrix(countData = cts_A,
                                  colData = node_A_t,
                                  design = ~ MicPCAm + spacePCA)
dds_pgA <- dds_pgA[keep,]
dds_gpA <- dds_gpA[keep,]
dds_pgA <- DESeq(dds_pgA)
dds_gpA <- DESeq(dds_gpA)
summary(results(dds_pgA))
summary(results(dds_gpA))

dds_pgB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ spacePCA + MicPCAm)
dds_gpB <- DESeqDataSetFromMatrix(countData = cts_B,
                                  colData = node_B_t,
                                  design = ~ MicPCAm + spacePCA)
dds_pgB <- dds_pgB[keep,]
dds_gpB <- dds_gpB[keep,]
dds_pgB <- DESeq(dds_pgB)
dds_gpB <- DESeq(dds_gpB)
summary(results(dds_pgB))
summary(results(dds_gpB))

dds_pgC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ spacePCA + MicPCAm)
dds_gpC <- DESeqDataSetFromMatrix(countData = cts_C,
                                  colData = node_C_t,
                                  design = ~ MicPCAm + spacePCA)
dds_pgC <- dds_pgC[keep,]
dds_gpC <- dds_gpC[keep,]
dds_pgC <- DESeq(dds_pgC)
dds_gpC <- DESeq(dds_gpC)
summary(results(dds_pgC))
summary(results(dds_gpC))


dds_pgO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ spacePCA + MicPCAm)
dds_gpO <- DESeqDataSetFromMatrix(countData = cts_O_t,
                                  colData = node_O_t,
                                  design = ~ MicPCAm + spacePCA)
dds_pgO <- dds_pgO[keep,]
dds_gpO <- dds_gpO[keep,]
dds_pgO <- DESeq(dds_pgO)
dds_gpO <- DESeq(dds_gpO)
summary(results(dds_pgO))
summary(results(dds_gpO))



sum(na.omit(results(dds_gpO)$padj) < 0.1)

SummaryDEGs <- as.data.frame(matrix(NA, nrow = 25, ncol = 4))
colnames(SummaryDEGs) <- c("A", "B", "C", "O")
rownames(SummaryDEGs) <- c("Social",
                           "Social controlling for behavior",
                           "Social controlling for age",
                           "Social controlling for physical",
                           "Social controlling for microbiota",
                           "Behavior",
                           "Behavior controlling for social",
                           "Behavior controlling for age",
                           "Behavior controlling for physical",
                           "Behavior controlling for microbiota",
                           "Age",
                           "Age controlling for social",
                           "Age controlling for behavior",
                           "Age controlling for physical",
                           "Age controlling for microbiota",
                           "Physical",
                           "Physical controlling for social",
                           "Physical controlling for behavior",
                           "Physical controlling for age",
                           "Physical controlling for microbiota",
                           "Microbiota",
                           "Microbiota controlling for social",
                           "Microbiota controlling for behavior",
                           "Microbiota controlling for age",
                           "Microbiota controlling for physical")
SummaryDEGs[1,1] <- sum(na.omit(results(ddsMat_A)$padj) < 0.05)
SummaryDEGs[1,2] <- sum(na.omit(results(ddsMat_B)$padj) < 0.05)
SummaryDEGs[1,3] <- sum(na.omit(results(ddsMat_C)$padj) < 0.05)
SummaryDEGs[1,4] <- sum(na.omit(results(ddsMat_O)$padj) < 0.05)
SummaryDEGs[2,1] <- sum(na.omit(results(dds_tmA)$padj) < 0.05)
SummaryDEGs[2,2] <- sum(na.omit(results(dds_tmB)$padj) < 0.05)
SummaryDEGs[2,3] <- sum(na.omit(results(dds_tmC)$padj) < 0.05)
SummaryDEGs[2,4] <- sum(na.omit(results(dds_tmO)$padj) < 0.05)
SummaryDEGs[3,1] <- sum(na.omit(results(dds_amA)$padj) < 0.05)
SummaryDEGs[3,2] <- sum(na.omit(results(dds_amB)$padj) < 0.05)
SummaryDEGs[3,3] <- sum(na.omit(results(dds_amC)$padj) < 0.05)
SummaryDEGs[3,4] <- sum(na.omit(results(dds_amO)$padj) < 0.05)
SummaryDEGs[4,1] <- sum(na.omit(results(dds_pmA)$padj) < 0.05)
SummaryDEGs[4,2] <- sum(na.omit(results(dds_pmB)$padj) < 0.05)
SummaryDEGs[4,3] <- sum(na.omit(results(dds_pmC)$padj) < 0.05)
SummaryDEGs[4,4] <- sum(na.omit(results(dds_pmO)$padj) < 0.05)
SummaryDEGs[5,1] <- sum(na.omit(results(dds_gmA)$padj) < 0.05)
SummaryDEGs[5,2] <- sum(na.omit(results(dds_gmB)$padj) < 0.05)
SummaryDEGs[5,3] <- sum(na.omit(results(dds_gmC)$padj) < 0.05)
SummaryDEGs[5,4] <- sum(na.omit(results(dds_gmO)$padj) < 0.05)
SummaryDEGs[6,1] <- sum(na.omit(results(ddsTask_A)$padj) < 0.05)
SummaryDEGs[6,2] <- sum(na.omit(results(ddsTask_B)$padj) < 0.05)
SummaryDEGs[6,3] <- sum(na.omit(results(ddsTask_C)$padj) < 0.05)
SummaryDEGs[6,4] <- sum(na.omit(results(ddsTask_O)$padj) < 0.05)
SummaryDEGs[7,1] <- sum(na.omit(results(dds_mtA)$padj) < 0.05)
SummaryDEGs[7,2] <- sum(na.omit(results(dds_mtB)$padj) < 0.05)
SummaryDEGs[7,3] <- sum(na.omit(results(dds_mtC)$padj) < 0.05)
SummaryDEGs[7,4] <- sum(na.omit(results(dds_mtO)$padj) < 0.05)
SummaryDEGs[8,1] <- sum(na.omit(results(dds_atA)$padj) < 0.05)
SummaryDEGs[8,2] <- sum(na.omit(results(dds_atB)$padj) < 0.05)
SummaryDEGs[8,3] <- sum(na.omit(results(dds_atC)$padj) < 0.05)
SummaryDEGs[8,4] <- sum(na.omit(results(dds_atO)$padj) < 0.05)
SummaryDEGs[9,1] <- sum(na.omit(results(dds_ptA)$padj) < 0.05)
SummaryDEGs[9,2] <- sum(na.omit(results(dds_ptB)$padj) < 0.05)
SummaryDEGs[9,3] <- sum(na.omit(results(dds_ptC)$padj) < 0.05)
SummaryDEGs[9,4] <- sum(na.omit(results(dds_ptO)$padj) < 0.05)
SummaryDEGs[10,1] <- sum(na.omit(results(dds_gtA)$padj) < 0.05)
SummaryDEGs[10,2] <- sum(na.omit(results(dds_gtB)$padj) < 0.05)
SummaryDEGs[10,3] <- sum(na.omit(results(dds_gtC)$padj) < 0.05)
SummaryDEGs[10,4] <- sum(na.omit(results(dds_gtO)$padj) < 0.05)
SummaryDEGs[11,1] <- sum(na.omit(results(ddsAge_A)$padj) < 0.05)
SummaryDEGs[11,2] <- sum(na.omit(results(ddsAge_B)$padj) < 0.05)
SummaryDEGs[11,3] <- sum(na.omit(results(ddsAge_C)$padj) < 0.05)
SummaryDEGs[11,4] <- sum(na.omit(results(ddsAge_O)$padj) < 0.05)
SummaryDEGs[12,1] <- sum(na.omit(results(dds_maA)$padj) < 0.05)
SummaryDEGs[12,2] <- sum(na.omit(results(dds_maB)$padj) < 0.05)
SummaryDEGs[12,3] <- sum(na.omit(results(dds_maC)$padj) < 0.05)
SummaryDEGs[12,4] <- sum(na.omit(results(dds_maO)$padj) < 0.05)
SummaryDEGs[13,1] <- sum(na.omit(results(dds_taA)$padj) < 0.05)
SummaryDEGs[13,2] <- sum(na.omit(results(dds_taB)$padj) < 0.05)
SummaryDEGs[13,3] <- sum(na.omit(results(dds_taC)$padj) < 0.05)
SummaryDEGs[13,4] <- sum(na.omit(results(dds_taO)$padj) < 0.05)
SummaryDEGs[14,1] <- sum(na.omit(results(dds_paA)$padj) < 0.05)
SummaryDEGs[14,2] <- sum(na.omit(results(dds_paB)$padj) < 0.05)
SummaryDEGs[14,3] <- sum(na.omit(results(dds_paC)$padj) < 0.05)
SummaryDEGs[14,4] <- sum(na.omit(results(dds_paO)$padj) < 0.05)
SummaryDEGs[15,1] <- sum(na.omit(results(dds_gaA)$padj) < 0.05)
SummaryDEGs[15,2] <- sum(na.omit(results(dds_gaB)$padj) < 0.05)
SummaryDEGs[15,3] <- sum(na.omit(results(dds_gaC)$padj) < 0.05)
SummaryDEGs[15,4] <- sum(na.omit(results(dds_gaO)$padj) < 0.05)
SummaryDEGs[16,1] <- sum(na.omit(results(ddsSpace_A)$padj) < 0.05)
SummaryDEGs[16,2] <- sum(na.omit(results(ddsSpace_B)$padj) < 0.05)
SummaryDEGs[16,3] <- sum(na.omit(results(ddsSpace_C)$padj) < 0.05)
SummaryDEGs[16,4] <- sum(na.omit(results(ddsSpace_O)$padj) < 0.05)
SummaryDEGs[17,1] <- sum(na.omit(results(dds_mpA)$padj) < 0.05)
SummaryDEGs[17,2] <- sum(na.omit(results(dds_mpB)$padj) < 0.05)
SummaryDEGs[17,3] <- sum(na.omit(results(dds_mpC)$padj) < 0.05)
SummaryDEGs[17,4] <- sum(na.omit(results(dds_mpO)$padj) < 0.05)
SummaryDEGs[18,1] <- sum(na.omit(results(dds_tpA)$padj) < 0.05)
SummaryDEGs[18,2] <- sum(na.omit(results(dds_tpB)$padj) < 0.05)
SummaryDEGs[18,3] <- sum(na.omit(results(dds_tpC)$padj) < 0.05)
SummaryDEGs[18,4] <- sum(na.omit(results(dds_tpO)$padj) < 0.05)
SummaryDEGs[19,1] <- sum(na.omit(results(dds_apA)$padj) < 0.05)
SummaryDEGs[19,2] <- sum(na.omit(results(dds_apB)$padj) < 0.05)
SummaryDEGs[19,3] <- sum(na.omit(results(dds_apC)$padj) < 0.05)
SummaryDEGs[19,4] <- sum(na.omit(results(dds_apO)$padj) < 0.05)
SummaryDEGs[20,1] <- sum(na.omit(results(dds_gpA)$padj) < 0.05)
SummaryDEGs[20,2] <- sum(na.omit(results(dds_gpB)$padj) < 0.05)
SummaryDEGs[20,3] <- sum(na.omit(results(dds_gpC)$padj) < 0.05)
SummaryDEGs[20,4] <- sum(na.omit(results(dds_gpO)$padj) < 0.05)
SummaryDEGs[21,1] <- sum(na.omit(results(ddsMic_A)$padj) < 0.05)
SummaryDEGs[21,2] <- sum(na.omit(results(ddsMic_B)$padj) < 0.05)
SummaryDEGs[21,3] <- sum(na.omit(results(ddsMic_C)$padj) < 0.05)
SummaryDEGs[21,4] <- sum(na.omit(results(ddsMic_O)$padj) < 0.05)
SummaryDEGs[22,1] <- sum(na.omit(results(dds_mgA)$padj) < 0.05)
SummaryDEGs[22,2] <- sum(na.omit(results(dds_mgB)$padj) < 0.05)
SummaryDEGs[22,3] <- sum(na.omit(results(dds_mgC)$padj) < 0.05)
SummaryDEGs[22,4] <- sum(na.omit(results(dds_mgO)$padj) < 0.05)
SummaryDEGs[23,1] <- sum(na.omit(results(dds_tgA)$padj) < 0.05)
SummaryDEGs[23,2] <- sum(na.omit(results(dds_tgB)$padj) < 0.05)
SummaryDEGs[23,3] <- sum(na.omit(results(dds_tgC)$padj) < 0.05)
SummaryDEGs[23,4] <- sum(na.omit(results(dds_tgO)$padj) < 0.05)
SummaryDEGs[24,1] <- sum(na.omit(results(dds_agA)$padj) < 0.05)
SummaryDEGs[24,2] <- sum(na.omit(results(dds_agB)$padj) < 0.05)
SummaryDEGs[24,3] <- sum(na.omit(results(dds_agC)$padj) < 0.05)
SummaryDEGs[24,4] <- sum(na.omit(results(dds_agO)$padj) < 0.05)
SummaryDEGs[25,1] <- sum(na.omit(results(dds_pgA)$padj) < 0.05)
SummaryDEGs[25,2] <- sum(na.omit(results(dds_pgB)$padj) < 0.05)
SummaryDEGs[25,3] <- sum(na.omit(results(dds_pgC)$padj) < 0.05)
SummaryDEGs[25,4] <- sum(na.omit(results(dds_pgO)$padj) < 0.05)

write.csv(SummaryDEGs, "DEGs_Controlling.csv")