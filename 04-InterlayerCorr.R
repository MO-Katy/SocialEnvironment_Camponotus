# Combine all
library('stringr');library('DESeq2');library('ggplot2');library('readr');
library('grid');library('gridExtra');library('igraph'); library('lme4')

# Set directories
MAINDIR  <- "/Users/tkay/Desktop/SNG/data"
FIGDIR   <- "/Users/tkay/Desktop/SNG/figures"

setwd(MAINDIR)
COLONIES <- c("A", "B", "C", "O")
for (col in COLONIES){
  assign(paste0("node_", col), read.csv(paste0("md_", col, "4.csv"), row.names=NULL))
}

DFA = DFB = DFC = DFO = as.data.frame(matrix(nrow = 6,ncol = 6,NA))
colnames(DFA) = colnames(DFB) = colnames(DFC) = colnames(DFO) =
  rownames(DFA) = rownames(DFB) = rownames(DFC) = rownames(DFO) = 
  c("Age", "Behaviour", "Physical environment", "Social", "Gene expression", "Microbiota")

DEGs <- read.csv("DEGs.csv", row.names = 1)[-5,]

DFA[1,2] = DFA[2,1] = round(summary(lm(node_A$Age~ node_A$taskPCA))$r.squared,3)
DFA[1,3] = DFA[3,1] = round(summary(lm(node_A$Age~ node_A$spacePCA))$r.squared,3)
DFA[2,3] = DFA[3,2] = round(summary(lm(node_A$taskPCA~ node_A$spacePCA))$r.squared,3)
DFA[1,4] = DFA[4,1] = round(summary(lm(node_A$Age~ node_A$maturity))$r.squared,3)
DFA[2,4] = DFA[4,2] = round(summary(lm(node_A$taskPCA~ node_A$maturity))$r.squared,3)
DFA[3,4] = DFA[4,3] = round(summary(lm(node_A$Age~ node_A$maturity))$r.squared,3)
DFA[1,5] = DFA[5,1] = round(summary(lm(node_A$Age~ node_A$genePC2))$r.squared,3)
DFA[2,5] = DFA[5,2] = round(summary(lm(node_A$taskPCA~ node_A$genePC2))$r.squared,3)
DFA[3,5] = DFA[5,3] = round(summary(lm(node_A$spacePCA~ node_A$genePC2))$r.squared,3)
DFA[4,5] = DFA[5,4] = round(summary(lm(node_A$maturity~ node_A$genePC2))$r.squared,3)
DFA[1,6] = DFA[6,1] = round(summary(lm(node_A$Age~ node_A$MicPCAm))$r.squared,3)
DFA[2,6] = DFA[6,2] = round(summary(lm(node_A$taskPCA~ node_A$MicPCAm))$r.squared,3)
DFA[3,6] = DFA[6,3] = round(summary(lm(node_A$spacePCA~ node_A$MicPCAm))$r.squared,3)
DFA[4,6] = DFA[6,4] = round(summary(lm(node_A$maturity~ node_A$MicPCAm))$r.squared,3)
DFA[5,6] = DFA[6,5] = round(summary(lm(node_A$genePC2~ node_A$MicPCAm))$r.squared,3)

DFB[1,2] = DFB[2,1] = round(summary(lm(node_B$Age~ node_B$taskPCA))$r.squared,3)
DFB[1,3] = DFB[3,1] = round(summary(lm(node_B$Age~ node_B$spacePCA))$r.squared,3)
DFB[2,3] = DFB[3,2] = round(summary(lm(node_B$taskPCA~ node_B$spacePCA))$r.squared,3)
DFB[1,4] = DFB[4,1] = round(summary(lm(node_B$Age~ node_B$maturity))$r.squared,3)
DFB[2,4] = DFB[4,2] = round(summary(lm(node_B$taskPCA~ node_B$maturity))$r.squared,3)
DFB[3,4] = DFB[4,3] = round(summary(lm(node_B$Age~ node_B$maturity))$r.squared,3)
DFB[1,5] = DFB[5,1] = round(summary(lm(node_B$Age~ node_B$genePC2))$r.squared,3)
DFB[2,5] = DFB[5,2] = round(summary(lm(node_B$taskPCA~ node_B$genePC2))$r.squared,3)
DFB[3,5] = DFB[5,3] = round(summary(lm(node_B$spacePCA~ node_B$genePC2))$r.squared,3)
DFB[4,5] = DFB[5,4] = round(summary(lm(node_B$maturity~ node_B$genePC2))$r.squared,3)
DFB[1,6] = DFB[6,1] = round(summary(lm(node_B$Age~ node_B$MicPCAm))$r.squared,3)
DFB[2,6] = DFB[6,2] = round(summary(lm(node_B$taskPCA~ node_B$MicPCAm))$r.squared,3)
DFB[3,6] = DFB[6,3] = round(summary(lm(node_B$spacePCA~ node_B$MicPCAm))$r.squared,3)
DFB[4,6] = DFB[6,4] = round(summary(lm(node_B$maturity~ node_B$MicPCAm))$r.squared,3)
DFB[5,6] = DFB[6,5] = round(summary(lm(node_B$genePC2~ node_B$MicPCAm))$r.squared,3)

DFC[1,2] = DFC[2,1] = round(summary(lm(node_C$Age~ node_C$taskPCA))$r.squared,3)
DFC[1,3] = DFC[3,1] = round(summary(lm(node_C$Age~ node_C$spacePCA))$r.squared,3)
DFC[2,3] = DFC[3,2] = round(summary(lm(node_C$taskPCA~ node_C$spacePCA))$r.squared,3)
DFC[1,4] = DFC[4,1] = round(summary(lm(node_C$Age~ node_C$maturity))$r.squared,3)
DFC[2,4] = DFC[4,2] = round(summary(lm(node_C$taskPCA~ node_C$maturity))$r.squared,3)
DFC[3,4] = DFC[4,3] = round(summary(lm(node_C$Age~ node_C$maturity))$r.squared,3)
DFC[1,5] = DFC[5,1] = round(summary(lm(node_C$Age~ node_C$genePC2))$r.squared,3)
DFC[2,5] = DFC[5,2] = round(summary(lm(node_C$taskPCA~ node_C$genePC2))$r.squared,3)
DFC[3,5] = DFC[5,3] = round(summary(lm(node_C$spacePCA~ node_C$genePC2))$r.squared,3)
DFC[4,5] = DFC[5,4] = round(summary(lm(node_C$maturity~ node_C$genePC2))$r.squared,3)
DFC[1,6] = DFC[6,1] = round(summary(lm(node_C$Age~ node_C$MicPCAm))$r.squared,3)
DFC[2,6] = DFC[6,2] = round(summary(lm(node_C$taskPCA~ node_C$MicPCAm))$r.squared,3)
DFC[3,6] = DFC[6,3] = round(summary(lm(node_C$spacePCA~ node_C$MicPCAm))$r.squared,3)
DFC[4,6] = DFC[6,4] = round(summary(lm(node_C$maturity~ node_C$MicPCAm))$r.squared,3)
DFC[5,6] = DFC[6,5] = round(summary(lm(node_C$genePC2~ node_C$MicPCAm))$r.squared,3)

DFO[1,2] = DFO[2,1] = round(summary(lm(node_O$Age~ node_O$taskPCA))$r.squared,3)
DFO[1,3] = DFO[3,1] = round(summary(lm(node_O$Age~ node_O$spacePCA))$r.squared,3)
DFO[2,3] = DFO[3,2] = round(summary(lm(node_O$taskPCA~ node_O$spacePCA))$r.squared,3)
DFO[1,4] = DFO[4,1] = round(summary(lm(node_O$Age~ node_O$maturity))$r.squared,3)
DFO[2,4] = DFO[4,2] = round(summary(lm(node_O$taskPCA~ node_O$maturity))$r.squared,3)
DFO[3,4] = DFO[4,3] = round(summary(lm(node_O$Age~ node_O$maturity))$r.squared,3)
DFO[1,5] = DFO[5,1] = round(summary(lm(node_O$Age~ node_O$genePC2))$r.squared,3)
DFO[2,5] = DFO[5,2] = round(summary(lm(node_O$taskPCA~ node_O$genePC2))$r.squared,3)
DFO[3,5] = DFO[5,3] = round(summary(lm(node_O$spacePCA~ node_O$genePC2))$r.squared,3)
DFO[4,5] = DFO[5,4] = round(summary(lm(node_O$maturity~ node_O$genePC2))$r.squared,3)
DFO[1,6] = DFO[6,1] = round(summary(lm(node_O$Age~ node_O$MicPCAm))$r.squared,3)
DFO[2,6] = DFO[6,2] = round(summary(lm(node_O$taskPCA~ node_O$MicPCAm))$r.squared,3)
DFO[3,6] = DFO[6,3] = round(summary(lm(node_O$spacePCA~ node_O$MicPCAm))$r.squared,3)
DFO[4,6] = DFO[6,4] = round(summary(lm(node_O$maturity~ node_O$MicPCAm))$r.squared,3)
DFO[5,6] = DFO[6,5] = round(summary(lm(node_O$genePC2~ node_O$MicPCAm))$r.squared,3)


DFALL <- round((DFO + DFA + DFC + DFB)/4, 2)
DFALL[is.na(DFALL)] <- 0


gr <- graph_from_adjacency_matrix(as.matrix(DFALL), weighted = TRUE, mode = c("undirected"))
layout_gr <- layout_with_fr(gr)
#E(gr)$color <- gray(1-(E(gr)$weight/max(E(gr)$weight, na.rm = 1)))
E(gr)$color <- "darkgray"

setwd(FIGDIR)
pdf('Multinet.pdf', width=30, height=20)
par(mfrow = c(1,1), mar = c(0,0,0,0), oma = c(0,0,0,0), bty="n", mgp = c(0,0,0), mai = c(0,0,0,0), family = "serif")
plot(gr, layout = layout_gr,
     vertex.size=strength(gr)*5, vertex.color = c("blue", "green", "orange", "pink", "purple", "red"),
     edge.width = ((E(gr)$weight/max(E(gr)$weight)))*40,
     edge.label = E(gr)$weight, edge.label.cex = 5,
     #vertex.label = NA,
     vertex.label.cex = 5,
     label.color = "black",
     edge.label.color = "black")
dev.off()


# Supplementary comparisons of correlations between PC1 and PC2 for different variables
mean(c(summary(lm(node_A$Age~ node_A$genePC1))$r.squared,
       summary(lm(node_B$Age~ node_B$genePC1))$r.squared,
       summary(lm(node_C$Age~ node_C$genePC1))$r.squared,
       summary(lm(node_O$Age~ node_O$genePC1))$r.squared))

mean(c(summary(lm(node_A$Age~ node_A$genePC2))$r.squared,
       summary(lm(node_B$Age~ node_B$genePC2))$r.squared,
       summary(lm(node_C$Age~ node_C$genePC2))$r.squared,
       summary(lm(node_O$Age~ node_O$genePC2))$r.squared))

mean(c(summary(lm(node_A$spacePCA~ node_A$genePC1))$r.squared,
       summary(lm(node_B$spacePCA~ node_B$genePC1))$r.squared,
       summary(lm(node_C$spacePCA~ node_C$genePC1))$r.squared,
       summary(lm(node_O$spacePCA~ node_O$genePC1))$r.squared))

mean(c(summary(lm(node_A$spacePCA~ node_A$genePC2))$r.squared,
       summary(lm(node_B$spacePCA~ node_B$genePC2))$r.squared,
       summary(lm(node_C$spacePCA~ node_C$genePC2))$r.squared,
       summary(lm(node_O$spacePCA~ node_O$genePC2))$r.squared))

mean(c(summary(lm(node_A$MicPCAm~ node_A$genePC1))$r.squared,
       summary(lm(node_B$MicPCAm~ node_B$genePC1))$r.squared,
       summary(lm(node_C$MicPCAm~ node_C$genePC1))$r.squared,
       summary(lm(node_O$MicPCAm~ node_O$genePC1))$r.squared))

mean(c(summary(lm(node_A$MicPCAm~ node_A$genePC2))$r.squared,
       summary(lm(node_B$MicPCAm~ node_B$genePC2))$r.squared,
       summary(lm(node_C$MicPCAm~ node_C$genePC2))$r.squared,
       summary(lm(node_O$MicPCAm~ node_O$genePC2))$r.squared))

mean(c(summary(lm(node_A$maturity~ node_A$genePC1))$r.squared,
       summary(lm(node_B$maturity~ node_B$genePC1))$r.squared,
       summary(lm(node_C$maturity~ node_C$genePC1))$r.squared,
       summary(lm(node_O$maturity~ node_O$genePC1))$r.squared))

mean(c(summary(lm(node_A$maturity~ node_A$genePC2))$r.squared,
       summary(lm(node_B$maturity~ node_B$genePC2))$r.squared,
       summary(lm(node_C$maturity~ node_C$genePC2))$r.squared,
       summary(lm(node_O$maturity~ node_O$genePC2))$r.squared))

mean(c(summary(lm(node_A$taskPCA~ node_A$genePC1))$r.squared,
       summary(lm(node_B$taskPCA~ node_B$genePC1))$r.squared,
       summary(lm(node_C$taskPCA~ node_C$genePC1))$r.squared,
       summary(lm(node_O$taskPCA~ node_O$genePC1))$r.squared))

mean(c(summary(lm(node_A$taskPCA~ node_A$genePC2))$r.squared,
       summary(lm(node_B$taskPCA~ node_B$genePC2))$r.squared,
       summary(lm(node_C$taskPCA~ node_C$genePC2))$r.squared,
       summary(lm(node_O$taskPCA~ node_O$genePC2))$r.squared))


mean(c(summary(lm(node_A$genePC1~ node_A$Extraction))$r.squared,
       summary(lm(node_B$genePC1~ node_B$Extraction))$r.squared,
       summary(lm(node_C$genePC1~ node_C$Extraction))$r.squared))

mean(c(summary(lm(node_A$genePC2~ node_A$Extraction))$r.squared,
       summary(lm(node_B$genePC2~ node_B$Extraction))$r.squared,
       summary(lm(node_C$genePC2~ node_C$Extraction))$r.squared))
     
t.test(DEGs$Social, DEGs$Task, paired = 1)

# Mutlinet plots for each of the four colonies separately

DFA[is.na(DFA)] <- 0
DFB[is.na(DFB)] <- 0
DFC[is.na(DFC)] <- 0
DFO[is.na(DFO)] <- 0
gr_A <- graph_from_adjacency_matrix(as.matrix(DFA), weighted = TRUE, mode = c("undirected"))
gr_B <- graph_from_adjacency_matrix(as.matrix(DFB), weighted = TRUE, mode = c("undirected"))
gr_C <- graph_from_adjacency_matrix(as.matrix(DFC), weighted = TRUE, mode = c("undirected"))
gr_O <- graph_from_adjacency_matrix(as.matrix(DFO), weighted = TRUE, mode = c("undirected"))

layout_gr_A <- layout_with_fr(gr_A)
layout_gr_B <- layout_with_fr(gr_B)
layout_gr_C <- layout_with_fr(gr_C)
layout_gr_O <- layout_with_fr(gr_O)

E(gr_A)$color = E(gr_B)$color = E(gr_C)$color = E(gr_O)$color = "darkgray"

setwd(FIGDIR)
pdf('Multinet_A.pdf', width=30, height=20)
par(mfrow = c(1,1), mar = c(0,0,0,0), oma = c(0,0,0,0), bty="n", mgp = c(0,0,0), mai = c(0,0,0,0), family = "serif")
plot(gr_A, layout = layout_gr_A,
     vertex.size=1, vertex.color = "blue",
     edge.width = ((E(gr_A)$weight/max(E(gr_A)$weight)))*40,
     edge.label = E(gr_A)$weight, edge.label.cex = 5,
     vertex.label.cex = 7,
     label.color = "black",
     edge.label.color = "black")
dev.off()
pdf('Multinet_B.pdf', width=30, height=20)
par(mfrow = c(1,1), mar = c(0,0,0,0), oma = c(0,0,0,0), bty="n", mgp = c(0,0,0), mai = c(0,0,0,0), family = "serif")
plot(gr_B, layout = layout_gr_B,
     vertex.size=1, vertex.color = "blue",
     edge.width = ((E(gr_B)$weight/max(E(gr_B)$weight)))*40,
     edge.label = E(gr_B)$weight, edge.label.cex = 5,
     vertex.label.cex = 7,
     label.color = "black",
     edge.label.color = "black")
dev.off()
pdf('Multinet_C.pdf', width=30, height=20)
par(mfrow = c(1,1), mar = c(0,0,0,0), oma = c(0,0,0,0), bty="n", mgp = c(0,0,0), mai = c(0,0,0,0), family = "serif")
plot(gr_C, layout = layout_gr_C,
     vertex.size=1, vertex.color = "blue",
     edge.width = ((E(gr_C)$weight/max(E(gr_C)$weight)))*40,
     edge.label = E(gr_C)$weight, edge.label.cex = 5,
     vertex.label.cex = 7,
     label.color = "black",
     edge.label.color = "black")
dev.off()
pdf('Multinet_O.pdf', width=30, height=20)
par(mfrow = c(1,1), mar = c(0,0,0,0), oma = c(0,0,0,0), bty="n", mgp = c(0,0,0), mai = c(0,0,0,0), family = "serif")
plot(gr_O, layout = layout_gr_O,
     vertex.size=1, vertex.color = "blue",
     edge.width = ((E(gr_O)$weight/max(E(gr_O)$weight)))*40,
     edge.label = E(gr_O)$weight, edge.label.cex = 5,
     vertex.label.cex = 7,
     label.color = "black",
     edge.label.color = "black")
dev.off()

# Each correlation value when controlling for each other variable
names(node_A)[15] <- "social"
names(node_B)[15] <- "social"
names(node_C)[15] <- "social"
names(node_O)[14] <- "social"
names(node_A)[9] <- "behavior"
names(node_B)[9] <- "behavior"
names(node_C)[9] <- "behavior"
names(node_O)[9] <- "behavior"
names(node_A)[10] <- "physical"
names(node_B)[10] <- "physical"
names(node_C)[10] <- "physical"
names(node_O)[10] <- "physical"
names(node_A)[12] <- "age"
names(node_B)[12] <- "age"
names(node_C)[12] <- "age"
names(node_O)[12] <- "age"
names(node_A)[17] <- "microbiota"
names(node_B)[17] <- "microbiota"
names(node_C)[17] <- "microbiota"
names(node_O)[16] <- "microbiota"

anova(lm(node_A$social ~ node_A$behavior))$`Mean Sq`[1]
anova(lm(node_A$social ~ node_A$age))$`Mean Sq`[1]
anova(lm(node_A$social ~ node_A$physical))$`Mean Sq`[1]
anova(lm(node_A$social ~ node_A$microbiota))$`Mean Sq`[1]
anova(lm(node_A$social ~ node_A$genePC2))$`Mean Sq`[1]