# Predict aspects of behavior/ age/ social position from brain gene expression data

# Set directories and load libraries
MAINDIR  <- "/Users/tkay/Desktop/Work/SNG/data"
FIGDIR   <- "/Users/tkay/Desktop/Work/SNG/figures"
FIGDATADIR   <- "/Users/tkay/Desktop/Work/SNG/figures_data"
library('stringr');library('DESeq2');library('ggplot2');library('readr');
library('grid');library('gridExtra');library('igraph'); library('lme4');
library('randomForest'); library('e1071')

# Set parameters
SPLIT <- 50 #T he percentage of ants used for training
P <- 0.05 # P-value cut off for DESeq2
LOW.COUNT.FILTER <- 100 # Low counts filtering of transcriptomes
ITERATIONS <- 100 # How many interations to run the script for

# Import data
setwd(MAINDIR)
cts  <- read.csv("brains_ABC.csv")[-1,]
cts_O <- read.csv("brains_O.csv", check.names = F)
COLONIES <- c("A", "B", "C", "O")
for (col in COLONIES){
  assign(paste0("node_", col), read.csv(paste0("md_", col, "4.csv"), row.names=NULL))
}

# Format transcriptomic data
newnames <- c() # Rename
for (i in 2:length(colnames(cts))){
  name1 <- unlist(str_split(colnames(cts)[[i]], "_"))[3]
  ID <- parse_number(name1)
  Col <- substr(name1,1,1)
  Ind <- paste(Col, ID, sep = "_")
  newnames <- c(newnames, Ind)
}
colnames(cts) <- c("geneID", newnames) # Change row names to geneID
rownames(cts) <- cts[,1]; cts <- cts[,-1]

ColID <- data.frame(ID = colnames(cts), col = substr(colnames(cts),1,1))

# Separate transcriptomes by colony
for (col in COLONIES[1:3]){
  assign(paste0("cts_", col), cts[,which(ColID$col == col)])
}

# Rename individuals for comparison across data-sets
colnames(cts_A) <- parse_number(colnames(cts_A))-1
colnames(cts_B) <- parse_number(colnames(cts_B))-1
colnames(cts_C) <- parse_number(colnames(cts_C))-1

# Remove IDs which are not present in the individual meta-data
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

# reformat old brains data
cts_O <- cts_O[,!colnames(cts_O) == "x"]
cts_O <- cts_O[,order(as.numeric(colnames(cts_O)))]
cts_O <- cts_O[,colnames(cts_O) %in% node_O$antID]
rownames(cts_O) <- rownames(cts_A)

# Name rows
row.names(node_A) <- node_A$antID
row.names(node_B) <- node_B$antID
row.names(node_C) <- node_C$antID
row.names(node_O) <- node_O$antID

# Remove one from A where there is no microbiota data
node_A <- na.omit(node_A)

# Trim and uniformize node data
cols <- c("antID", "Age", "cleaning", "broodcare", "guarding", "queenattending",
           "maturity", "prop", "taskPCA_ALL", "spacePCA", "MicPCAm")
node_A <- node_A[,colnames(node_A) %in% cols]
node_B <- node_B[,colnames(node_B) %in% cols]
node_C <- node_C[,colnames(node_C) %in% cols]
node_O <- node_O[,colnames(node_O) %in% cols]

# give sample names with colony letter
rownames(node_A) <- paste0("A_", rownames(node_A))
rownames(node_B) <- paste0("B_", rownames(node_B))
rownames(node_C) <- paste0("C_", rownames(node_C))
rownames(node_O) <- paste0("O_", rownames(node_O))

node_A$colony <- "A"
node_B$colony <- "B"
node_C$colony <- "C"
node_O$colony <- "O"

# bind node data
node_ALL <- rbind(node_A, node_B, node_C, node_O)

# bind all counts data
cts_A <- cts_A[,order(as.numeric(colnames(cts_A)))]
colnames(cts_A) <- paste0("A_", colnames(cts_A))
colnames(cts_B) <- paste0("B_", colnames(cts_B))
colnames(cts_C) <- paste0("C_", colnames(cts_C))
colnames(cts_O) <- paste0("O_", colnames(cts_O))

cts_ALL <- cbind(cts_A, cts_B, cts_C, cts_O)

maturityR <- c()
maturityO  <- c()
for (it in 1:ITERATIONS){
  print(it)
  # Define row numbers for the training and test dataset
  train <- sample(1:ncol(cts_ALL), ncol(cts_ALL)/(100/SPLIT), replace = FALSE)
  test <- seq(from = 1, to = ncol(cts_ALL), by = 1)
  test <- test[!test%in%train]
  
  # Extract these rows from the count data and metadata
  cts.train <- cts_ALL[,train]
  cts.test <- cts_ALL[,test]
  coldata.train <- node_ALL[train,]
  coldata.test <- node_ALL[test,]
  
  # Keep genes differentially expressed by module in training data set
  forDESeq <- cts.train
  DESeqColdata <- coldata.train
  
  dds <- DESeqDataSetFromMatrix(countData = forDESeq,
                                colData = DESeqColdata,
                                design = ~ colony + cluster_1)
  keep <- rowSums(counts(dds)) >= LOW.COUNT.FILTER
  dds <- dds[keep,]
  vsd <- vst(dds, blind = FALSE)
  DESeqNorm <- as.data.frame(assay(vsd))
  
  dds <- DESeq(dds)
  res <- na.omit(as.data.frame(results(dds)))
  ressig <- res[res$padj<= P,]
  DEGenes <- rownames(ressig)
  
  # Test and training data-sets with only these genes
  cts.train <- cts.train[rownames(cts.train)%in% DEGenes,]
  cts.test <- cts.test[rownames(cts.test)%in% DEGenes,]
  
  # Support vector machine model
  # Create training and test datasets for svm input
  SVM.train <- rbind(cts.train, maturity = coldata.train$cluster_1)
  SVM.test <- rbind(cts.test, maturity = coldata.test$cluster_1)
  
  SVM.train <- as.data.frame(t(SVM.train))
  SVM.test <- as.data.frame(t(SVM.test))
  
  # SVM model
  classifier = svm(formula = maturity ~ .,
                   data = SVM.train,
                   type = 'eps-regression',
                   kernel = 'linear')
  
  # Predict for test data-set
  y_pred <- predict(classifier, newdata = SVM.test[-ncol(SVM.test)])
  
  # Evaluate accuracy
  cm <- as.data.frame(cbind(real = SVM.test[,ncol(SVM.test)], pred = y_pred))
  out <- 0
  if (max(cm$pred) > (mean(cm$pred) + sd(cm$pred)*2)) {
     out <- out + sum(cm$pred > (mean(cm$pred) + sd(cm$pred)*2))
    cm <- cm[-which(cm$pred > (mean(cm$pred) + sd(cm$pred)*2)),]
  }
  if (min(cm$pred) < (mean(cm$pred) - sd(cm$pred)*2)){
    out <- out + sum(cm$pred < (mean(cm$pred) - sd(cm$pred)*2))
    cm <- cm[-which(cm$pred < (mean(cm$pred) - sd(cm$pred)*2)),]
  }
  maturityO <- c(maturityO, out)
  
  plot(cm$pred, cm$real)
  
  maturityR <- c(maturityR, summary(lm(cm$real~ cm$pred))$r.squared)
}

write.csv(maturityR, "PredictionMaturity.csv")
write.csv(maturityO, "OutlierMaturity.csv")

propR <- c()
propO  <- c()
for (it in 1:ITERATIONS){
  print(it)
  train <- sample(1:ncol(cts_ALL), ncol(cts_ALL)/(100/SPLIT), replace = FALSE)
  test <- seq(from = 1, to = ncol(cts_ALL), by = 1)
  test <- test[!test%in%train]
  cts.train <- cts_ALL[,train]
  cts.test <- cts_ALL[,test]
  coldata.train <- node_ALL[train,]
  coldata.test <- node_ALL[test,]
  forDESeq <- cts.train
  DESeqColdata <- coldata.train
  
  dds <- DESeqDataSetFromMatrix(countData = forDESeq,
                                colData = DESeqColdata,
                                design = ~ colony + prop)
  
  keep <- rowSums(counts(dds)) >= LOW.COUNT.FILTER
  dds <- dds[keep,]
  vsd <- vst(dds, blind = FALSE)
  DESeqNorm <- as.data.frame(assay(vsd))
  dds <- DESeq(dds)
  res <- na.omit(as.data.frame(results(dds)))
  ressig <- res[res$padj<= P,]
  DEGenes <- rownames(ressig)
  cts.train <- cts.train[rownames(cts.train)%in% DEGenes,]
  cts.test <- cts.test[rownames(cts.test)%in% DEGenes,]
  SVM.train <- rbind(cts.train, prop = coldata.train$prop)
  SVM.test <- rbind(cts.test, prop = coldata.test$prop)
  SVM.train <- as.data.frame(t(SVM.train))
  SVM.test <- as.data.frame(t(SVM.test))
  classifier = svm(formula = prop ~ .,
                   data = SVM.train,
                   type = 'eps-regression',
                   kernel = 'linear')
  y_pred <- predict(classifier, newdata = SVM.test[-ncol(SVM.test)])
  cm <- as.data.frame(cbind(real = SVM.test[,ncol(SVM.test)], pred = y_pred))
  out <- 0
  if (max(cm$pred) > (mean(cm$pred) + sd(cm$pred)*2)) {
    out <- out + sum(cm$pred > (mean(cm$pred) + sd(cm$pred)*2))
    cm <- cm[-which(cm$pred > (mean(cm$pred) + sd(cm$pred)*2)),]
  }
  if (min(cm$pred) < (mean(cm$pred) - sd(cm$pred)*2)){
    out <- out + sum(cm$pred < (mean(cm$pred) - sd(cm$pred)*2))
    cm <- cm[-which(cm$pred < (mean(cm$pred) - sd(cm$pred)*2)),]
  }
  propO <- c(propO, out)
  plot(cm$pred, cm$real)
  propR <- c(propR, summary(lm(cm$real~ cm$pred))$r.squared)
}
write.csv(propR, "PredictionProp.csv")
write.csv(propO, "OutlierProp.csv")

ageR <- c()
ageO  <- c()
for (it in 1:ITERATIONS){
  print(it)
  train <- sample(1:ncol(cts_ALL), ncol(cts_ALL)/(100/SPLIT), replace = FALSE)
  test <- seq(from = 1, to = ncol(cts_ALL), by = 1)
  test <- test[!test%in%train]
  cts.train <- cts_ALL[,train]
  cts.test <- cts_ALL[,test]
  coldata.train <- node_ALL[train,]
  coldata.test <- node_ALL[test,]
  forDESeq <- cts.train
  DESeqColdata <- coldata.train
  
  dds <- DESeqDataSetFromMatrix(countData = forDESeq,
                                colData = DESeqColdata,
                                design = ~ colony + Age)
  
  keep <- rowSums(counts(dds)) >= LOW.COUNT.FILTER
  dds <- dds[keep,]
  vsd <- vst(dds, blind = FALSE)
  DESeqNorm <- as.data.frame(assay(vsd))
  dds <- DESeq(dds)
  res <- na.omit(as.data.frame(results(dds)))
  ressig <- res[res$padj<= P,]
  DEGenes <- rownames(ressig)
  cts.train <- cts.train[rownames(cts.train)%in% DEGenes,]
  cts.test <- cts.test[rownames(cts.test)%in% DEGenes,]
  SVM.train <- rbind(cts.train, age = coldata.train$Age)
  SVM.test <- rbind(cts.test, age = coldata.test$Age)
  SVM.train <- as.data.frame(t(SVM.train))
  SVM.test <- as.data.frame(t(SVM.test))
  classifier = svm(formula = age ~ .,
                   data = SVM.train,
                   type = 'eps-regression',
                   kernel = 'linear')
  y_pred <- predict(classifier, newdata = SVM.test[-ncol(SVM.test)])
  cm <- as.data.frame(cbind(real = SVM.test[,ncol(SVM.test)], pred = y_pred))
  out <- 0
  if (max(cm$pred) > (mean(cm$pred) + sd(cm$pred)*2)) {
    out <- out + sum(cm$pred > (mean(cm$pred) + sd(cm$pred)*2))
    cm <- cm[-which(cm$pred > (mean(cm$pred) + sd(cm$pred)*2)),]
  }
  if (min(cm$pred) < (mean(cm$pred) - sd(cm$pred)*2)){
    out <- out + sum(cm$pred < (mean(cm$pred) - sd(cm$pred)*2))
    cm <- cm[-which(cm$pred < (mean(cm$pred) - sd(cm$pred)*2)),]
  }
  ageO <- c(ageO, out)
  plot(cm$pred, cm$real)
  ageR <- c(ageR, summary(lm(cm$real~ cm$pred))$r.squared)
}
write.csv(ageR, "PredictionAge.csv")
write.csv(ageO, "OutlierAge.csv")

nursingR <- c()
nursingO  <- c()
for (it in 1:ITERATIONS){
  print(it)
  train <- sample(1:ncol(cts_ALL), ncol(cts_ALL)/(100/SPLIT), replace = FALSE)
  test <- seq(from = 1, to = ncol(cts_ALL), by = 1)
  test <- test[!test%in%train]
  cts.train <- cts_ALL[,train]
  cts.test <- cts_ALL[,test]
  coldata.train <- node_ALL[train,]
  coldata.test <- node_ALL[test,]
  forDESeq <- cts.train
  DESeqColdata <- coldata.train
  
  dds <- DESeqDataSetFromMatrix(countData = forDESeq,
                                colData = DESeqColdata,
                                design = ~ colony + broodcare)
  
  keep <- rowSums(counts(dds)) >= LOW.COUNT.FILTER
  dds <- dds[keep,]
  vsd <- vst(dds, blind = FALSE)
  DESeqNorm <- as.data.frame(assay(vsd))
  dds <- DESeq(dds)
  res <- na.omit(as.data.frame(results(dds)))
  ressig <- res[res$padj<= P,]
  DEGenes <- rownames(ressig)
  cts.train <- cts.train[rownames(cts.train)%in% DEGenes,]
  cts.test <- cts.test[rownames(cts.test)%in% DEGenes,]
  SVM.train <- rbind(cts.train, nursing = coldata.train$broodcare)
  SVM.test <- rbind(cts.test, nursing = coldata.test$broodcare)
  SVM.train <- as.data.frame(t(SVM.train))
  SVM.test <- as.data.frame(t(SVM.test))
  classifier = svm(formula = nursing ~ .,
                   data = SVM.train,
                   type = 'eps-regression',
                   kernel = 'linear')
  y_pred <- predict(classifier, newdata = SVM.test[-ncol(SVM.test)])
  cm <- as.data.frame(cbind(real = SVM.test[,ncol(SVM.test)], pred = y_pred))
  out <- 0
  if (max(cm$pred) > (mean(cm$pred) + sd(cm$pred)*2)) {
    out <- out + sum(cm$pred > (mean(cm$pred) + sd(cm$pred)*2))
    cm <- cm[-which(cm$pred > (mean(cm$pred) + sd(cm$pred)*2)),]
  }
  if (min(cm$pred) < (mean(cm$pred) - sd(cm$pred)*2)){
    out <- out + sum(cm$pred < (mean(cm$pred) - sd(cm$pred)*2))
    cm <- cm[-which(cm$pred < (mean(cm$pred) - sd(cm$pred)*2)),]
  }
  nursingO <- c(nursingO, out)
  plot(cm$pred, cm$real)
  nursingR <- c(nursingR, summary(lm(cm$real~ cm$pred))$r.squared)
}
write.csv(nursingR, "PredictionNursing.csv")
write.csv(nursingO, "OutlierNursing.csv")

cleaningR <- c()
cleaningO  <- c()
for (it in 1:ITERATIONS){
  print(it)
  train <- sample(1:ncol(cts_ALL), ncol(cts_ALL)/(100/SPLIT), replace = FALSE)
  test <- seq(from = 1, to = ncol(cts_ALL), by = 1)
  test <- test[!test%in%train]
  cts.train <- cts_ALL[,train]
  cts.test <- cts_ALL[,test]
  coldata.train <- node_ALL[train,]
  coldata.test <- node_ALL[test,]
  forDESeq <- cts.train
  DESeqColdata <- coldata.train
  
  dds <- DESeqDataSetFromMatrix(countData = forDESeq,
                                colData = DESeqColdata,
                                design = ~ colony + cleaning)
  
  keep <- rowSums(counts(dds)) >= LOW.COUNT.FILTER
  dds <- dds[keep,]
  vsd <- vst(dds, blind = FALSE)
  DESeqNorm <- as.data.frame(assay(vsd))
  dds <- DESeq(dds)
  res <- na.omit(as.data.frame(results(dds)))
  ressig <- res[res$padj<= P,]
  DEGenes <- rownames(ressig)
  cts.train <- cts.train[rownames(cts.train)%in% DEGenes,]
  cts.test <- cts.test[rownames(cts.test)%in% DEGenes,]
  SVM.train <- rbind(IDs = cts.train, cleaning = coldata.train$cleaning)
  SVM.test <- rbind(IDs =cts.test, cleaning = coldata.test$cleaning)
  SVM.train <- as.data.frame(t(SVM.train))
  SVM.test <- as.data.frame(t(SVM.test))
  if(length(DEGenes > 0)){
    classifier = svm(formula = cleaning ~ .,
                     data = SVM.train,
                     type = 'eps-regression',
                     kernel = 'linear')
    y_pred <- predict(classifier, newdata = SVM.test[-ncol(SVM.test)])
    cm <- as.data.frame(cbind(real = SVM.test[,ncol(SVM.test)], pred = y_pred))
    out <- 0
    if (max(cm$pred) > (mean(cm$pred) + sd(cm$pred)*2)) {
      out <- out + sum(cm$pred > (mean(cm$pred) + sd(cm$pred)*2))
      cm <- cm[-which(cm$pred > (mean(cm$pred) + sd(cm$pred)*2)),]
    }
    if (min(cm$pred) < (mean(cm$pred) - sd(cm$pred)*2)){
      out <- out + sum(cm$pred < (mean(cm$pred) - sd(cm$pred)*2))
      cm <- cm[-which(cm$pred < (mean(cm$pred) - sd(cm$pred)*2)),]
    }
    cleaningO <- c(cleaningO, out)
    plot(cm$pred, cm$real)
    cleaningR <- c(cleaningR, summary(lm(cm$real~ cm$pred))$r.squared)
  } else{
    cleaningR <- c(cleaningR, NA)
    cleaningO <- c(cleaningO, NA)
  }
}
write.csv(cleaningR, "PredictionCleaning.csv")
write.csv(cleaningO, "OutlierCleaning.csv")

guardingR <- c()
guardingO  <- c()
for (it in 1:ITERATIONS){
  print(it)
  train <- sample(1:ncol(cts_ALL), ncol(cts_ALL)/(100/SPLIT), replace = FALSE)
  test <- seq(from = 1, to = ncol(cts_ALL), by = 1)
  test <- test[!test%in%train]
  cts.train <- cts_ALL[,train]
  cts.test <- cts_ALL[,test]
  coldata.train <- node_ALL[train,]
  coldata.test <- node_ALL[test,]
  forDESeq <- cts.train
  DESeqColdata <- coldata.train
  
  dds <- DESeqDataSetFromMatrix(countData = forDESeq,
                                colData = DESeqColdata,
                                design = ~ colony + guarding)
  
  keep <- rowSums(counts(dds)) >= LOW.COUNT.FILTER
  dds <- dds[keep,]
  vsd <- vst(dds, blind = FALSE)
  DESeqNorm <- as.data.frame(assay(vsd))
  dds <- DESeq(dds)
  res <- na.omit(as.data.frame(results(dds)))
  ressig <- res[res$padj<= P,]
  DEGenes <- rownames(ressig)
  cts.train <- cts.train[rownames(cts.train)%in% DEGenes,]
  cts.test <- cts.test[rownames(cts.test)%in% DEGenes,]
  SVM.train <- rbind(IDs = cts.train, guarding = coldata.train$guarding)
  SVM.test <- rbind(IDs =cts.test, guarding = coldata.test$guarding)
  SVM.train <- as.data.frame(t(SVM.train))
  SVM.test <- as.data.frame(t(SVM.test))
  if(length(DEGenes > 0)){
    classifier = svm(formula = guarding ~ .,
                     data = SVM.train,
                     type = 'eps-regression',
                     kernel = 'linear')
    y_pred <- predict(classifier, newdata = SVM.test[-ncol(SVM.test)])
    cm <- as.data.frame(cbind(real = SVM.test[,ncol(SVM.test)], pred = y_pred))
    out <- 0
    if (max(cm$pred) > (mean(cm$pred) + sd(cm$pred)*2)) {
      out <- out + sum(cm$pred > (mean(cm$pred) + sd(cm$pred)*2))
      cm <- cm[-which(cm$pred > (mean(cm$pred) + sd(cm$pred)*2)),]
    }
    if (min(cm$pred) < (mean(cm$pred) - sd(cm$pred)*2)){
      out <- out + sum(cm$pred < (mean(cm$pred) - sd(cm$pred)*2))
      cm <- cm[-which(cm$pred < (mean(cm$pred) - sd(cm$pred)*2)),]
    }
    guardingO <- c(guardingO, out)
    plot(cm$pred, cm$real)
    guardingR <- c(guardingR, summary(lm(cm$real~ cm$pred))$r.squared)
  } else{
    guardingR <- c(guardingR, NA)
    guardingO <- c(guardingO, NA)
  }
}
write.csv(guardingR, "Predictionguarding.csv")
write.csv(guardingO, "Outlierguarding.csv")

queenattendingR <- c()
queenattendingO  <- c()
for (it in 1:ITERATIONS){
  print(it)
  train <- sample(1:ncol(cts_ALL), ncol(cts_ALL)/(100/SPLIT), replace = FALSE)
  test <- seq(from = 1, to = ncol(cts_ALL), by = 1)
  test <- test[!test%in%train]
  cts.train <- cts_ALL[,train]
  cts.test <- cts_ALL[,test]
  coldata.train <- node_ALL[train,]
  coldata.test <- node_ALL[test,]
  forDESeq <- cts.train
  DESeqColdata <- coldata.train
  
  dds <- DESeqDataSetFromMatrix(countData = forDESeq,
                                colData = DESeqColdata,
                                design = ~ colony + queenattending)
  
  keep <- rowSums(counts(dds)) >= LOW.COUNT.FILTER
  dds <- dds[keep,]
  vsd <- vst(dds, blind = FALSE)
  DESeqNorm <- as.data.frame(assay(vsd))
  dds <- DESeq(dds)
  res <- na.omit(as.data.frame(results(dds)))
  ressig <- res[res$padj<= P,]
  DEGenes <- rownames(ressig)
  cts.train <- cts.train[rownames(cts.train)%in% DEGenes,]
  cts.test <- cts.test[rownames(cts.test)%in% DEGenes,]
  SVM.train <- rbind(IDs = cts.train, queenattending = coldata.train$queenattending)
  SVM.test <- rbind(IDs =cts.test, queenattending = coldata.test$queenattending)
  SVM.train <- as.data.frame(t(SVM.train))
  SVM.test <- as.data.frame(t(SVM.test))
  if(length(DEGenes > 0)){
    classifier = svm(formula = queenattending ~ .,
                     data = SVM.train,
                     type = 'eps-regression',
                     kernel = 'linear')
    y_pred <- predict(classifier, newdata = SVM.test[-ncol(SVM.test)])
    cm <- as.data.frame(cbind(real = SVM.test[,ncol(SVM.test)], pred = y_pred))
    out <- 0
    if (max(cm$pred) > (mean(cm$pred) + sd(cm$pred)*2)) {
      out <- out + sum(cm$pred > (mean(cm$pred) + sd(cm$pred)*2))
      cm <- cm[-which(cm$pred > (mean(cm$pred) + sd(cm$pred)*2)),]
    }
    if (min(cm$pred) < (mean(cm$pred) - sd(cm$pred)*2)){
      out <- out + sum(cm$pred < (mean(cm$pred) - sd(cm$pred)*2))
      cm <- cm[-which(cm$pred < (mean(cm$pred) - sd(cm$pred)*2)),]
    }
    queenattendingO <- c(queenattendingO, out)
    plot(cm$pred, cm$real)
    queenattendingR <- c(queenattendingR, summary(lm(cm$real~ cm$pred))$r.squared)
  } else{
    queenattendingR <- c(queenattendingR, NA)
    queenattendingO <- c(queenattendingO, NA)
  }
}
write.csv(queenattendingR, "PredictionQieen.csv")
write.csv(queenattendingO, "OutlierQueen.csv")

taskR <- c()
taskO  <- c()
for (it in 1:ITERATIONS){
  print(it)
  train <- sample(1:ncol(cts_ALL), ncol(cts_ALL)/(100/SPLIT), replace = FALSE)
  test <- seq(from = 1, to = ncol(cts_ALL), by = 1)
  test <- test[!test%in%train]
  cts.train <- cts_ALL[,train]
  cts.test <- cts_ALL[,test]
  coldata.train <- node_ALL[train,]
  coldata.test <- node_ALL[test,]
  forDESeq <- cts.train
  DESeqColdata <- coldata.train
  
  dds <- DESeqDataSetFromMatrix(countData = forDESeq,
                                colData = DESeqColdata,
                                design = ~ colony + taskPCA_ALL)
  
  keep <- rowSums(counts(dds)) >= LOW.COUNT.FILTER
  dds <- dds[keep,]
  vsd <- vst(dds, blind = FALSE)
  DESeqNorm <- as.data.frame(assay(vsd))
  dds <- DESeq(dds)
  res <- na.omit(as.data.frame(results(dds)))
  ressig <- res[res$padj<= P,]
  DEGenes <- rownames(ressig)
  cts.train <- cts.train[rownames(cts.train)%in% DEGenes,]
  cts.test <- cts.test[rownames(cts.test)%in% DEGenes,]
  SVM.train <- rbind(cts.train, task = coldata.train$taskPCA_ALL)
  SVM.test <- rbind(cts.test, task = coldata.test$taskPCA_ALL)
  SVM.train <- as.data.frame(t(SVM.train))
  SVM.test <- as.data.frame(t(SVM.test))
  classifier = svm(formula = task ~ .,
                   data = SVM.train,
                   type = 'eps-regression',
                   kernel = 'linear')
  y_pred <- predict(classifier, newdata = SVM.test[-ncol(SVM.test)])
  cm <- as.data.frame(cbind(real = SVM.test[,ncol(SVM.test)], pred = y_pred))
  out <- 0
  if (max(cm$pred) > (mean(cm$pred) + sd(cm$pred)*2)) {
    out <- out + sum(cm$pred > (mean(cm$pred) + sd(cm$pred)*2))
    cm <- cm[-which(cm$pred > (mean(cm$pred) + sd(cm$pred)*2)),]
  }
  if (min(cm$pred) < (mean(cm$pred) - sd(cm$pred)*2)){
    out <- out + sum(cm$pred < (mean(cm$pred) - sd(cm$pred)*2))
    cm <- cm[-which(cm$pred < (mean(cm$pred) - sd(cm$pred)*2)),]
  }
  taskO <- c(taskO, out)
  plot(cm$pred, cm$real)
  taskR <- c(taskR, summary(lm(cm$real~ cm$pred))$r.squared)
}
write.csv(taskR, "PredictionTaskPCA.csv")
write.csv(taskO, "OutlierTaskPCA.csv")

# Read back in predictions
mat_to_plot <- read.csv("PredictionMaturity.csv", row.names = 1)
box_to_plot <- read.csv("PredictionProp.csv", row.names = 1)
nur_to_plot <- read.csv("PredictionNursing.csv", row.names = 1)
cle_to_plot <- read.csv("PredictionCleaning.csv", row.names = 1)
gua_to_plot <- read.csv("PredictionGuarding.csv", row.names = 1)
que_to_plot <- read.csv("PredictionQueen.csv", row.names = 1)
age_to_plot <- read.csv("PredictionAge.csv", row.names = 1)
beh_to_plot <- read.csv("PredictionTaskPCA.csv", row.names = 1)

timeline <- data.frame(Maturity = mat_to_plot$x, Behaviour_PC1 = beh_to_plot$x,
                       Age = age_to_plot$x, Foraging = box_to_plot$x,
                       Nursing = nur_to_plot$x,Queen_contact = que_to_plot$x,
                       Guarding = gua_to_plot$x, Cleaning = cle_to_plot$x)
labels <- names(timeline)
labels[2] <- "PC1 of behavior"
labels[6] <- "Queen interaction"
labels[1] <- "Social maturity"

setwd(FIGDATADIR)
write.csv(timeline, "Fig4A.csv", row.names = FALSE)

setwd(FIGDIR)
jpeg('Quantitative_Prediction_Comparison.jpg', width=3000, height=2000, unit='px')
par(oma=c(0,0,0,0), family="serif", mgp=c(15,4,0),  cex.axis=6, cex.lab=9, tcl=-0.2, lend=2, mar=c(7,5,4,1), mai=c(4,9.5,2,2))
boxplot(timeline, horizontal = TRUE, yaxt = "n", col = "gray", notch = TRUE, xlab = expression("Prediction vs. observation (R"^2*")"),
        lwd = 10)
axis(side = 2, las = 2, labels = labels, at = 1:length(names(timeline)))
dev.off()

# 10 random iterations with social maturity for pedagogical plot
it.cols <- c("cyan", "chartreuse", "blue", "darkmagenta", "firebrick2", "gray0",
             "gray", "darkgoldenrod1", "hotpink", "darkgreen")
valid_verify_BEHdf_full <- data.frame(real = NA, pred = NA, col = NA)
for (it in 1:10){
  print(it)
  train <- sample(1:ncol(cts_ALL), ncol(cts_ALL)/(100/SPLIT), replace = FALSE)
  test <- seq(from = 1, to = ncol(cts_ALL), by = 1)
  test <- test[!test%in%train]
  cts.train <- cts_ALL[,train]
  cts.test <- cts_ALL[,test]
  coldata.train <- node_ALL[train,]
  coldata.test <- node_ALL[test,]
  forDESeq <- cts.train
  DESeqColdata <- coldata.train
  dds <- DESeqDataSetFromMatrix(countData = forDESeq,
                                colData = DESeqColdata,
                                design = ~ colony + maturity)
  keep <- rowSums(counts(dds)) >= LOW.COUNT.FILTER
  dds <- dds[keep,]
  vsd <- vst(dds, blind = FALSE)
  DESeqNorm <- as.data.frame(assay(vsd))
  
  dds <- DESeq(dds)
  res <- na.omit(as.data.frame(results(dds)))
  ressig <- res[res$padj<= P,]
  DEGenes <- rownames(ressig)
  cts.train <- cts.train[rownames(cts.train)%in% DEGenes,]
  cts.test <- cts.test[rownames(cts.test)%in% DEGenes,]
  SVM.train <- rbind(cts.train, maturity = coldata.train$maturity)
  SVM.test <- rbind(cts.test, maturity = coldata.test$maturity)
  SVM.train <- as.data.frame(t(SVM.train))
  SVM.test <- as.data.frame(t(SVM.test))
  classifier = svm(formula = maturity ~ .,
                   data = SVM.train,
                   type = 'eps-regression',
                   kernel = 'linear')
  y_pred <- predict(classifier, newdata = SVM.test[-ncol(SVM.test)])
  cm <- as.data.frame(cbind(real = SVM.test[,ncol(SVM.test)], pred = y_pred))
  out <- 0
  if (max(cm$pred) > (mean(cm$pred) + sd(cm$pred)*2)) {
    out <- out + sum(cm$pred > (mean(cm$pred) + sd(cm$pred)*2))
    cm <- cm[-which(cm$pred > (mean(cm$pred) + sd(cm$pred)*2)),]
  }
  if (min(cm$pred) < (mean(cm$pred) - sd(cm$pred)*2)){
    out <- out + sum(cm$pred < (mean(cm$pred) - sd(cm$pred)*2))
    cm <- cm[-which(cm$pred < (mean(cm$pred) - sd(cm$pred)*2)),]
  }
  
  plot(cm$pred, cm$real)
  cm$col <- it.cols[it]
  
  valid_verify_BEHdf_full <- rbind(valid_verify_BEHdf_full, cm)
}

jpeg('MAT_Predictions.jpg', width=2000, height=2000, unit='px')
par(oma=c(0,0,0,0), family="serif", mgp=c(8,8,0),  cex.axis=1, cex.lab=1, tcl=-0.2, lend=2)
par(mfrow = c(1,1), mar=c(3,3,4,1), mai=c(5,5,2,2))
plot(valid_verify_BEHdf_full$pred ~ valid_verify_BEHdf_full$real, xlab = "", 
     ylab = "", pch = 19, cex = 10, col = as.character(valid_verify_BEHdf_full$col), cex.axis = 7,
     xlim = c(0,1), ylim = c(-.5,1.5), bty='l')
title(xlab="Actual position",mgp=c(18,1,0), cex.lab = 10)
title(ylab="Predicted position",mgp=c(18,1,0), cex.lab = 10)
abline(lm(valid_verify_BEHdf_full$pred ~ valid_verify_BEHdf_full$real), col = "darkgray", lwd = 20)
dev.off()
