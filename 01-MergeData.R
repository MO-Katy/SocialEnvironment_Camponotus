# SNG Network and Behavioral analysis
MAINDIR  <- "/Users/tkay/Desktop/Work/SNG/data" # Set data directory
FIGDIR   <- "/Users/tkay/Desktop/Work/SNG/figures"
FACETIN <- "/Users/tkay/Desktop/Work/SNG/FacetNet_Input"
FACETOUT <- "/Users/tkay/Desktop/Work/SNG/FacetNet_Output"
FACETNETDIR <- "/Users/tkay/Desktop/Work/facet_unil"

library(igraph); library(viridis);library(jpeg); library(assertthat);library(devtools);
library(ggfortify); library(pals); library(raster)

setwd(MAINDIR)
net_A  <- read.csv("EdgeData_A.csv")
node_A <- read.csv("NodeData_A.csv")
net_B  <- read.csv("EdgeData_B.csv")
node_B <- read.csv("NodeData_B.csv")
net_C  <- read.csv("EdgeData_C.csv")
node_C <- read.csv("NodeData_C.csv")

vid_A <- read.csv("SNG_vidAnnotation_A.csv")
vid_A[is.na(vid_A)] <- 0
vid_B <- read.csv("SNG_vidAnnotation_B.csv")
vid_B[is.na(vid_B)] <- 0
vid_C <- read.csv("SNG_vidAnnotation_C.csv")
vid_C[is.na(vid_C)] <- 0

age_A <- read.csv("TagsFile_A.csv")
age_B <- read.csv("TagsFile_B.csv")
age_C <- read.csv("TagsFile_C.csv")

md_A <- merge(node_A, age_A, by = "Tag")
md_A <- md_A[,!colnames(md_A) == "antID.y"]
colnames(md_A)[colnames(md_A) == "antID.x"] <- "antID"
md_B <- merge(node_B, age_B, by = "Tag")
md_B <- md_B[,!colnames(md_B) == "antID.y"]
colnames(md_B)[colnames(md_B) == "antID.x"] <- "antID"
md_C <- merge(node_C, age_C, by = "Tag")
md_C <- md_C[,!colnames(md_C) == "antID.y"]
colnames(md_C)[colnames(md_C) == "antID.x"] <- "antID"

md_A1 <- merge(md_A, vid_A, by = "antID")
md_B1 <- merge(md_B, vid_B, by = "antID")
md_C1 <- merge(md_C, vid_C, by = "antID")

to_exclude_A <- md_A1$antID[which(md_A1$count < mean(md_A1$count) - sd(md_A1$count)*2)]
to_exclude_B <- md_B1$antID[which(md_B1$count < mean(md_B1$count) - sd(md_B1$count)*2)] 
to_exclude_C <- md_C1$antID[which(md_C1$count < mean(md_C1$count) - sd(md_C1$count)*2)]

md_A1 <- md_A1[!(md_A1$antID %in% to_exclude_A),]
md_B1 <- md_B1[!(md_B1$antID %in% to_exclude_B),]
md_C1 <- md_C1[!(md_C1$antID %in% to_exclude_C),]

net_A1 <- net_A[net_A$ant1 %in% md_A1$antID & net_A$ant2 %in% md_A1$antID,]
net_B1 <- net_B[net_B$ant1 %in% md_B1$antID & net_B$ant2 %in% md_B1$antID,]
net_C1 <- net_C[net_C$ant1 %in% md_C1$antID & net_C$ant2 %in% md_C1$antID,]

net_A1 <- aggregate(net_A1$weightHH ~ net_A1$ant1 + net_A1$ant2, FUN = sum)
colnames(net_A1) <- c("ant1", "ant2", "weightHH")
net_B1 <- aggregate(net_B1$weightHH ~ net_B1$ant1 + net_B1$ant2, FUN = sum)
colnames(net_B1) <- c("ant1", "ant2", "weightHH")
net_C1 <- aggregate(net_C1$weightHH ~ net_C1$ant1 + net_C1$ant2, FUN = sum)
colnames(net_C1) <- c("ant1", "ant2", "weightHH")

# net_A1 <- aggregate(net_A1$weight ~ net_A1$ant1 + net_A1$ant2, FUN = sum)
# colnames(net_A1) <- c("ant1", "ant2", "weightHH")
# net_B1 <- aggregate(net_B1$weight ~ net_B1$ant1 + net_B1$ant2, FUN = sum)
# colnames(net_B1) <- c("ant1", "ant2", "weightHH")
# net_C1 <- aggregate(net_C1$weight ~ net_C1$ant1 + net_C1$ant2, FUN = sum)
# colnames(net_C1) <- c("ant1", "ant2", "weightHH")

node_O <- read.csv("NodeData_O.csv")
net_O  <- read.csv("EdgeData_O.csv")
net_O <- net_O[net_O$ant1 %in% node_O$antID & net_O$ant2 %in% node_O$antID,]
node_O$count <- node_O$nursebox + node_O$foragebox
to_exclude_O <- node_O$antID[which(node_O$count < mean(na.omit(node_O$count)) - sd(na.omit(node_O$count))*2)]
net_O  <- net_O[!(net_O$ant1 %in% to_exclude_O | net_O$ant2 %in% to_exclude_O),]


### Write data for FacetNET
setwd(paste(FACETIN, "contacts_A", sep = "/"))
write.table(net_A1, "0.edgelist", row.names = FALSE, col.names=FALSE)
setwd(paste(FACETIN, "contacts_B", sep = "/"))
write.table(net_B1, "0.edgelist", row.names = FALSE, col.names=FALSE)
setwd(paste(FACETIN, "contacts_C", sep = "/"))
write.table(net_C1, "0.edgelist", row.names = FALSE, col.names=FALSE)
setwd(paste(FACETIN, "contacts_O", sep = "/"))
write.table(net_O, "0.edgelist", row.names = FALSE, col.names=FALSE)

### run FacetNET
alpha <- 0.9
m <- 2
t_steps <- 1
command <- paste("python3", paste(FACETNETDIR, "facetnet_step.py", sep = "/"), paste(FACETIN,"contacts_A", "0.edgelist", sep = "/"), alpha, m, paste(FACETOUT, "contacts_A", sep = "/"), t_steps, sep = " ")
system(command)
command <- paste("python3", paste(FACETNETDIR, "facetnet_step.py", sep = "/"), paste(FACETIN,"contacts_B", "0.edgelist", sep = "/"), alpha, m, paste(FACETOUT, "contacts_B", sep = "/"), t_steps, sep = " ")
system(command)
command <- paste("python3", paste(FACETNETDIR, "facetnet_step.py", sep = "/"), paste(FACETIN,"contacts_C", "0.edgelist", sep = "/"), alpha, m, paste(FACETOUT, "contacts_C", sep = "/"), t_steps, sep = " ")
system(command)
command <- paste("python3", paste(FACETNETDIR, "facetnet_step.py", sep = "/"), paste(FACETIN,"contacts_O", "0.edgelist", sep = "/"), alpha, m, paste(FACETOUT, "contacts_O", sep = "/"), t_steps, sep = " ")
system(command)

setwd(FACETOUT)
for (folder in list.files()){
  setwd(paste(FACETOUT, folder, sep = "/"))
  assign(sub("\\..*", "", paste(folder, "SoftMod", sep = "_")), read.csv(list.files()[which(grepl("soft_comm_step", list.files(), fixed = TRUE))]))
}

# Merge to metadata
colnames(contacts_A_SoftMod)[1] <- "antID"
colnames(contacts_B_SoftMod)[1] <- "antID"
colnames(contacts_C_SoftMod)[1] <- "antID"

md_A2 <- merge(md_A1, contacts_A_SoftMod, by = "antID")
md_B2 <- merge(md_B1, contacts_B_SoftMod, by = "antID")
md_C2 <- merge(md_C1, contacts_C_SoftMod, by = "antID")
md_A2$prop <- md_A2$forage / (md_A2$nest + md_A2$forage)
md_B2$prop <- md_B2$forage / (md_B2$nest + md_B2$forage)
md_C2$prop <- md_C2$forage / (md_C2$nest + md_C2$forage)
md_A2$trophallaxis <- md_A2$trophallaxis + md_A2$trophallaxisout
md_A2$cleaning <- md_A2$cleaning + md_A2$cleaningout
md_A2$grooming <- md_A2$grooming + md_A2$groomingout
md_A2$out <- md_A2$trophallaxisout + md_A2$cleaningout + md_A2$groomingout
md_B2$trophallaxis <- md_B2$trophallaxis + md_B2$trophallaxisout
md_B2$grooming <- md_B2$grooming + md_B2$groomingout
md_B2$out <- md_B2$trophallaxisout + md_B2$cleaningout + md_B2$groomingout
md_C2$trophallaxis <- md_C2$trophallaxis + md_C2$trophallaxisout
md_C2$grooming <- md_C2$grooming + md_C2$groomingout
md_C2$out <- md_C2$trophallaxisout + md_C2$cleaningout + md_C2$groomingout

net_A2 <- net_A1[net_A1$ant1 %in% md_A2$antID & net_A1$ant2 %in% md_A2$antID,]
net_B2 <- net_B1[net_B1$ant1 %in% md_B2$antID & net_B1$ant2 %in% md_B2$antID,]
net_C2 <- net_C1[net_C1$ant1 %in% md_C2$antID & net_C1$ant2 %in% md_C2$antID,]

# set social maturity
if (cor.test(md_A2$cluster_0, md_A2$forage)$estimate > 0){
  names(md_A2)[which(names(md_A2) == "cluster_0")] <- "maturity"
} else {
  names(md_A2)[which(names(md_A2) == "cluster_1")] <- "maturity"
}
if (cor.test(md_B2$cluster_0, md_B2$forage)$estimate > 0){
  names(md_B2)[which(names(md_B2) == "cluster_0")] <- "maturity"
} else {
  names(md_B2)[which(names(md_B2) == "cluster_1")] <- "maturity"
}
if (cor.test(md_C2$cluster_0, md_C2$forage)$estimate > 0){
  names(md_C2)[which(names(md_C2) == "cluster_0")] <- "maturity"
} else {
  names(md_C2)[which(names(md_C2) == "cluster_1")] <- "maturity"
}

setwd(MAINDIR)
write.csv(net_A2, "Net_A2.csv", row.names = F)
write.csv(net_B2, "Net_B2.csv", row.names = F)
write.csv(net_C2, "Net_C2.csv", row.names = F)
write.csv(md_A2, "md_A2.csv", row.names = F)
write.csv(md_B2, "md_B2.csv", row.names = F)
write.csv(md_C2, "md_C2.csv", row.names = F)

