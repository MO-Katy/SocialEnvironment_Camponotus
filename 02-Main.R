# SNG Network and Behavioral analysis
MAINDIR  <- "/Users/tkay/Desktop/Work/SNG/data"
FIGDIR   <- "/Users/tkay/Desktop/Work/SNG/figures"

library(igraph); library(viridis);library(jpeg); library(assertthat);library(devtools);
library(ggfortify); library(pals); library(lmerTest); library(r2glmm);library(stringr); 
library(raster);

setwd(MAINDIR)
net_A  <- read.csv("Net_A2.csv")
net_B  <- read.csv("Net_B2.csv")
net_C  <- read.csv("Net_C2.csv")
node_A <- read.csv("md_A2.csv")
node_B <- read.csv("md_B2.csv")
node_C <- read.csv("md_C2.csv")
node_O <- read.csv("NodeData_O.csv")
net_O  <- read.csv("EdgeData_O.csv")
net_O <- net_O[net_O$ant1 %in% node_O$antID & net_O$ant2 %in% node_O$antID,]
Nest_A <- read.csv("normalized_counts_space_sng1_Nest_space_2021-11-02 00_01_00+01_00_2021-11-08 23_59_00+01_00.csv")
Fora_A <- read.csv("normalized_counts_space_sng1_Forage_space_2021-11-02 00_01_00+01_00_2021-11-08 23_59_00+01_00.csv")
Nest_B <- read.csv("normalized_counts_space_SNG_B_Nest_space_2021-11-25 00_01_00+01_00_2021-12-01 23_59_00+01_00.csv")
Fora_B <- read.csv("normalized_counts_space_SNG_B_Forage_space_2021-11-25 00_01_00+01_00_2021-12-01 23_59_00+01_00.csv")
Nest_C <- read.csv("normalized_counts_space_SNG colony C_Nest_space_2021-12-10 00_01_00+01_00_2021-12-16 23_59_00+01_00.csv")
Fora_C <- read.csv("normalized_counts_space_SNG colony C_Forage_space_2021-12-10 00_01_00+01_00_2021-12-16 23_59_00+01_00.csv")
space_O <- read.csv("Space_O.csv")
MiBi_O   <- read.csv("ASVtable_O.csv", row.names=1)
MiBi_ABC <- read.csv("ASVtable_ABC.csv", row.names = 1)

# Exclusion of lowly detected individuals from colony 1
node_O$count <- node_O$nursebox + node_O$foragebox
to_exclude_O <- node_O$antID[which(node_O$count < mean(na.omit(node_O$count)) - sd(na.omit(node_O$count))*2)]
net_O  <- net_O[!(net_O$ant1 %in% to_exclude_O | net_O$ant2 %in% to_exclude_O),]
node_O <- node_O[!(node_O$antID %in% to_exclude_O),]
space_O <- space_O[!(space_O$Ant %in% to_exclude_O),]
MiBi_O <- MiBi_O[!(rownames(MiBi_O) %in% to_exclude_O),]

# Age distributions
node_B$Age <- as.numeric(node_B$Age)
node_C$Age <- as.numeric(node_C$Age)

# Plot age distribution
AGE <- c("0-5", "5-10",  "10-15", "15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50")
FREQ <- hist(c(node_A$Age, node_B$Age, node_C$Age, na.omit(node_O$age)))$counts
AGEqplot <- data.frame(bin = AGE, count = FREQ)
AGEqplot$bin <- factor(AGEqplot$bin, levels = c("0-5", "5-10",  "10-15", "15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50"))

setwd(FIGDIR)
jpeg('Age_Distribution.jpg', width=4000, height=2000, unit='px')
ggplot(data=AGEqplot, aes(x=bin, y=count)) +
  ylim(0,130) +
  geom_bar(stat="identity", width = c(rep(1,10)), col = "darkgray", fill = "darkgray") + xlab("\nAge (weeks)") + ylab("Frequency\n") + 
  theme(axis.text=element_text(size=110,family="serif"), axis.title = element_text(size=170,family="serif"),
        panel.background = element_blank())
dev.off()
jpeg('Age_Distribution_PerCol.jpg', width=4000, height=2000, unit='px')
par(mfrow = c(1,4), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
hist(node_O$age, yaxt="n", xaxt="n", cex.lab = 10, cex.main = 15, cex = 10,
     xlab = "Age", main = "Colony 1", xlim = c(0,55), ylim = c(0,45), col = "black")
axis(2, at = c(0,45), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,55), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
hist(node_A$Age, yaxt="n", xaxt="n", cex.lab = 10, cex.main = 15, cex = 10,
     xlab = "Age", main = "Colony 2", xlim = c(0,55), ylim = c(0,45), col = "black")
axis(2, at = c(0,45), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,55), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
hist(node_B$Age, yaxt="n", xaxt="n", cex.lab = 10, cex.main = 15, cex = 10,
     xlab = "Age", main = "Colony 3", xlim = c(0,55), ylim = c(0,45), col = "black")
axis(2, at = c(0,45), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,55), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
hist(node_C$Age, yaxt="n", xaxt="n", cex.lab = 10, cex.main = 15, cex = 10,
     xlab = "Age", main = "Colony 4", xlim = c(0,55), ylim = c(0,45), col = "black")
axis(2, at = c(0,45), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,55), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
dev.off()

# Construct iGraph object
net_A$ant1 <- as.character(net_A$ant1)
net_A$ant2 <- as.character(net_A$ant2)
net_A1 <- graph_from_edgelist(as.matrix(net_A[,1:2]), directed = FALSE)
E(net_A1)$weight <- net_A[,3]
E(net_A1)$weight <- round((E(net_A1)$weight/max(E(net_A1)$weight))*1000) + 1
net_layout_A1 <- layout_with_fr(net_A1)
E(net_A1)$color <- gray(1-(E(net_A1)$weight/max(E(net_A1)$weight)))
net_shape_A1 <- rep("circle", length(V(net_A1)))
net_shape_A1[which(node_A$queen == 1)] <- "square"
V(net_A1)$shape <- net_shape_A1
net_color_A1 <- c(parula(116)[rank((node_A$nest / (node_A$nest + node_A$forage))[-117])], "magenta")
net_color_A2 <- c(parula(116)[rank(1-node_A$maturity[-116])], "magenta")
net_color_A3 <- c(parula(116)[rank(1/node_A$Age[-117])], "magenta")

net_B$ant1 <- as.character(net_B$ant1)
net_B$ant2 <- as.character(net_B$ant2)
net_B1 <- graph_from_edgelist(as.matrix(net_B[,1:2]), directed = FALSE)
E(net_B1)$weight <- net_B[,3]
E(net_B1)$weight <- round((E(net_B1)$weight/max(E(net_B1)$weight))*1000) + 1
net_layout_B1 <- layout_with_fr(net_B1)
E(net_B1)$color <- gray(1-(E(net_B1)$weight/max(E(net_B1)$weight)))
net_shape_B1 <- rep("circle", length(V(net_B1)))
net_shape_B1[which(node_B$queen == 1)] <- "square"
V(net_B1)$shape <- net_shape_B1
net_color_B1 <- c(parula(126)[rank((node_B$nest / (node_B$nest + node_B$forage))[-127])], "magenta")
net_color_B2 <- c(parula(126)[rank(1-node_B$maturity[-127])], "magenta")
net_color_B3 <- node_B$Age[-127]
net_color_B3 <- c(parula(126)[rank(1/as.numeric(net_color_B3))], "magenta")

net_C$ant1 <- as.character(net_C$ant1)
net_C$ant2 <- as.character(net_C$ant2)
net_C1 <- graph_from_edgelist(as.matrix(net_C[,1:2]), directed = FALSE)
E(net_C1)$weight <- net_C[,3]

# Remove biggest outlier for plotting
E(net_C1)$weight[which(E(net_C1)$weight == max(E(net_C1)$weight))] <- round(mean(E(net_C1)$weight))
E(net_C1)$weight <- round((E(net_C1)$weight/max(E(net_C1)$weight))*1000)  + 1

net_layout_C1 <- layout_with_fr(net_C1)
E(net_C1)$color <- gray(1-(E(net_C1)$weight/max(E(net_C1)$weight)))
net_shape_C1 <- rep("circle", length(V(net_C1)))
net_shape_C1[which(node_C$queen == 1)] <- "square"
V(net_C1)$shape <- net_shape_C1
net_color_C1 <- c(parula(130)[rank((node_C$nest / (node_C$nest + node_C$forage))[-121])], "magenta")
net_color_C2 <- c(parula(130)[rank(1-node_C$maturity[-131])], "magenta")
net_color_C3 <- node_C$Age[-131]
net_color_C3 <- c(parula(130)[rank(1/as.numeric(net_color_C3))], "magenta")

net_O$ant1 <- as.character(net_O$ant1)
net_O$ant2 <- as.character(net_O$ant2)
net_O$weightHH <- round(net_O$weightHH/3,0) # For plotting only
net_Ox <- data.frame(ant1 = net_O$ant2, ant2 = net_O$ant1)
net_O1 <- graph_from_edgelist(as.matrix(net_Ox[,1:2]), directed = FALSE)
E(net_O1)$weight <- net_O[,3]
E(net_O1)$weight <- round((E(net_O1)$weight/max(E(net_O1)$weight))*1000)  + 1
net_layout_O1 <- layout_with_fr(net_O1)
E(net_O1)$color <- gray(1-(E(net_O1)$weight/max(E(net_O1)$weight)))
net_shape_O1 <- rep("circle", length(V(net_O1)))
net_shape_O1[which(node_O$queen == 1)] <- "square"
V(net_O1)$shape <- net_shape_O1
net_color_O1 <- c(parula(80)[rank(1-node_O$boxratio)[-81]], "magenta")
net_color_O2 <- c(parula(80)[rank(1-node_O$soft)[-81]], "magenta")
net_color_O3 <- node_O$age
net_color_O3 <- c(parula(80)[rank(1/as.numeric(net_color_O3))[-81]], "magenta")
net_color_O3[is.na(node_O$age)] <-"gray"
net_color_O3[length(net_color_O3)] <-"magenta"

##########################
# Visual network inspection
##########################

setwd(FIGDIR)
jpeg('Network_A.jpg', width=4000, height=2000, unit='px')
par(mfrow = c(1,3), mar = c(0,0,0,0), oma = c(0,0,0,0), bty="n", mgp = c(0,0,0), mai = c(0,0,0,0), family = "serif")
plot(net_A1, layout = net_layout_A1,
     vertex.size=7, vertex.label=NA, vertex.color = net_color_A1, # Vertex size = 7
     edge.width = ((E(net_A1)$weight/max(E(net_A1)$weight)))*40)
text(0.15,1.25,"Time foraging", cex = 20)
plot(net_A1, layout = net_layout_A1,
     vertex.size=7, vertex.label=NA, vertex.color = net_color_A2,
     edge.width = ((E(net_A1)$weight/max(E(net_A1)$weight)))*40)
text(0.15,1.25,"Social maturity", cex = 20)
plot(net_A1, layout = net_layout_A1,
     vertex.size=7, vertex.label=NA, vertex.color = net_color_A3,
     edge.width = ((E(net_A1)$weight/max(E(net_A1)$weight)))*40)
text(0.15,1.25,"Age", cex = 20)
dev.off()

jpeg('Network_B.jpg', width=4000, height=2000, unit='px')
par(mfrow = c(1,3), mar = c(0,0,0,0), oma = c(0,0,0,0), bty="n", mgp = c(0,0,0), mai = c(0,0,0,0), family = "serif")
plot(net_B1, layout = net_layout_B1,
     vertex.size=7, vertex.label=NA, vertex.color = net_color_B1,
     edge.width = ((E(net_B1)$weight/max(E(net_B1)$weight)))*40)
text(0.15,1.25,"Time foraging", cex = 20)
plot(net_B1, layout = net_layout_B1,
     vertex.size=7, vertex.label=NA, vertex.color = net_color_B2,
     edge.width = ((E(net_B1)$weight/max(E(net_B1)$weight)))*40)
text(0.15,1.25,"Social maturity", cex = 20)
plot(net_B1, layout = net_layout_B1,
     vertex.size=7, vertex.label=NA, vertex.color = net_color_B3,
     edge.width = ((E(net_B1)$weight/max(E(net_B1)$weight)))*40)
text(0.15,1.25,"Age", cex = 20)
dev.off()

jpeg('Network_C.jpg', width=4000, height=2000, unit='px')
par(mfrow = c(1,3), mar = c(0,0,0,0), oma = c(0,0,0,0), bty="n", mgp = c(0,0,0), mai = c(0,0,0,0), family = "serif")
plot(net_C1, layout = net_layout_C1,
     vertex.size=7, vertex.label=NA, vertex.color = net_color_C1,
     edge.width = ((E(net_C1)$weight/max(E(net_C1)$weight)))*40)
text(0.15,1.25,"Time foraging", cex = 20)
plot(net_C1, layout = net_layout_C1,
     vertex.size=7, vertex.label=NA, vertex.color = net_color_C2,
     edge.width = ((E(net_C1)$weight/max(E(net_C1)$weight)))*40)
text(0.15,1.25,"Social maturity", cex = 20)
plot(net_C1, layout = net_layout_C1,
     vertex.size=7, vertex.label=NA, vertex.color = net_color_C3,
     edge.width = ((E(net_C1)$weight/max(E(net_C1)$weight)))*40)
text(0.15,1.25,"Age", cex = 20)
dev.off()

jpeg('Network_O.jpg', width=4000, height=2000, unit='px')
par(mfrow = c(1,3), mar = c(0,0,0,0), oma = c(0,0,0,0), bty="n", mgp = c(0,0,0), mai = c(0,0,0,0), family = "serif")
plot(net_O1, layout = net_layout_O1,
     vertex.size=7, vertex.label=NA, vertex.color = net_color_O1,
     edge.width = ((E(net_O1)$weight/max(E(net_O1)$weight)))*40)
text(0.15,1.25,"Time foraging", cex = 20)
plot(net_O1, layout = net_layout_O1,
     vertex.size=7, vertex.label=NA, vertex.color = net_color_O2,
     edge.width = ((E(net_O1)$weight/max(E(net_O1)$weight)))*40)
text(0.15,1.25,"Social maturity", cex = 20)
plot(net_O1, layout = net_layout_O1,
     vertex.size=7, vertex.label=NA, vertex.color = net_color_O3,
     edge.width = ((E(net_O1)$weight/max(E(net_O1)$weight)))*40)
text(0.15,1.25,"Age", cex = 20)
dev.off()

##########################
# Social Maturity distributions
##########################
# remove queens
node_A <- node_A[-nrow(node_A),]
node_B <- node_B[-nrow(node_B),]
node_C <- node_C[-nrow(node_C),]
node_O <- node_O[-nrow(node_O),]

jpeg('Maturity_Distribution.jpg', width=5000, height=2000, unit='px')
par(mfrow = c(1,4), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
hist(node_A$maturity,
     xlab = "Social Maturity", col = "black", main = "",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 15, cex = 5)
axis(2, at = c(0,140), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
hist(node_B$maturity,
     xlab = "Social Maturity", col = "black", main = "",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 15, cex = 5)
axis(2, at = c(0,140), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
hist(node_C$maturity,
     xlab = "Social Maturity", col = "black", main = "",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 15, cex = 5)
axis(2, at = c(0,140), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
hist(node_O$soft,
     xlab = "Social Maturity", col = "black", main = "",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 15, cex = 5)
axis(2, at = c(0,140), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
dev.off()

##########################
# Social Maturity v Foraging v Age
##########################

Forage_ALL <- data.frame(Age = c(node_A$Age, node_B$Age, node_C$Age, node_O$age),
                         Forage = c(node_A$prop, node_B$prop, node_C$prop, node_O$boxratio),
                         Nurse = c(node_A$broodcare, node_B$broodcare, node_C$broodcare, node_O$nursing),
                         Maturity = c(node_A$maturity, node_B$maturity, node_C$maturity, node_O$soft),
                         colony = c(rep("A", nrow(node_A)),rep("B", nrow(node_B)),rep("C", nrow(node_C)),rep("O", nrow(node_O))))

r2beta(lmerTest::lmer(Age ~ Maturity + (1|colony), Forage_ALL))
summary(lmerTest::lmer(Age ~ Maturity + (1|colony), Forage_ALL))
r2beta(lmerTest::lmer(Forage ~ Maturity + (1|colony), Forage_ALL))
summary(lmerTest::lmer(Forage ~ Maturity + (1|colony), Forage_ALL))


Forage_ALL$Nurse <- (Forage_ALL$Nurse - min(Forage_ALL$Nurse, na.rm = 1)) / (max(Forage_ALL$Nurse, na.rm = 1) - min(Forage_ALL$Nurse, na.rm = 1))
jpeg('ForMatAge.jpg', width=5000, height=2000, unit='px')
par(mfrow = c(1,3), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
plot(Forage_ALL$Forage~Forage_ALL$Age,
     #col = c(rep("orange", nrow(node_A)), rep("purple", nrow(node_B)), rep("lightblue", nrow(node_C)), rep("green3", nrow(node_O))),
     col = "darkblue",
     pch = 16,
     xlab = "Age (weeks)", ylab = "Prop. of time spent foraging",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 15, cex = 10)
axis(2, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,45), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
plot(Forage_ALL$Age~Forage_ALL$Maturity,
     #col = c(rep("orange", nrow(node_A)), rep("purple", nrow(node_B)), rep("lightblue", nrow(node_C)), rep("green3", nrow(node_O))),
     col = "darkblue",
     pch = 16,
     xlab = "Social Maturity", ylab = "Age (weeks)",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 15, cex = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(2, at = c(0,45), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
plot(Forage_ALL$Forage~Forage_ALL$Maturity,
     #col = c(rep("orange", nrow(node_A)), rep("purple", nrow(node_B)), rep("lightblue", nrow(node_C)), rep("green3", nrow(node_O))),
     col = "darkblue",
     pch = 16,
     xlab = "Social Maturity", ylab = "Prop. of time spent foraging",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 15, cex = 10)
axis(2, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
dev.off()

jpeg('ForNur.jpg', width=5000, height=2000, unit='px')
par(mfrow = c(1,2), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
plot(Forage_ALL$Forage~Forage_ALL$Age,
     col = "darkblue",
     pch = 16,
     xlab = "Age (weeks)", ylab = "Prop. time spent foraging",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 15, cex = 10)
abline(lm(Forage_ALL$Forage~Forage_ALL$Age), lwd = 25, col = "darkgray")
axis(2, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,45), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
plot(Forage_ALL$Nurse~Forage_ALL$Age,
     col = "darkblue",
     pch = 16,
     xlab = "Age (weeks)", ylab = "Prop. time spent nursing",
     yaxt="n", xaxt="n", cex.lab = 10, cex.main = 15, cex = 10)
abline(lm(Forage_ALL$Nurse~Forage_ALL$Age), lwd = 25, col = "darkgray")
axis(2, at = c(0,1), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
axis(1, at = c(0,45), labels = TRUE, tick = TRUE, lwd = 8, cex.axis = 10)
dev.off()

### Analysis of annotated behaviors
tasks <- c("cleaning", "trophallaxis", "broodcare",
           "guarding", "queenattending",
           "prop")
task_A <- node_A[,colnames(node_A) %in% tasks]
task_B <- node_B[,colnames(node_B) %in% tasks]
task_C <- node_C[,colnames(node_C) %in% tasks]
task_O <- node_O[,colnames(node_O) %in% c("cleaning", "tropholaxis", "nursing", "Prop_time_guarding_2BL_FPT1000sec",
                                          "N_Q_contacts", "boxratio")]

# In one colony one individual is missing task performance data. Replace with population average.
task_O[rownames(task_O)=="74",1:4] <- as.numeric(colMeans(na.omit(task_O))[1:4])

# Uniformize old colony with new colony
task_O <- task_O[,c(3,2,1,6,5,4)]
colnames(task_O) <- tasks

for (col in 1:ncol(task_A)){
  task_A[,col] <- (task_A[,col] - min(task_A[,col])) / ((max(task_A[,col]) - min(task_A[,col])))
}
for (col in 1:ncol(task_B)){
  task_B[,col] <- (task_B[,col] - min(task_B[,col])) / ((max(task_B[,col]) - min(task_B[,col])))
}
for (col in 1:ncol(task_C)){
  task_C[,col] <- (task_C[,col] - min(task_C[,col])) / ((max(task_C[,col]) - min(task_C[,col])))
}
for (col in 1:ncol(task_O)){
  task_O[,col] <- (task_O[,col] - min(task_O[,col])) / ((max(task_O[,col]) - min(task_O[,col])))
}

rownames(task_A) <- paste0("A_", rownames(task_A))
rownames(task_B) <- paste0("B_", rownames(task_B))
rownames(task_C) <- paste0("C_", rownames(task_C))
rownames(task_O) <- paste0("O_", rownames(task_O))

task_PCA_A <- prcomp(task_A)
task_PCA_B <- prcomp(task_B)
task_PCA_C <- prcomp(task_C)
task_PCA_O <- prcomp(task_O)

task_ALL <- rbind(task_A, task_B, task_C, task_O)

jpeg('Task_PCA_A.jpg', width=2000, height=2000, unit='px')
par(mfrow = c(1,1), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
autoplot(task_PCA_A,
         loadings = TRUE, loadings.colour = 'darkgray', loadings.label.colour = "darkgreen",
         loadings.size = 30, loadings.label = TRUE, loadings.label.size = 15) + 
  geom_point(size = 10) + 
  theme(text = element_text(size=80),
        axis.text = element_text(size=0))
dev.off()

jpeg('Task_PCA_B.jpg', width=2000, height=2000, unit='px')
par(mfrow = c(1,1), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
autoplot(task_PCA_B,
         loadings = TRUE, loadings.colour = 'darkgray', loadings.label.colour = "darkgreen",
         loadings.size = 30, loadings.label = TRUE, loadings.label.size = 15) + 
  geom_point(size = 10) + 
  theme(text = element_text(size=80),
        axis.text = element_text(size=0))
dev.off()

jpeg('Task_PCA_C.jpg', width=2000, height=2000, unit='px')
par(mfrow = c(1,1), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
autoplot(task_PCA_C,
         loadings = TRUE, loadings.colour = 'darkgray', loadings.label.colour = "darkgreen",
         loadings.size = 30, loadings.label = TRUE, loadings.label.size = 15) + 
  geom_point(size = 10) + 
  theme(text = element_text(size=80),
        axis.text = element_text(size=0))
dev.off()

jpeg('Task_PCA_O.jpg', width=2000, height=2000, unit='px')
par(mfrow = c(1,1), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
autoplot(task_PCA_O,
         loadings = TRUE, loadings.colour = 'darkgray', loadings.label.colour = "darkgreen",
         loadings.size = 30, loadings.label = TRUE, loadings.label.size = 15) + 
  geom_point(size = 10) + 
  theme(text = element_text(size=80),
        axis.text = element_text(size=0))
dev.off()

task_ALL$colony <- c(rep("orange", nrow(task_A)), rep("purple", nrow(task_B)), rep("lightblue", nrow(task_C)), rep("green3", nrow(task_O)))
task_ALL$age <- c(node_A$Age, node_B$Age, node_C$Age, node_O$age)
task_ALL$age[is.na(task_ALL$age)] <- mean(na.omit(task_ALL$age))
task_ALL$SM <- c(node_A$maturity, node_B$maturity, node_C$maturity, node_O$soft)

task_PCA_ALL <- prcomp(task_ALL[,1:6])

jpeg('Task_PCA_ALL.jpg', width=2000, height=2000, unit='px')
par(mfrow = c(1,1), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
autoplot(task_PCA_ALL,
         loadings = TRUE, loadings.colour = 'darkgray', loadings.label.colour = "darkgreen",
         loadings.size = 50, loadings.label = TRUE, loadings.label.size = 30) + 
  geom_point(size = 10, col = "darkblue") + 
  theme(text = element_text(size=80),
        axis.text = element_text(size=0))
dev.off()
jpeg('Task_PCA_ALL_age.jpg', width=2000, height=2000, unit='px')
par(mfrow = c(1,2), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
autoplot(task_PCA_ALL,
         loadings = TRUE, loadings.colour = 'darkgray', loadings.label.colour = "darkgreen",
         loadings.size = 50, loadings.label = TRUE, loadings.label.size = 30) + 
  geom_point(size = (task_ALL$age^2)/100, col = "darkblue") + 
  theme(text = element_text(size=80),
        axis.text = element_text(size=0))
dev.off()
jpeg('Task_PCA_ALL_sm.jpg', width=2000, height=2000, unit='px')
par(mfrow = c(1,2), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
autoplot(task_PCA_ALL,
         loadings = TRUE, loadings.colour = 'darkgray', loadings.label.colour = "darkgreen",
         loadings.size = 50, loadings.label = TRUE, loadings.label.size = 30) + 
  geom_point(size = 20, col = parula(452)[rank(1/task_ALL$SM)]) + 
  theme(text = element_text(size=80),
        axis.text = element_text(size=0))
dev.off()

summary(lm(task_PCA_A$x[,1]~ node_A$maturity))
summary(lm(task_PCA_B$x[,1]~ node_B$maturity))
summary(lm(task_PCA_C$x[,1]~ node_C$maturity))
summary(lm(task_PCA_O$x[,1]~ node_O$soft))

summary(lm(task_PCA_A$x[,1]~ node_A$Age))
summary(lm(task_PCA_B$x[,1]~ node_B$Age))
summary(lm(task_PCA_C$x[,1]~ node_C$Age))
summary(lm(task_PCA_O$x[,1]~ node_O$age))

summary(lm(node_A$maturity ~ node_A$Age))
summary(lm(node_B$maturity ~ node_B$Age))
summary(lm(node_C$maturity ~ node_C$Age))
summary(lm(node_O$soft ~ node_O$age))

# Network by behavioral caste
task_A$classification <- NA
for (ant in 1:nrow(task_A)){
  if(length(colnames(task_A)[which(task_A[ant,1:6] == max(task_A[ant,1:6]))]) == 1){
    task_A$classification[ant] <- colnames(task_A)[which(task_A[ant,1:6] == max(task_A[ant,1:6]))]
  }
}
task_A_col <- c(task_A[,7], "magenta")
task_A_col[task_A_col ==  "trophallaxis"] <- "darkgray"
task_A_col[task_A_col ==  "cleaning"] <- "yellow"
task_A_col[task_A_col ==  "broodcare"] <- "orange"
task_A_col[task_A_col ==  "guarding"] <- "red"
task_A_col[task_A_col ==  "prop"] <- "blue"
task_A_col[task_A_col ==  "queenattending"] <- "purple"

task_B$classification <- NA
for (ant in 1:nrow(task_B)){
  if(length(colnames(task_B)[which(task_B[ant,1:6] == max(task_B[ant,1:6]))]) == 1){
    task_B$classification[ant] <- colnames(task_B)[which(task_B[ant,1:6] == max(task_B[ant,1:6]))]
  }
}
task_B_col <- c(task_B[,7], "magenta")
task_B_col[task_B_col ==  "trophallaxis"] <- "darkgray"
task_B_col[task_B_col ==  "cleaning"] <- "yellow"
task_B_col[task_B_col ==  "broodcare"] <- "orange"
task_B_col[task_B_col ==  "guarding"] <- "red"
task_B_col[task_B_col ==  "prop"] <- "blue"
task_B_col[task_B_col ==  "queenattending"] <- "purple"

task_C$classification <- NA
for (ant in 1:nrow(task_C)){
  if(length(colnames(task_C)[which(task_C[ant,1:6] == max(task_C[ant,1:6]))]) == 1){
    task_C$classification[ant] <- colnames(task_C)[which(task_C[ant,1:6] == max(task_C[ant,1:6]))]
  }
}
task_C_col <- c(task_C[,7], "magenta")
task_C_col[task_C_col ==  "trophallaxis"] <- "darkgray"
task_C_col[task_C_col ==  "cleaning"] <- "yellow"
task_C_col[task_C_col ==  "broodcare"] <- "orange"
task_C_col[task_C_col ==  "guarding"] <- "red"
task_C_col[task_C_col ==  "prop"] <- "blue"
task_C_col[task_C_col ==  "queenattending"] <- "purple"

task_O$classification <- NA
for (ant in 1:nrow(task_O)){
  if(length(colnames(task_O)[which(task_O[ant,1:6] == max(task_O[ant,1:6]))]) == 1){
    task_O$classification[ant] <- colnames(task_O)[which(task_O[ant,1:6] == max(task_O[ant,1:6]))]
  }
}
task_O_col <- c(task_O[,7], "magenta")
task_O_col[task_O_col ==  "trophallaxis"] <- "darkgray"
task_O_col[task_O_col ==  "cleaning"] <- "yellow"
task_O_col[task_O_col ==  "broodcare"] <- "orange"
task_O_col[task_O_col ==  "guarding"] <- "red"
task_O_col[task_O_col ==  "prop"] <- "blue"
task_O_col[task_O_col ==  "queenattending"] <- "purple"

jpeg('Network_task.jpg', width=4000, height=2000, unit='px')
par(mfrow = c(1,5), mar = c(0,0,0,0), oma = c(0,0,0,0), bty="n", mgp = c(0,0,0), mai = c(0,0,0,0), family = "serif")
plot(net_A1, layout = net_layout_A1,
     vertex.size=7, vertex.label=NA, vertex.color = task_A_col,
     edge.width = ((E(net_A1)$weight/max(E(net_A1)$weight)))*40)
plot(net_B1, layout = net_layout_B1,
     vertex.size=7, vertex.label=NA, vertex.color = task_B_col,
     edge.width = ((E(net_B1)$weight/max(E(net_B1)$weight)))*40)
plot(net_C1, layout = net_layout_C1,
     vertex.size=7, vertex.label=NA, vertex.color = task_C_col,
     edge.width = ((E(net_C1)$weight/max(E(net_C1)$weight)))*40)
plot(net_O1, layout = net_layout_O1,
     vertex.size=7, vertex.label=NA, vertex.color = task_O_col,
     edge.width = ((E(net_O1)$weight/max(E(net_O1)$weight)))*10)
plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10))
legend(1, 8, legend=c("queen", "trophallaxis", "cleaning", "broodcare", "guarding", "foraging", "queen attending"),
              fill=c("magenta", "darkgray", "yellow", "orange", "red", "blue", "purple"), lty=1:2, cex=8 , bty = "n")
dev.off()

# Spatial similarity
space_A <- cbind(Nest_A, Fora_A[,-1])
space_B <- cbind(Nest_B, Fora_B[,-1])
space_C <- cbind(Nest_C, Fora_C[,-1])

# Spatial similarity
space_A <- space_A[space_A$X %in% node_A$antID,]
space_B <- space_B[space_B$X %in% node_B$antID,]
space_C <- space_C[space_C$X %in% node_C$antID,]

space_PCA_A <- prcomp(space_A[,-1])
space_PCA_B <- prcomp(space_B[,-1])
space_PCA_C <- prcomp(space_C[,-1])
space_PCA_O <- prcomp(space_O[,-1])

jpeg('Space_PCA_A.jpg', width=2000, height=2000, unit='px')
par(mfrow = c(1,1), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
autoplot(space_PCA_A) + 
  geom_point(size = 10) + 
  theme(text = element_text(size=80),
        axis.text = element_text(size=0))
dev.off()
jpeg('Space_PCA_B.jpg', width=2000, height=2000, unit='px')
par(mfrow = c(1,1), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
autoplot(space_PCA_B) + 
  geom_point(size = 10) + 
  theme(text = element_text(size=80),
        axis.text = element_text(size=0))
dev.off()
jpeg('Space_PCA_C.jpg', width=2000, height=2000, unit='px')
par(mfrow = c(1,1), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
autoplot(space_PCA_C) + 
  geom_point(size = 10) + 
  theme(text = element_text(size=80),
        axis.text = element_text(size=0))
dev.off()
jpeg('Space_PCA_O.jpg', width=2000, height=2000, unit='px')
par(mfrow = c(1,1), mar = c(15,15,10,10), oma = c(0,0,0,0), bty="n", mgp = c(8,8,0), family = "serif")
autoplot(space_PCA_O) + 
  geom_point(size = 10) + 
  theme(text = element_text(size=80),
        axis.text = element_text(size=0))
dev.off()

summary(lm(space_PCA_A$x[,1] ~ node_A$maturity))
summary(lm(space_PCA_B$x[,1] ~ node_B$maturity))
summary(lm(space_PCA_C$x[,1] ~ node_C$maturity))
summary(lm(space_PCA_O$x[,1] ~ node_O$soft))

summary(lm(space_PCA_A$x[,1] ~ node_A$Age))
summary(lm(space_PCA_B$x[,1] ~ node_B$Age))
summary(lm(space_PCA_C$x[,1] ~ node_C$Age))
summary(lm(space_PCA_O$x[,1] ~ node_O$age))

summary(lm(space_PCA_A$x[,1] ~ task_PCA_A$x[,1]))
summary(lm(space_PCA_B$x[,1] ~ task_PCA_B$x[,1]))
summary(lm(space_PCA_C$x[,1] ~ task_PCA_C$x[,1]))
summary(lm(space_PCA_O$x[,1] ~ task_PCA_O$x[,1]))

task_A$taskPCA <- task_PCA_A$x[,1]
task_B$taskPCA <- task_PCA_B$x[,1]
task_C$taskPCA <- task_PCA_C$x[,1]
task_O$taskPCA <- task_PCA_O$x[,1]

task_ALL$taskPCA_ALL <- task_PCA_ALL$x[,1]

task_A$spacePCA <- space_PCA_A$x[,1]
task_B$spacePCA <- space_PCA_B$x[,1]
task_C$spacePCA <- space_PCA_C$x[,1]
task_O$spacePCA <- space_PCA_O$x[,1]

task_A$taskPCA_ALL <- task_ALL[str_extract(rownames(task_ALL), "[A-Z]+" ) == "A",]$taskPCA_ALL
task_B$taskPCA_ALL <- task_ALL[str_extract(rownames(task_ALL), "[A-Z]+" ) == "B",]$taskPCA_ALL
task_C$taskPCA_ALL <- task_ALL[str_extract(rownames(task_ALL), "[A-Z]+" ) == "C",]$taskPCA_ALL
task_O$taskPCA_ALL <- task_ALL[str_extract(rownames(task_ALL), "[A-Z]+" ) == "O",]$taskPCA_ALL

task_A$Age <- node_A$Age
task_B$Age <- node_B$Age
task_C$Age <- node_C$Age
task_O$Age <- node_O$age

task_A$antID <- node_A$antID
task_B$antID <- node_B$antID
task_C$antID <- node_C$antID
task_O$antID <- node_O$antID

task_A$Tag <- node_A$Tag
task_B$Tag <- node_B$Tag
task_C$Tag <- node_C$Tag
task_O$Tag <- node_O$antID

task_A$Extraction <- node_A$Extraction
task_B$Extraction <- node_B$Extraction
task_C$Extraction <- node_C$Extraction

task_A$maturity <- node_A$maturity
task_B$maturity <- node_B$maturity
task_C$maturity <- node_C$maturity
task_O$maturity <- node_O$soft
  
  
# Microbiota
####################

# Separate microbiota data by colony
Mic_A <- MiBi_ABC[str_extract(rownames(MiBi_ABC), "[A-Z]+" ) == "A",]
Mic_B <- MiBi_ABC[str_extract(rownames(MiBi_ABC), "[A-Z]+" ) == "B",]
Mic_C <- MiBi_ABC[str_extract(rownames(MiBi_ABC), "[A-Z]+" ) == "C",]

rownames(Mic_A) <- str_extract(rownames(Mic_A), "[0-9]+")
rownames(Mic_B) <- str_extract(rownames(Mic_B), "[0-9]+")
rownames(Mic_C) <- str_extract(rownames(Mic_C), "[0-9]+")

# Convert sample names back to Tag IDs
rownames(Mic_A) <- as.numeric(rownames(Mic_A)) - 1
rownames(Mic_B) <- as.numeric(rownames(Mic_B)) - 1
rownames(Mic_C) <- as.numeric(rownames(Mic_C)) - 1

# Convert colnames from tagIDs to antIDs
rownamesMic_A <- c()
rownamesMic_B <- c()
rownamesMic_C <- c()
for(i in 1:nrow(Mic_A)){
  tag <- rownames(Mic_A)[i]
  rownamesMic_A <- c(rownamesMic_A, node_A$antID[which(node_A$Tag == tag)])
}
for(i in 1:nrow(Mic_B)){
  tag <- rownames(Mic_B)[i]
  rownamesMic_B <- c(rownamesMic_B, node_B$antID[which(node_B$Tag == tag)])
}
for(i in 1:nrow(Mic_C)){
  tag <- rownames(Mic_C)[i]
  rownamesMic_C <- c(rownamesMic_C, node_C$antID[which(node_C$Tag == tag)])
}
rownames(Mic_A) <- rownamesMic_A
rownames(Mic_B) <- rownamesMic_B
rownames(Mic_C) <- rownamesMic_C

# Trim to microbiota in metadata
Mic_O <- MiBi_O[rownames(MiBi_O) %in% node_O$antID,]

# Reorder rows
Mic_A <- Mic_A[order(as.numeric(rownames(Mic_A))),]
Mic_B <- Mic_B[order(as.numeric(rownames(Mic_B))),]
Mic_C <- Mic_C[order(as.numeric(rownames(Mic_C))),]
Mic_O <- Mic_O[order(as.numeric(rownames(Mic_O))),]

# Remove Blochmannia
Mic_Am <- Mic_A[,-1]
Mic_Bm <- Mic_B[,-1]
Mic_Cm <- Mic_C[,-1]
Mic_Om <- Mic_O[,-1]

# Normalise within rows
for(row in 1:nrow(Mic_A)){
  Mic_A[row,] <-  Mic_A[row,]/sum(Mic_A[row,])
}
for(row in 1:nrow(Mic_B)){
  Mic_B[row,] <-  Mic_B[row,]/sum(Mic_B[row,])
}
for(row in 1:nrow(Mic_C)){
  Mic_C[row,] <-  Mic_C[row,]/sum(Mic_C[row,])
}
for(row in 1:nrow(Mic_O)){
  Mic_O[row,] <-  Mic_O[row,]/sum(Mic_O[row,])
}
for(row in 1:nrow(Mic_Am)){
  if(sum(Mic_Am[row,]) > 0){
    Mic_Am[row,] <-  Mic_Am[row,]/sum(Mic_Am[row,])}
}
for(row in 1:nrow(Mic_Bm)){
  if(sum(Mic_Bm[row,]) > 0){
    Mic_Bm[row,] <-  Mic_Bm[row,]/sum(Mic_Bm[row,])}
}
for(row in 1:nrow(Mic_Cm)){
  if(sum(Mic_Cm[row,]) > 0){
    Mic_Cm[row,] <-  Mic_Cm[row,]/sum(Mic_Cm[row,])}
}
for(row in 1:nrow(Mic_Om)){
  if(sum(Mic_Om[row,]) > 0){
    Mic_Om[row,] <-  Mic_Om[row,]/sum(Mic_Om[row,])}
}

Mic_PCA_A <- prcomp(Mic_A)
Mic_PCA_B <- prcomp(Mic_B)
Mic_PCA_C <- prcomp(Mic_C)
Mic_PCA_O <- prcomp(Mic_O)
Mic_PCA_Am <- prcomp(Mic_Am)
Mic_PCA_Bm <- prcomp(Mic_Bm)
Mic_PCA_Cm <- prcomp(Mic_Cm)
Mic_PCA_Om <- prcomp(Mic_Om)

Mic_PC1_A <- sqrt((Mic_PCA_A$x[,1] - min(Mic_PCA_A$x[,1])) + 0.001)
Mic_PC1_B <- sqrt((Mic_PCA_B$x[,1] - min(Mic_PCA_B$x[,1])) + 0.001)
Mic_PC1_C <- sqrt((Mic_PCA_C$x[,1] - min(Mic_PCA_C$x[,1])) + 0.001)
Mic_PC1_O <- sqrt((Mic_PCA_O$x[,1] - min(Mic_PCA_O$x[,1])) + 0.001)

Mic_PC1_Am <- sqrt((Mic_PCA_Am$x[,1] - min(Mic_PCA_Am$x[,1])) + 0.001)
Mic_PC1_Bm <- sqrt((Mic_PCA_Bm$x[,1] - min(Mic_PCA_Bm$x[,1])) + 0.001)
Mic_PC1_Cm <- sqrt((Mic_PCA_Cm$x[,1] - min(Mic_PCA_Cm$x[,1])) + 0.001)
Mic_PC1_Om <- sqrt((Mic_PCA_Om$x[,1] - min(Mic_PCA_Om$x[,1])) + 0.001)

node_At <- node_A[node_A$antID %in% names(Mic_PC1_Am),]
node_Bt <- node_B[node_B$antID %in% names(Mic_PC1_Bm),]
node_Ct <- node_C[node_C$antID %in% names(Mic_PC1_Cm),]
node_Ot <- node_O[node_O$antID %in% names(Mic_PC1_Om),]

node_At$MicPCA <- Mic_PC1_A
node_Bt$MicPCA <- Mic_PC1_B
node_Ct$MicPCA <- Mic_PC1_C
node_Ot$MicPCA <- Mic_PC1_O

node_At$MicPCAm <- Mic_PC1_Am
node_Bt$MicPCAm <- Mic_PC1_Bm
node_Ct$MicPCAm <- Mic_PC1_Cm
node_Ot$MicPCAm <- Mic_PC1_Om

node_Atl <- node_At[,colnames(node_At) %in% c("antID", "MicPCA", "MicPCAm")]
node_Btl <- node_Bt[,colnames(node_Bt) %in% c("antID", "MicPCA", "MicPCAm")]
node_Ctl <- node_Ct[,colnames(node_Ct) %in% c("antID", "MicPCA", "MicPCAm")]
node_Otl <- node_Ot[,colnames(node_Ot) %in% c("antID", "MicPCA", "MicPCAm")]

task_A_full <- merge(task_A, node_Atl, by= "antID", all = TRUE)
task_B_full <- merge(task_B, node_Btl, by= "antID", all = TRUE)
task_C_full <- merge(task_C, node_Ctl, by= "antID", all = TRUE)
task_O_full <- merge(task_O, node_Otl, by= "antID", all = TRUE)


setwd(MAINDIR)
write.csv(task_A_full, "md_A3.csv", row.names=FALSE)
write.csv(task_B_full, "md_B3.csv", row.names=FALSE)
write.csv(task_C_full, "md_C3.csv", row.names=FALSE)
write.csv(task_O_full, "md_O3.csv", row.names=FALSE)
