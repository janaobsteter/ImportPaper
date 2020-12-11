setwd("/home/jana/Documents/PhD/Projects/inProgress/GenomicStrategies_Import/")
popSplit <- read.csv("PopulationSplit.txt")
colnames(popSplit) <- c("group", "Indiv")
ped <- read.table("PedigreeAndGeneticValues_cat.txt", header=TRUE)

pedHome <- read.csv("GenPed_EBVhome.txt")
pedImport <- read.csv("GenPed_EBVimport.txt")
genped <- read.csv("GenPed_EBV.txt")

PED <- merge(ped, popSplit, by="Indiv")  
PED <- merge(PED, genped[,c("Indiv", "EBV")], by="Indiv")  
  
PED <- PED[PED$Generation > 39,]
PED <- merge(PED, genped[,c("Indiv", "EBV")], by="Indiv")
table(PED$group)

library(AlphaPart)
PED$partGroup <- paste0(PED$group, PED$sex)

Part = AlphaPart(x = PED, sort = FALSE,
                  colId = "Indiv", colFid = "Father", colMid = "Mother",
                  colPath = "partGroup", colAGV = "gvNormUnres1")


PartSummary = summary(object = Part, by = "Generation")
plot(PartSummary)

PartHome <- Part
PartHome$gvNormUnres1$Generation <- as.factor(PartHome$gvNormUnres1$Generation)
PartHome$gvNormRestr1 <- PartHome$gvNormUnres1[PartHome$gvNormUnres1$group == "home",]
nrow(PartHome$gvNormRestr1)
table(PartHome$gvNormRestr1$Generation)

PartHomeSummary <- summary(object = PartHome, by = "Generation")
head(PartHomeSummary$gvNormUnres1$abs)
plot(PartHomeSummary)

#summary
homeF_h <- aggregate(PartHome$gvNormRestr1$gvNormUnres1_homeF ~ PartHome$gvNormRestr1$Generation, FUN="mean")
colnames(homeF_h) <- c("Generation", "Value")
homeF_h$group <- "homeF"

homeM_h <- aggregate(PartHome$gvNormRestr1$gvNormUnres1_homeM ~ PartHome$gvNormRestr1$Generation, FUN="mean")
colnames(homeM_h) <- c("Generation", "Value")
homeM_h$group <- "homeM"

importF_h <- aggregate(PartHome$gvNormRestr1$gvNormUnres1_importF ~ PartHome$gvNormRestr1$Generation, FUN="mean")
colnames(importF_h) <- c("Generation", "Value")
importF_h$group <- "importF"

importM_h <- aggregate(PartHome$gvNormRestr1$gvNormUnres1_importM ~ PartHome$gvNormRestr1$Generation, FUN="mean")
colnames(importM_h) <- c("Generation", "Value")
importM_h$group <- "importM"

homePart <- do.call(rbind, list(homeF_h, homeM_h, importF_h, importM_h))
homePart$Generation <- as.numeric(homePart$Generation)
ggplot(data=homePart, aes(x=Generation, y = Value, colour=group)) + geom_path()


PartImport <- Part
PartImport$gvNormUnres1$Generation <- as.factor(PartImport$gvNormUnres1$Generation)
PartImport$gvNormRestr1 <- PartImport$gvNormUnres1[PartImport$gvNormUnres1$group == "import",]
PartImportSummary <- summary(PartImport, by = "Generation")
plot(PartImportSummary)

homeF_i <- aggregate(PartImport$gvNormRestr1$gvNormUnres1_homeF ~ PartImport$gvNormRestr1$Generation, FUN="mean")
colnames(homeF_i) <- c("Generation", "Value")
homeF_i$group <- "homeF"

homeM_i <- aggregate(PartImport$gvNormRestr1$gvNormUnres1_homeM ~ PartImport$gvNormRestr1$Generation, FUN="mean")
colnames(homeM_i) <- c("Generation", "Value")
homeM_i$group <- "homeM"

importF_i <- aggregate(PartImport$gvNormRestr1$gvNormUnres1_importF ~ PartImport$gvNormRestr1$Generation, FUN="mean")
colnames(importF_i) <- c("Generation", "Value")
importF_i$group <- "importF"

importM_i <- aggregate(PartImport$gvNormRestr1$gvNormUnres1_importM ~ PartImport$gvNormRestr1$Generation, FUN="mean")
colnames(importM_i) <- c("Generation", "Value")
importM_i$group <- "importM"

importPart <- do.call(rbind, list(homeF_i, homeM_i, importF_i, importM_i))
importPart$Generation <- as.numeric(importPart$Generation)
ggplot(data=importPart, aes(x=Generation, y = Value, colour=group)) + geom_path()



##############################################################################
##############################################################################
setwd("/home/jana/Documents/PhD/Projects/inProgress/GenomicStrategies_Import/Import50BM/")
popSplit <- read.csv("PopulationSplit.txt")
colnames(popSplit) <- c("group", "Indiv")
ped <- read.table("PedigreeAndGeneticValues_cat.txt", header=TRUE)

pedHome <- read.csv("GenPed_EBVhome.txt")
pedImport <- read.csv("GenPed_EBVimport.txt")
genped <- read.csv("GenPed_EBV.txt")


PED <- merge(ped, popSplit, by="Indiv")  
PED <- merge(PED, genped[,c("Indiv", "EBV")], by="Indiv")  
  
PED <- PED[PED$Generation > 39,]
PED <- merge(PED, genped[,c("Indiv", "EBV")], by="Indiv")
table(PED$group)

library(AlphaPart)
PED$partGroup <- paste0(PED$group, PED$sex)

Part = AlphaPart(x = PED, sort = FALSE,
                  colId = "Indiv", colFid = "Father", colMid = "Mother",
                  colPath = "partGroup", colAGV = "gvNormUnres1")


PartSummary = summary(object = Part, by = "Generation")
plot(PartSummary)

PartHome <- Part
PartHome$gvNormUnres1$Generation <- as.factor(PartHome$gvNormUnres1$Generation)
PartHome$gvNormRestr1 <- PartHome$gvNormUnres1[PartHome$gvNormUnres1$group == "home",]
nrow(PartHome$gvNormRestr1)
table(PartHome$gvNormRestr1$Generation)

PartHomeSummary <- summary(object = PartHome, by = "Generation")
head(PartHomeSummary$gvNormUnres1$abs)
plot(PartHomeSummary)

#summary
homeF_h <- aggregate(PartHome$gvNormRestr1$gvNormUnres1_homeF ~ PartHome$gvNormRestr1$Generation, FUN="mean")
colnames(homeF_h) <- c("Generation", "Value")
homeF_h$group <- "homeF"

homeM_h <- aggregate(PartHome$gvNormRestr1$gvNormUnres1_homeM ~ PartHome$gvNormRestr1$Generation, FUN="mean")
colnames(homeM_h) <- c("Generation", "Value")
homeM_h$group <- "homeM"

importF_h <- aggregate(PartHome$gvNormRestr1$gvNormUnres1_importF ~ PartHome$gvNormRestr1$Generation, FUN="mean")
colnames(importF_h) <- c("Generation", "Value")
importF_h$group <- "importF"

importM_h <- aggregate(PartHome$gvNormRestr1$gvNormUnres1_importM ~ PartHome$gvNormRestr1$Generation, FUN="mean")
colnames(importM_h) <- c("Generation", "Value")
importM_h$group <- "importM"

homePart <- do.call(rbind, list(homeF_h, homeM_h, importF_h, importM_h))
homePart$Generation <- as.numeric(homePart$Generation)
ggplot(data=homePart, aes(x=Generation, y = Value, colour=group)) + geom_path()


PartImport <- Part
PartImport$gvNormUnres1$Generation <- as.factor(PartImport$gvNormUnres1$Generation)
PartImport$gvNormRestr1 <- PartImport$gvNormUnres1[PartImport$gvNormUnres1$group == "import",]
PartImportSummary <- summary(PartImport, by = "Generation")
plot(PartImportSummary)

homeF_i <- aggregate(PartImport$gvNormRestr1$gvNormUnres1_homeF ~ PartImport$gvNormRestr1$Generation, FUN="mean")
colnames(homeF_i) <- c("Generation", "Value")
homeF_i$group <- "homeF"

homeM_i <- aggregate(PartImport$gvNormRestr1$gvNormUnres1_homeM ~ PartImport$gvNormRestr1$Generation, FUN="mean")
colnames(homeM_i) <- c("Generation", "Value")
homeM_i$group <- "homeM"

importF_i <- aggregate(PartImport$gvNormRestr1$gvNormUnres1_importF ~ PartImport$gvNormRestr1$Generation, FUN="mean")
colnames(importF_i) <- c("Generation", "Value")
importF_i$group <- "importF"

importM_i <- aggregate(PartImport$gvNormRestr1$gvNormUnres1_importM ~ PartImport$gvNormRestr1$Generation, FUN="mean")
colnames(importM_i) <- c("Generation", "Value")
importM_i$group <- "importM"

importPart <- do.call(rbind, list(homeF_i, homeM_i, importF_i, importM_i))
importPart$Generation <- as.numeric(importPart$Generation)
ggplot(data=importPart, aes(x=Generation, y = Value, colour=group)) + geom_path()
