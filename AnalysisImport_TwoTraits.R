ped <- read.table("/home/jana/PedigreeAndGeneticValues_cat.txt", header=TRUE, sep=" ")


library(tidyverse)
ped <- read_table("/home/jana/PedigreeAndGeneticValues.txt")
# head(as.data.frame(ped[ped$Generation == 21,]))
# head(as.data.frame(ped[ped$Generation == 20,]))
# head(as.data.frame(ped[ped$Generation == 19,]))
# head(as.data.frame(ped[ped$Generation == 10,]))
# head(as.data.frame(ped[ped$Generation == 1,]))
# head(as.data.frame(ped[ped$Generation == 2,]))
# head(as.data.frame(ped[ped$Generation == 3,]))
# head(as.data.frame(ped[ped$Generation == 4,]))
# head(as.data.frame(ped[ped$Generation == 5,]))
# head(as.data.frame(ped[ped$Generation == 10,]))
# head(as.data.frame(ped[ped$Generation == 9,]))
# head(as.data.frame(ped[ped$Generation == 8,]))
# head(as.data.frame(ped[ped$Generation == 7,]))
# head(as.data.frame(ped[ped$Generation == 6,]))    
ped <- read_table("SimulatedData/PedigreeAndGeneticValues.txt")
popsplit <- read.csv("PopulationSplit.txt", header=TRUE)
colnames(popsplit) <- c("Group", "Indiv")

ped <- merge(ped, popsplit, by="Indiv")

library(ggplot2)

library(reshape)
PED <- ped[,c("Indiv", "Generation", "Group", "gvNormUnres1", "gvNormUnres2", "gvNormUnres3")]
pedm <- melt(PED, id.vars = c("Indiv", "Generation", "Group"))

library(dplyr)
pedma <- pedm %>% group_by(Generation, Group, variable) %>% summarize(Gain = mean(value))

pdf("PlotGain.pdf")
ggplot(data = pedma, aes(x = Generation, y = Gain, group = Group, colour = Group)) + geom_line() + facet_grid(. ~ variable)
dev.off()



###
library(ggplot2)
library(plyr)
library(dplyr)


setwd("/home/jana/Documents/PhD/Projects/inProgress/GenomicStrategies_Import/")
#Read in the TGVsAll tables with genetic trends and variances by generation 
#This is for BurnIn
TGVsAll_burnIn <- read.table("TGVsAll_import_BurnIn_28092020.csv", header=T)
head(TGVsAll_burnIn)
# Plot the genetic gain by Trait (this is actually correlation, since trait 2 has 0.9 and trait 3 0.8 correlation) by group
ggplot(data = TGVsAll_burnIn, aes(x=Generation, y=zMean, colour=Group, group=Group)) + geom_line() + 
  facet_grid(cols = vars(Trait), rows=vars(variable))
# extract just the one trait for one population 
# home population performs selection on trait 1 and import on either trait 2 (r=0.9) or trait 3 (r = 0.8)
TGVsAll_burnIn$scearioSpec <- paste(TGVsAll_burnIn$Group, TGVsAll_burnIn$variable, TGVsAll_burnIn$Trait, sep="_")
burnIn_oneTrait <- TGVsAll_burnIn[TGVsAll_burnIn$scearioSpec %in% c("home_gvNormUnres1_2", "home_gvNormUnres1_3", "import_gvNormUnres2_2", "import_gvNormUnres3_3"), ]
head(burnIn_oneTrait)
burnIn_oneTrait$Trait <- as.character(burnIn_oneTrait$Trait)
burnIn_oneTrait$Trait <- plyr::revalue(burnIn_oneTrait$Trait, c("2" = "r = 0.9", "3" = "r = 0.8"))
ggplot(data = burnIn_oneTrait, aes(x=Generation, y=zMean, colour=Group, group=Group)) + geom_line() + 
  facet_grid(cols = vars(Trait)) +
  theme_bw(base_size = 18) + ylab("Genetic gain")
burnIn_oneTrait[burnIn_oneTrait$Generation == 40,]



#This is for scenario
TGVsAll <- read.table("TGVsAll_import_06102020_noDelay.csv", header=T)
table(TGVsAll$Import)
head(TGVsAll)

#Plot the reps
TGVsAll$Rep <- as.factor(TGVsAll$Rep)
ggplot(data=TGVsAll[TGVsAll$Group == "home",], aes(x=Generation, y=zMean, group=Rep, colour=Rep)) + geom_line() + facet_grid(cols = vars(variable)) +
  theme_bw(base_size = 18) + ylab("Genetic gain") + facet_grid(rows = vars(variable), cols = vars(Import)) 
ggplot(data=TGVsAll[TGVsAll$Group == "import",], aes(x=Generation, y=zMean, group=Rep, colour=Rep)) + geom_line() + facet_grid(cols = vars(variable)) +
  theme_bw(base_size = 18) + ylab("Genetic gain") + facet_grid(rows = vars(variable), cols = vars(Import)) 

#Create aggregated value for mean true genetic value (value), standardised genetic value (zMean), standardised genetic variance (SDSt) and standardised genic variance (SDGenicSt)
TGVAvg <- TGVsAll %>% group_by(Scenario, Import, Generation, Group, variable) %>% 
          dplyr::summarise(meanTGV = mean(value), 
                           sdTGV = sd(value), 
                           meanZTGV = mean(zMean), 
                           sdZTGV = sd(zMean), 
                           meanGeneticSD = mean(sd),
                           meanzGeneticSD = mean(SDSt), 
                           meanGenicSD = mean(SDGenic), 
                           meanzGenicSD = mean(SDGenicSt))
#TGVAvg[TGVAvg$Generation == 40,]

# Revalue the variables so that TGV for trait 1 says Trait 1 and the one for trait 2 says Trait 2
TGVAvg$variable <- plyr::revalue(TGVAvg$variable, c("gvNormUnres1" = "Trait 1", "gvNormUnres3" = "Trait 2"))
# Correct order of the traits
TGVAvg$variable <- factor(TGVAvg$variable, c( "Trait 1", "Trait 2"))
# Correct order of the import
TGVAvg$Import <- factor(TGVAvg$Import, c("0_0", "10_10", "50_50", "100_100", "0_100"))
# Plot genetic gain
ggplot(data=TGVAvg, aes(x=Generation, y=meanZTGV, group=Group, colour=Group)) + geom_line() + facet_grid(cols = vars(variable)) +
  theme_bw(base_size = 18) + ylab("Genetic gain") + facet_grid(rows = vars(variable), cols = vars(Import)) + 
  geom_ribbon(aes(ymin = meanZTGV - sdZTGV, ymax = meanZTGV + sdZTGV, fill = Group),  linetype = 0, alpha = 0.3)
# Now plot the genetic gain for trait 1 in home and trait 2 in import - since this is condidered the same trait
head(TGVAvg)
TGVAvg_oneTrait <- rbind(TGVAvg[TGVAvg$Group == "import" & TGVAvg$variable == "Trait 2",], TGVAvg[TGVAvg$Group == "home" & TGVAvg$variable == "Trait 1",])
TGVAvg_oneTrait$variable <- "Correlated trait"
# Plot genetic gain for just the one correlated trait
ggplot(data=TGVAvg_oneTrait, aes(x=Generation, y=meanZTGV, group=Group, colour=Group)) + geom_line() + 
  theme_bw(base_size = 18) + ylab("Genetic gain") + facet_grid(cols = vars(Import)) + 
  scale_y_continuous(breaks = c(seq(3, 9, 0.5))) +
  geom_ribbon(aes(ymin = meanZTGV - sdZTGV, ymax = meanZTGV + sdZTGV, fill = Group),  linetype = 0, alpha = 0.3)
# Plot genetic gain for just the one correlated trait for generations 40-60
ggplot(data=TGVAvg_oneTrait[TGVAvg_oneTrait$Generation > 29,], aes(x=Generation, y=meanZTGV, group=Group, colour=Group)) + geom_line() + 
  theme_bw(base_size = 18) + ylab("Genetic gain") + facet_grid(cols = vars(Import)) + 
  scale_y_continuous(breaks = c(seq(3, 9, 0.5))) +
  #geom_vline(xintercept = 50) +
  geom_hline(yintercept = TGVAvg_oneTrait$meanZTGV[TGVAvg_oneTrait$Import == "100_100" & TGVAvg_oneTrait$Generation == 50 & TGVAvg_oneTrait$Group == "home"]) + 
  geom_ribbon(aes(ymin = meanZTGV - sdZTGV, ymax = meanZTGV + sdZTGV, fill = Group),  linetype = 0, alpha = 0.3)

compareDataH <- TGVAvg_oneTrait[TGVAvg_oneTrait$Generation > 29 & TGVAvg_oneTrait$Group =="home" & TGVAvg_oneTrait$Import == "0_0",]
compareDataI <- TGVAvg_oneTrait[TGVAvg_oneTrait$Generation > 29 & TGVAvg_oneTrait$Group =="import"& TGVAvg_oneTrait$Import == "0_0",]
modelH <- lm(data =compareDataH, formula = meanZTGV ~ Generation)
modelI <- lm(data =compareDataI, formula = meanZTGV ~ Generation)
modelH$coefficients
modelI$coefficients

#Plot genetic SD
ggplot(data=TGVAvg, aes(x=Generation, y=meanzGeneticSD, group=Group, colour=Group)) + geom_line() + facet_grid(cols = vars(variable)) +
  theme_bw(base_size = 18) + ylab("Genetic sd") + facet_grid(rows = vars(variable), cols = vars(Import))
#Plot genetic SD - just for the one correlated trait
ggplot(data=TGVAvg_oneTrait, aes(x=Generation, y=meanzGeneticSD, group=Group, colour=Group)) + geom_line() + facet_grid(cols = vars(variable)) +
  theme_bw(base_size = 18) + ylab("Genetic sd") + facet_grid(rows = vars(variable), cols = vars(Import))
#Plot genic sd
ggplot(data=TGVAvg[TGVAvg$Generation > 40,], aes(x=Generation, y=meanGenicSD, group=Group, colour=Group)) + geom_line() + facet_grid(cols = vars(variable)) +
  theme_bw(base_size = 18) + ylab("Genic sd") + facet_grid(rows = vars(variable), cols = vars(Import))


###############################################
###############################################
# Check the accuracies
###############################################
###############################################
acc <- read.csv("Accuracies.csv")
ggplot(data = acc, aes(x=Gen, y=Cor, group=Cat, colour=Cat)) + 
        geom_line() + 
        facet_grid(cols=vars(Scenario), rows = vars(Group))

#Only last generation accuracies
acc <- read.csv("Acc_LastGen.csv")
acc <- read.csv("/home/jana/EddieDir/10K/SU55_import/AccuraciesPost.csv")
nrow(acc)
head(acc)
# accA <- as.data.frame(acc %>% group_by(CatAge, Group, Scenario) %>% dplyr::summarise(accH = mean(accEBVh), accI = mean(accEBVi)))
# acc_m <- melt(accA, id.vars = c("CatAge", "Group", "Scenario"))
# ggplot(data = acc_m, aes(x=CatAge, y=value, group=variable, colour=variable)) + geom_point() + facet_grid(rows = vars(Group))
ggplot(data = acc[!acc$Category %in% c("Male candidates1", "izl"),], aes(x=Generation, y=Cor, group=Category, colour=Category)) + 
  geom_line() + 
  facet_grid(rows = vars(Group), cols=vars(Scenario)) + 
  theme_bw(base_size=16)

acc <- read.csv("/home/jana/EddieDir/10K/SU55_import/AccuraciesPost_rep.csv")
accA <- acc %>%  group_by(Category, Generation, Group, Scenario) %>% dplyr::summarise(meanAcc = mean(Cor), sdAcc=sd(Cor))
head(accA)
ggplot(data = accA[!accA$Category %in% c("Male candidates1", "izl"),], aes(x=Generation, y=meanAcc, group=Category, colour=Category)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = meanAcc - sdAcc, ymax = meanAcc + sdAcc, fill = Category),  linetype = 0, alpha = 0.3)+
  facet_grid(rows = vars(Group), cols=vars(Scenario)) + 
  theme_bw(base_size=16)


ggplot(data = accA[!accA$Category %in% c("Male candidates1", "izl"),], aes(x=Generation, y=meanAcc, group=Group, colour=Group)) + 
  geom_line() + 
  #geom_ribbon(aes(ymin = meanAcc - sdAcc, ymax = meanAcc + sdAcc, fill = Group),  linetype = 0, alpha = 0.3)+
  facet_grid(rows = vars(Category), cols=vars(Scenario)) + 
  theme_bw(base_size=16)
###############################################
###############################################
# Check the EBVs (last estimation)
###############################################
###############################################
ebvs <- read.csv("EBVs.csv")
ebvs_m <- melt(ebvs, id.vars = c("Generation", "Group", "scenario"))
ggplot(data = ebvs_m, aes(x=Generation, y=value, group=variable, colour=variable)) + 
        geom_line() + 
        facet_grid(cols=vars(scenario), rows = vars(Group))


#Check the TGVs of home-imported bulls
#Check the TGVs for all genomically tested bulls
gpb <- read.csv("AllBulls_TGVs.csv")
gpb_m <- melt(gpb, id.vars = c("scenario", "Group", "Generation"))
ggplot(data = gpb_m, aes(x=Generation, y = value, group=variable, fill=variable)) + geom_bar(stat="identity", position="dodge") + facet_grid(rows = vars(Group), cols = vars(scenario))
# Combine home-trait 1 and import-trait 3
gpb_m_oneT <- rbind(gpb_m[gpb_m$Group == "home" & gpb_m$variable == "gvNormUnres1",], gpb_m[gpb_m$Group == "import" & gpb_m$variable == "gvNormUnres3",])
# Plot just the one correlated trait
ggplot(data = gpb_m_oneT, aes(x=Generation, y = value, group=variable, fill=variable)) + geom_bar(stat="identity", position="dodge") + facet_grid(cols = vars(scenario))


#Check the TGVs of bulls used in the home population
home_gpb <- read.csv("HomeFathers_TGVs.csv")
home_gpb_m <- melt(home_gpb, id.vars = c("scenario", "Group", "Generation"))
home_gpb_m$variable <- plyr::revalue(home_gpb_m$variable, c("gvNormUnres1" = "Trait 1","gvNormUnres3" = "Trait 2"))
ggplot(data = home_gpb_m, aes(x=Generation, y = value, group=Group, fill=Group)) + 
        geom_bar(stat="identity", position="dodge") + facet_grid(rows = vars(variable), cols = vars(scenario)) +
        ylab("True genetic value") + 
        scale_fill_manual("", labels = c("Domestic", "Import"), values = c("#199fd4", "#e05a80")) + 
        theme_bw(base_size = 16)
# Combine home-trait 1 and import-trait 3
home_gpb_m_oneT <- rbind(home_gpb_m[home_gpb_m$Group == "home" & home_gpb_m$variable == "gvNormUnres1",], 
                    home_gpb_m[home_gpb_m$Group == "import" & home_gpb_m$variable == "gvNormUnres3",])
table(home_gpb_m_oneT$Group, home_gpb_m_oneT$variable)
# Plot just the one correlated trait, you have only one trait per Group, so trait = Group
ggplot(data = home_gpb_m_oneT, aes(x=Generation, y = value, group=Group, fill=Group)) + 
        geom_bar(stat="identity", position="dodge") + facet_grid(cols = vars(scenario)) + 
        ylab("True genetic value") + 
        scale_fill_manual("", labels = c("Domestic", "Import"), values = c("#199fd4", "#e05a80")) + 
        theme_bw(base_size = 16)

#Read in the dataframe with the EBVs  
# For home population
gh <- read.csv("GenPed_EBVhome.txt")
gh$Group <- "Home"
# For import population
gi <- read.csv("GenPed_EBVimport.txt")
gi$Group <- "Import"
# Bind them together
g <- rbind(gh, gi)
#Aggregate EBVs by Generation and Group 
gs <- g %>% group_by(Generation, Group) %>% summarise(meanEBV = mean(EBV))
#Plot the EBVs
ggplot(data = gs, aes(x = Generation, y = meanEBV, group = Group, colour = Group)) + geom_line()



#Scenario
ped <- read.table("PedigreeAndGeneticValues_cat.txt", header=T, sep=" ")
ps <- read.csv("PopulationSplit.txt")
colnames(ps) <- c("Group", "Indiv")

PED <- merge(ped, ps, by="Indiv")
#check fathers --> where do they come from
table(PED$Group[PED$Indiv %in% PED$Father[PED$cat == "potomciNP" & PED$Group == "home"]])
table(PED$Generation[PED$Indiv %in% PED$Father[PED$cat == "potomciNP" & PED$Group == "home"]])
table(PED$Group[PED$Indiv %in% PED$Father[PED$cat == "nr" & PED$Group == "home"]])
table(PED$Group[PED$Indiv %in% PED$Father[PED$cat == "potomciNP" & PED$Group == "import"]])
table(PED$Group[PED$Indiv %in% PED$Father[PED$cat == "nr" & PED$Group == "import"]])
table(PED$Group[PED$Indiv %in% PED$Father[PED$Group == "home"]])
######

PEDs <- PED[PED$Generation > 39,]

library(reshape)
PEDs <- PEDs[,c("Indiv", "Generation", "Group", "gvNormUnres1", "gvNormUnres2", "gvNormUnres3")]
pedm <- melt(PEDs, id.vars = c("Indiv", "Generation", "Group"))

pedma <- pedm %>% group_by(Generation, Group, variable) %>% summarize(Gain = mean(value), GainSD = sd(value))

pedma <- pedma[pedma$variable %in% c("gvNormUnres1", "gvNormUnres3"),]
pedma$variable <- plyr::revalue(pedma$variable, c("gvNormUnres1" = "Trait 1", "gvNormUnres3" = "Trait 2"))
ggplot(data = pedma, aes(x = Generation, y = Gain, group = Group, colour = Group)) + geom_line() + facet_grid(. ~ variable) +
  theme_bw(base_size = 18) + ylab("Genetic gain")
ggplot(data = pedma, aes(x = Generation, y = GainSD, group = Group, colour = Group)) + geom_line() + facet_grid(. ~ variable) +
  theme_bw(base_size = 18) + ylab("Genetic gain")


###check the EBVs from BurnIn
meanEBV <- read.csv("MeanEBVs.txt")
head(meanEBV)
ggplot(data = meanEBV, aes(x=Generation, y = EBVi, colour=Group, group=Group)) + geom_line() +facet_grid(. ~YearEBV)
ggplot(data = meanEBV, aes(x=Generation, y = EBVh, colour=Group, group=Group)) + geom_line() +facet_grid(. ~YearEBV)


#Check EBVs in the BurnIn
bpedH <- read.csv("BurnIn/GenPed_EBVhome.txt")
bpedI <- read.csv("BurnIn/GenPed_EBVimport.txt")
bpedH$Group <- "home"
bpedI$Group <- "import"
head(bpedI)
head(bpedH)
tail(bpedI)
tail(bpedH)

bped <- rbind(bpedH, bpedI)
library(dplyr)
bpedM <- bped %>% group_by(Generation, Group) %>% summarise(meanTGV1=mean(gvNormUnres1), meanTGV3 = mean(gvNormUnres3), meanEBV = mean(EBV))
ggplot(data=bpedM, aes(x=Generation, y=meanTGV1, group=Group, colour=Group)) + geom_line()
ggplot(data=bpedM, aes(x=Generation, y=meanTGV3, group=Group, colour=Group)) + geom_line()
ggplot(data=bpedM, aes(x=Generation, y=meanEBV, group=Group, colour=Group)) + geom_line()


#Add EBVs from generation 40
bped <- bped[order(bped$Indiv),]
ebvH <- read.table("BurnIn/renumbered_Solutions_home_40")[,-1]
ebvI <- read.table("BurnIn/renumbered_Solutions_import_40")[,-1]
colnames(ebvH) <- c("Indiv", "EBV_T1")
colnames(ebvI) <- c("Indiv", "EBV_T3")
head(ebvH)
tail(ebvH)

bped <- merge(bped, ebvH, by="Indiv")
bped <- merge(bped, ebvI, by="Indiv")
bpedM <- bped %>% group_by(Generation, Group) %>% summarise(meanTGV1=mean(gvNormUnres1), meanTGV3 = mean(gvNormUnres3), meanEBV_T1 = mean(EBV_T1), meanEBV_T2 = mean(EBV_T3))
head(bpedM)
bpedMM <- melt(as.data.frame(bpedM), id.vars = c("Group", "Generation"))
ggplot(data=bpedMM, aes(x=Generation, y=value, group=Group, colour=Group)) + geom_line() + facet_grid(rows = vars(variable))


#Pedigree categories
pedC <- read.table("BurnIn/PedigreeAndGeneticValues_cat.txt", sep=" ", header=TRUE)
pedC <- merge(pedC, bped[,c("Indiv", "Group", "EBV_T1", "EBV_T3")], by="Indiv")
pedC_bulls <- pedC[pedC$cat %in% c("gpb", "pb"),]
table(pedC_bulls$Group)
table(pedC_bulls$Group, pedC_bulls$Generation)
table(pedC_bulls$cat)

pedBulls <- pedC_bulls %>% group_by(Generation, Group) %>% summarise(meanTGV1=mean(gvNormUnres1), meanTGV3 = mean(gvNormUnres3), meanEBV_T1 = mean(EBV_T1), meanEBV_T2 = mean(EBV_T3))
pedBulls[pedBulls$Generation == 34,]
pedBullsM <- melt(as.data.frame(pedBulls), id.vars = c("Group", "Generation"))
ggplot(data=pedBullsM, aes(x=Generation, y=value, group=Group, colour=Group)) + geom_line() + facet_grid(rows = vars(variable))


####################################3
# Check generation intervals
#####################################
gi <- read.csv("GENINTS_all_25092020.csv")
gi$lineSex <- paste0(gi$line, gi$sex)
head(gi)
table(gi$Gen)
gi[gi$Gen==55 & gi$scenario=="10_10",]
ggplot(data=gi, aes(x=Gen, y=genInt, group=lineSex, colour=lineSex)) + geom_line() +
  facet_grid(cols=vars(Group), rows=vars(scenario))

             