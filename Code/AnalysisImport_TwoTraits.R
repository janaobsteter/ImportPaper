library(lme4)
library(tidyr)
library(ggplot2)
library(reshape)
library(lme4)
library(viridis)
library(plyr)
library(dplyr)


setwd("/home/jana/Documents/PhD/Projects/inProgress/GenomicStrategies_Import/")
#Read in the TGVsAll tables with genetic trends and variances by generation 
#This is for BurnIn
TGVsAll_burnIn <- read.table("TGVsAll_import_BurnIn_29092020.csv", header=T)
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


########################################################################
########################################################################
#This is for scenario
TGVsAll <- read.table("Results/TGVsAll_import_14122020_noDelay.csv", header=T)
# Revalue the variables so that TGV for trait 1 says Trait 1 and the one for trait 2 says Trait 2
TGVsAll$variable <- plyr::revalue(TGVsAll$variable, c("gvNormUnres1" = "Trait home", "gvNormUnres2" = "Trait import", "gvNormUnres3" = "Trait import"))
TGVsAll$Trait <- as.factor(TGVsAll$Trait)
TGVsAll$Trait <- plyr::revalue(TGVsAll$Trait, c("2" = "cor = 0.9", "3" = "cor = 0.8"))
# Correct order of the traits
TGVsAll$variable <- factor(TGVsAll$variable, c( "Trait home", "Trait import"))
# Correct order of the import
TGVsAll$Import <- factor(TGVsAll$Import, c("0_0", "10_10", "25_25", "50_50", "100_100", "0_100"))
TGVsAll$Import <- plyr::revalue(TGVsAll$Import, c("0_0" = "0", "10_10" = "10", "25_25" = "25", "50_50" = "50", "100_100" = "100", "0_100" = "100BD"))
TGVsAll$PlotGroup <- ifelse(TGVsAll$Group == "home", paste0(TGVsAll$Group, "_", TGVsAll$Import), TGVsAll$Group)
table(TGVsAll$PlotGroup)
table(TGVsAll$Import)
table(TGVsAll$Rep)
head(TGVsAll)

#Plot the reps
TGVsAll$Rep <- as.factor(TGVsAll$Rep)
ggplot(data=TGVsAll[TGVsAll$Group == "home" & TGVsAll$Scenario == "GenGen" & TGVsAll$Trait == 3,], aes(x=Generation, y=zMean, group=Rep, colour=Rep)) + 
  geom_line() + facet_grid(cols = vars(variable)) +
  theme_bw(base_size = 18) + ylab("Genetic gain") + 
  facet_grid(rows = vars(variable), cols = vars(Import)) 
ggplot(data=TGVsAll[TGVsAll$Group == "home" & TGVsAll$Scenario == "ClassGen" & TGVsAll$Trait == 3,], aes(x=Generation, y=zMean, group=Rep, colour=Rep)) + 
  geom_line() + facet_grid(cols = vars(variable)) +
  theme_bw(base_size = 18) + ylab("Genetic gain") + 
  facet_grid(rows = vars(variable), cols = vars(Import)) 
ggplot(data=TGVsAll[TGVsAll$Group == "home" & TGVsAll$Scenario == "GenGen" & TGVsAll$Trait == 2,], aes(x=Generation, y=zMean, group=Rep, colour=Rep)) + 
  geom_line() + facet_grid(cols = vars(variable)) +
  theme_bw(base_size = 18) + ylab("Genetic gain") + 
  facet_grid(rows = vars(variable), cols = vars(Import)) 
ggplot(data=TGVsAll[TGVsAll$Group == "home" & TGVsAll$Scenario == "ClassGen" & TGVsAll$Trait == 2,], aes(x=Generation, y=zMean, group=Rep, colour=Rep)) + 
  geom_line() + facet_grid(cols = vars(variable)) +
  theme_bw(base_size = 18) + ylab("Genetic gain") + 
  facet_grid(rows = vars(variable), cols = vars(Import)) 


#Create aggregated value for mean true genetic value (value), standardised genetic value (zMean), standardised genetic variance (SDSt) and standardised genic variance (SDGenicSt)
TGVsAll[TGVsAll$Generation == 30 & TGVsAll$Group == "home" & TGVsAll$Import == "10_10" & TGVsAll$variable == "Trait home",]
head(TGVsAll)
head(TGVsAll[TGVsAll$Generation > 29,])
summary(TGVsAll$zMean30)
TGVAvg <- TGVsAll %>% group_by(Scenario, Generation, PlotGroup, Group, variable, Trait) %>%
                 summarise(meanTGV = mean(value), 
                           sdTGV = sd(value), 
                           meanZTGV = mean(zMean), 
                           sdZTGV = sd(zMean), 
                           zMean20 = mean(zMean_noC), 
                           zSd20 = sd(zMean_noC), 
                           zSd30 = sd(zMean30),
                           zMean30 = mean(zMean30), 
                           meanGeneticSD = mean(sd),
                           meanzGeneticSD = mean(SDSt), 
                           meanGenicSD = mean(SDGenic), 
                           meanzGenicSD = mean(SDGenicSt))
#TGVAvg[TGVAvg$Generation == 40,]

# Revalue the variables so that TGV for trait 1 says Trait 1 and the one for trait 2 says Trait 2
head(TGVAvg)
head(TGVAvg[TGVAvg$Generation > 29,])
TGVAvg$Trait <- as.factor(TGVAvg$Trait)
table(TGVAvg$Trait)
TGVsAll$Trait <- as.character(TGVsAll$Trait)

# Plot genetic gain
# ggplot(data=TGVAvg[TGVAvg$Scenario == "GenGen" & TGVAvg$Generation > 29,], aes(x=Generation, y=zMean30, group=Group, colour=Group)) + geom_line() + facet_grid(cols = vars(variable)) +
#   theme_bw(base_size = 18) + ylab("Genetic gain") + 
#   facet_grid(rows = vars(variable, Trait), cols = vars(Import)) + 
#   geom_ribbon(aes(ymin = zMean30 - zSd30, ymax = zMean30 + zSd30, fill = Group),  linetype = 0, alpha = 0.3)
# Now plot the genetic gain for trait 1 in home and trait 2 in import - since this is condidered the same trait
head(TGVAvg)
table(TGVAvg$Group, TGVAvg$variable)
TGVAvg_oneTrait <- rbind(TGVAvg[TGVAvg$Group == "import" & TGVAvg$variable == "Trait import",], 
                         TGVAvg[TGVAvg$Group == "home" & TGVAvg$variable == "Trait home",])
table(TGVAvg_oneTrait$Group)
TGVAvg_oneTrait$variable <- "Correlated trait"

# Plot genetic gain for just the one correlated trait
table(TGVAvg_oneTrait$PlotGroup)
TGVAvg_oneTrait$PlotGroup <- factor(TGVAvg_oneTrait$PlotGroup, c("home_0", "home_10",   "home_25", "home_50", "home_100", "home_100BD", "import"))
table(TGVAvg_oneTrait$PlotGroup)
TGVAvg_oneTrait$Scenario <- plyr::revalue(TGVAvg_oneTrait$Scenario, c("ClassGen" = "10y_delay", "GenGen" = "No_delay"))
# ggplot(data=TGVAvg_oneTrait[TGVAvg_oneTrait$Scenario == "GenGen" & TGVAvg_oneTrait$Generation > 29,], aes(x=Generation, y=meanZTGV, group=Group, colour=Group)) + geom_line() + 
#   theme_bw(base_size = 18) + ylab("Genetic gain") + facet_grid(cols = vars(Import), rows = vars(Trait)) + 
#   #scale_y_continuous(breaks = c(seq(3, 9, 0.5))) +
#   geom_ribbon(aes(ymin = meanZTGV - sdZTGV, ymax = meanZTGV + sdZTGV, fill = Group),  linetype = 0, alpha = 0.3) + ggtitle("Both populations run genomic selection")
ggplot(data=TGVAvg_oneTrait[TGVAvg_oneTrait$Generation > 29,], aes(x=Generation, y=meanZTGV, group=PlotGroup, colour=PlotGroup)) + 
  geom_line(size=0.9) + 
  theme_bw(base_size = 18) + ylab("Genetic gain") + facet_grid(cols = vars(Trait), rows = vars(Scenario)) + 
  #scale_y_continuous(breaks = c(seq(3, 9, 0.5))) +
  geom_ribbon(aes(ymin = meanZTGV - sdZTGV, ymax = meanZTGV + sdZTGV, fill = PlotGroup),  linetype = 0, alpha = 0.12) + 
  theme(legend.position = "right") + 
  scale_colour_manual("", values = c(viridis::plasma(18)[rev(seq(1, 18, 3))][1:6], 'black')) +
  scale_fill_manual("", values = c(viridis::plasma(18)[rev(seq(1, 18, 3))][1:6], 'black')) 

# Standardized onto generation 30
head(TGVAvg_oneTrait)
ggplot(data=TGVAvg_oneTrait[TGVAvg_oneTrait$Generation > 29,], aes(x=Generation, y=zMean30, group=PlotGroup, colour=PlotGroup)) + 
  geom_line(size=0.9) + 
  theme_bw(base_size = 18) + ylab("Genetic gain") + facet_grid(cols = vars(Trait), rows = vars(Scenario)) + 
  #scale_y_continuous(breaks = c(seq(3, 9, 0.5))) +
  geom_ribbon(aes(ymin = zMean30 - zSd30, ymax = zMean30 + zSd30, fill = PlotGroup),  linetype = 0, alpha = 0.12) + 
  theme(legend.position = "right") + 
  scale_colour_manual("", values = c(viridis::plasma(18)[rev(seq(1, 18, 3))][1:6], 'black')) +
  scale_fill_manual("", values = c(viridis::plasma(18)[rev(seq(1, 18, 3))][1:6], 'black')) 


# Cmpute the slope of the regression lines, you need data by reps
TGV_oneTrait <- rbind(TGVsAll[TGVsAll$Group == "import" & TGVsAll$variable == "Trait import",], 
                         TGVsAll[TGVsAll$Group == "home" & TGVsAll$variable == "Trait home",])
head(TGV_oneTrait)
TGV_oneTrait$regressionGroup <- paste(TGV_oneTrait$Scenario, TGV_oneTrait$PlotGroup, TGV_oneTrait$Rep, sep="-")
head(TGV_oneTrait)
tail(TGV_oneTrait)

gain0.9 <- lmList(zMean30 ~ Generation | regressionGroup, data=TGV_oneTrait[TGV_oneTrait$Trait == "cor = 0.9",])
gain0.8 <- lmList(zMean30 ~ Generation | regressionGroup, data=TGV_oneTrait[TGV_oneTrait$Trait == "cor = 0.8",])
gain0.8
gain <- data.frame(Group = c(names(gain0.8), names(gain0.9)), 
                   Slope = c(coefficients(gain0.8)$Generation, coefficients(gain0.9)$Generation), 
                   Cor = rep(c(0.8, 0.9), each = length(names(gain0.8)))
                   )

head(gain)
gain <- gain %>%  separate(Group, into = c("Scenario",  "Import","Rep"), sep="-")
head(gain)
gainAvg <- gain %>% group_by(Cor, Scenario, Import) %>% dplyr::summarize(meanSlope = mean(Slope), sd = sd(Slope))
head(gainAvg)
table(gainAvg$Import)
gainAvg$Import <- factor(gainAvg$Import, c("home_0", "home_10",   "home_25", "home_50", "home_100", "home_100BD", "import"))
gainAvg[gainAvg$Scenario == "GenGen" & gainAvg$Cor == 0.8,]
gainAvg[gainAvg$Scenario == "GenGen" & gainAvg$Cor == 0.9,]
head(gainAvg)



#For Class Gen you need two slopes, each for each scheme / era
head(TGV_oneTrait)
gain0.9_a <- lmList(zMean30 ~ Generation | regressionGroup, data=TGV_oneTrait[TGV_oneTrait$Trait == "cor = 0.9" & TGV_oneTrait$Generation %in% 31:40,])
gain0.8_a <- lmList(zMean30 ~ Generation | regressionGroup, data=TGV_oneTrait[TGV_oneTrait$Trait == "cor = 0.8" & TGV_oneTrait$Generation %in% 31:40,])
gain_a <- data.frame(Group = c(names(gain0.8_a), names(gain0.9_a)), 
                   Slope = c(coefficients(gain0.8_a)$Generation, coefficients(gain0.9_a)$Generation), 
                   Cor = rep(c(0.8, 0.9), each = length(names(gain0.8_a)))
)

# Second - genomics - era
head(TGV_oneTrait)
gain0.9_b <- lmList(zMean30 ~ Generation | regressionGroup, data=TGV_oneTrait[TGV_oneTrait$Trait == "cor = 0.9" & TGV_oneTrait$Generation %in% 41:50,])
gain0.8_b <- lmList(zMean30 ~ Generation | regressionGroup, data=TGV_oneTrait[TGV_oneTrait$Trait == "cor = 0.8" & TGV_oneTrait$Generation %in% 41:50,])
gain_b <- data.frame(Group = c(names(gain0.8_b), names(gain0.9_b)), 
                   Slope = c(coefficients(gain0.8_b)$Generation, coefficients(gain0.9_b)$Generation), 
                   Cor = rep(c(0.8, 0.9), each = length(names(gain0.8_b)))
)

gain_a <- gain_a %>%  separate(Group, into = c("Scenario",  "Import","Rep"), sep="-")
gainAvg_a <- gain_a %>% group_by(Cor, Scenario, Import) %>% dplyr::summarize(meanSlope = mean(Slope), sd = sd(Slope))
gain_b <- gain_b %>%  separate(Group, into = c("Scenario",  "Import","Rep"), sep="-")
gainAvg_b <- gain_b %>% group_by(Cor, Scenario, Import) %>% dplyr::summarize(meanSlope = mean(Slope), sd = sd(Slope))

gainAvg_a[gainAvg_a$Scenario == "ClassGen" & gainAvg$Cor == 0.8,]
gainAvg_a[gainAvg_a$Scenario == "ClassGen" & gainAvg$Cor == 0.9,]
gainAvg_b[gainAvg_a$Scenario == "ClassGen" & gainAvg$Cor == 0.8,]
gainAvg_b[gainAvg_a$Scenario == "ClassGen" & gainAvg$Cor == 0.9,]



# gain$Cor <- as.factor(gain$Cor)
# gain$Scenario <- vapply(strsplit(gain$Group,"_"), `[`, 1, FUN.VALUE=character(1))
# gain$Population <- vapply(strsplit(gain$Group,"_"), `[`, 2, FUN.VALUE=character(1))
# gain$Import <- vapply(strsplit(gain$Group,"_"), `[`, 3, FUN.VALUE=character(1))
# gain$Scenario <- plyr::revalue(gain$Scenario, c("ClassGen" = "10y_delay", "GenGen" = "No_delay"))

gainAvg$meanSlope <- as.numeric(gainAvg$meanSlope)
gainAvg$sd <- as.numeric(gainAvg$sd)
gainAvg$Cor <- as.factor(gainAvg$Cor)
gainAvg$Era <- "No delay"

# Regerssion slopes for ClassGen - two era
gainAvg_a$Era <- "Delay - PT"
gainAvg_b$Era <- "Delay - GS"
gainAvgC <- rbind(gainAvg_a[gainAvg_a$Scenario == "ClassGen", ], gainAvg_b[gainAvg_b$Scenario == "ClassGen", ])
gainAvgC$sd <- as.numeric(gainAvgC$sd)
gainAvgC$Cor <- as.factor(gainAvgC$Cor)
gainAvgC <- rbind(gainAvgC, gainAvg[gainAvg$Scenario == "GenGen",])
gainAvgC$meanSlope <- as.numeric(gainAvgC$meanSlope)
gainAvgC$sd <- as.numeric(gainAvgC$sd)
gainAvgC$Cor <- as.factor(gainAvgC$Cor)
gainAvgC$Import <- factor(gainAvgC$Import, c("home_0", "home_10","home_25", "home_50", "home_100", "home_100BD", "import"))

table(gainAvgC$Scenario, gainAvgC$Era)
gainAvgC$Era <- factor(gainAvgC$Era, c("No delay", "Delay - PT", "Delay - GS"))
ggplot(data = gainAvgC, aes(x=Import, y = meanSlope, group=Cor, fill=Cor)) + 
  geom_bar(stat="identity", position="dodge") + 
  geom_errorbar(aes(ymin = meanSlope - sd, ymax = meanSlope + sd), position ="dodge", size=0.8) + 
  scale_fill_manual("Correlation", values = viridis::plasma(8)[c(3, 6)]) + 
  theme_bw(base_size=16) + 
  facet_grid(cols=vars(Era)) 

library(Rmisc)
multiplot(slopeGenGen, slopeClassGen, cols=2)


# Plot genetic variance
head(TGVAvg_oneTrait[TGVAvg_oneTrait$Generation > 29,])
ggplot(data=TGVAvg_oneTrait[TGVAvg_oneTrait$Generation > 29,], aes(x=Generation, y=zSd30, group=PlotGroup, colour=PlotGroup)) + 
  geom_line(size=0.9) + 
  theme_bw(base_size = 18) + ylab("Genetic variance") + facet_grid(cols = vars(Trait), rows = vars(Scenario)) + 
  #scale_y_continuous(breaks = c(seq(3, 9, 0.5))) +
 # geom_ribbon(aes(ymin = sdG - sdsdG, ymax = sdG + sdsdG, fill = PlotGroup),  linetype = 0, alpha = 0.3) + 
  theme(legend.position = "right") + 
  scale_colour_manual(values = c(viridis::plasma(18)[rev(seq(1, 21, 3))][2:7], 'black')) 


# Plot genic variance
summary(TGVAvg_oneTrait$meanGenicSD)
summary(TGVAvg_oneTrait$meanzGenicSD)
ggplot(data=TGVAvg_oneTrait[TGVAvg_oneTrait$Generation > 29,], aes(x=Generation, y=meanGenicSD, group=PlotGroup, colour=PlotGroup)) + 
  geom_line(size=1) + 
  theme_bw(base_size = 18) + ylab("Genetic variance") + facet_grid(cols = vars(Trait), rows = vars(Scenario)) + 
  #scale_y_continuous(breaks = c(seq(3, 9, 0.5))) +
  #geom_ribbon(aes(ymin = sdG - sdsdG, ymax = sdG + sdsdG, fill = PlotGroup),  linetype = 0, alpha = 0.3) + 
  theme(legend.position = "right") + 
  scale_colour_manual(values = c(viridis::plasma(18)[rev(seq(1, 21, 3))][2:7], 'black')) 


# Plot genetic gain for just the one correlated trait for generations 30:50
ggplot(data=TGVAvg_oneTrait[TGVAvg_oneTrait$Generation > 39,], aes(x=Generation, y=meanZTGV, group=Group, colour=Group)) + geom_line() + 
  theme_bw(base_size = 18) + ylab("Genetic gain") + facet_grid(cols = vars(Import)) + 
  scale_y_continuous(breaks = c(seq(3, 9, 0.5))) +
  #geom_vline(xintercept = 50) +
  geom_hline(yintercept = TGVAvg_oneTrait$meanZTGV[TGVAvg_oneTrait$Import == "100_100" & TGVAvg_oneTrait$Generation == 60 & TGVAvg_oneTrait$Group == "home"]) + 
  geom_ribbon(aes(ymin = meanZTGV - sdZTGV, ymax = meanZTGV + sdZTGV, fill = Group),  linetype = 0, alpha = 0.3)

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
# acc <- read.csv("Accuracies.csv")
# head(acc)
# ggplot(data = acc[acc$Strategy == "GenGen",], aes(x=Gen, y=Cor, group=Cat, colour=Cat)) + 
#         geom_line() + 
#         facet_grid(cols=c(vars(Scenario)), rows = vars(Group))

#Only last generation accuracies

acc <- read.csv("Results/AccuraciesPost_rep.csv")
head(acc)

accA <- acc %>%  group_by(Category, Generation, Group, Scenario, Strategy, Trait) %>% dplyr::summarise(meanAcc = mean(Cor), sdAcc=sd(Cor))
head(accA)
ggplot(data = accA[!accA$Category %in% c("Male candidates1", "izl") & accA$Strategy == "GenGen",], aes(x=Generation, y=meanAcc, group=Category, colour=Category)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = meanAcc - sdAcc, ymax = meanAcc + sdAcc, fill = Category),  linetype = 0, alpha = 0.3)+
  facet_grid(rows = vars(Scenario), cols=vars(Group, Trait)) + 
  theme_bw(base_size=16)
#Compare home and import
ggplot(data = accA[!accA$Category %in% c("Male candidates1", "izl") & accA$Strategy == "GenGen",], aes(x=Generation, y=meanAcc, group=Group, colour=Group)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = meanAcc - sdAcc, ymax = meanAcc + sdAcc, fill = Group),  linetype = 0, alpha = 0.3)+
  facet_grid(rows = vars(Scenario), cols=vars(Category, Trait)) + 
  theme_bw(base_size=16)

# For delay scenarios
ggplot(data = accA[!accA$Category %in% c("Male candidates1", "izl") & accA$Strategy == "ClassGen",], aes(x=Generation, y=meanAcc, group=Category, colour=Category)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = meanAcc - sdAcc, ymax = meanAcc + sdAcc, fill = Category),  linetype = 0, alpha = 0.3)+
  facet_grid(rows = vars(Scenario), cols=vars(Group, Trait)) + 
  theme_bw(base_size=16)
#Compare home and import
ggplot(data = accA[!accA$Category %in% c("Male candidates1", "izl") & accA$Strategy == "ClassGen",], aes(x=Generation, y=meanAcc, group=Group, colour=Group)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = meanAcc - sdAcc, ymax = meanAcc + sdAcc, fill = Group),  linetype = 0, alpha = 0.3)+
  facet_grid(rows = vars(Scenario), cols=vars(Category, Trait)) + 
  theme_bw(base_size=16)

#Compare home and import
ggplot(data = accA[!accA$Category %in% c("Male candidates1", "izl") & accA$Trait == 2,], aes(x=Generation, y=meanAcc, group=Strategy, colour=Strategy)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = meanAcc - sdAcc, ymax = meanAcc + sdAcc, fill = Strategy),  linetype = 0, alpha = 0.3)+
  facet_grid(rows = vars(Scenario), cols=vars(Category, Group)) + 
  theme_bw(base_size=16)



# Now compare the one value for accuracy
# But separately for butn-in, first and second ten years
head(acc)
acc$Era <- ifelse(acc$Generation %in% 21:30, "Burn-in", ifelse(acc$Generation %in% 31:40, "FirstTenY", "SecondTenY"))
accAE <- acc %>%  group_by(Category, Era, Group, Scenario, Strategy, Trait) %>% dplyr::summarise(meanAcc = mean(Cor), sdAcc=sd(Cor))

#BurnIN
ggplot(data = accAE[!accAE$Category %in% c("Male candidates1", "izl") & accAE$Strategy == "GenGen" & accAE$Era =="Burn-in",], aes(x=Category, y=meanAcc, group=Group, fill=Group)) + 
  geom_bar(stat="identity", position = "dodge") + 
  facet_grid(rows = vars(Scenario), cols=vars(Trait)) + 
  theme_bw(base_size=16)

# First ten years
ggplot(data = accAE[!accAE$Category %in% c("Male candidates1", "izl") & accAE$Strategy == "GenGen" & accAE$Era =="FirstTenY",], aes(x=Category, y=meanAcc, group=Group, fill=Group)) + 
  geom_bar(stat="identity", position = "dodge") + 
  facet_grid(rows = vars(Scenario), cols=vars(Trait)) + 
  theme_bw(base_size=16)
# Second ten years
ggplot(data = accAE[!accAE$Category %in% c("Male candidates1", "izl") & accAE$Strategy == "GenGen" & accAE$Era =="SecondTenY",], aes(x=Category, y=meanAcc, group=Group, fill=Group)) + 
  geom_bar(stat="identity", position = "dodge") + 
  facet_grid(rows = vars(Scenario), cols=vars(Trait)) + 
  theme_bw(base_size=16)

#CLASS-GEN
#BurnIN
ggplot(data = accAE[!accAE$Category %in% c("Male candidates1", "izl") & accAE$Strategy == "ClassGen" & accAE$Era =="Burn-in",], aes(x=Category, y=meanAcc, group=Group, fill=Group)) + 
  geom_bar(stat="identity", position = "dodge") + 
  facet_grid(rows = vars(Scenario), cols=vars(Trait)) + 
  theme_bw(base_size=16)

# First ten years
ggplot(data = accAE[!accAE$Category %in% c("Male candidates1", "izl") & accAE$Strategy == "ClassGen" & accAE$Era =="FirstTenY",], aes(x=Category, y=meanAcc, group=Group, fill=Group)) + 
  geom_bar(stat="identity", position = "dodge") + 
  facet_grid(rows = vars(Scenario), cols=vars(Trait)) + 
  theme_bw(base_size=16)
# Second ten years
ggplot(data = accAE[!accAE$Category %in% c("Male candidates1", "izl") & accAE$Strategy == "ClassGen" & accAE$Era =="SecondTenY",], aes(x=Category, y=meanAcc, group=Group, fill=Group)) + 
  geom_bar(stat="identity", position = "dodge") + 
  facet_grid(rows = vars(Scenario), cols=vars(Trait)) + 
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
setwd("/home/jana/EddieDir/10K/SU55_import/GenGen0_10_1013/")
gh <- read.csv("GenPed_EBVhome.txt")
gh$Group <- "Home"
# For import population
gi <- read.csv("GenPed_EBVimport.txt")
gi$Group <- "Import"
# Bind them together
g <- rbind(gh, gi)
#Aggregate EBVs by Generation and Group 
gs <- g %>% group_by(Generation, Group) %>% summarise(meanEBV = mean(EBV), meanTGV1 = mean(gvNormUnres1), meanTGV2 = mean(gvNormUnres2)) %>% 
  pivot_longer(cols = c(meanEBV, meanTGV1, meanTGV2))
#Plot the EBVs
ggplot(data = gs, aes(x = Generation, y = value, group = name, colour = name)) + geom_line() + facet_grid(rows = vars(Group))



#Scenario
setwd("/home/jana/EddieDir/10K/SU55_import/GenGen0_10_1012/")
ped <- read.table("SimulatedData/PedigreeAndGeneticValues_cat.txt", header=T, sep=" ")
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
PEDs <- PED

library(reshape)
PEDs <- PEDs[,c("Indiv", "Generation", "Group", "gvNormUnres1", "gvNormUnres2", "gvNormUnres3")]
pedm <- melt(PEDs, id.vars = c("Indiv", "Generation", "Group"))

pedma <- pedm %>% group_by(Generation, Group, variable) %>% summarize(Gain = mean(value), GainSD = sd(value))

pedma <- pedma[pedma$variable %in% c("gvNormUnres1", "gvNormUnres2", "gvNormUnres3"),]
pedma$variable <- plyr::revalue(pedma$variable, c("gvNormUnres1" = "Trait 1", "gvNormUnres2" = "Trait 2", "gvNormUnres3" = "Trait 3"))
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
gi <- read.csv("GENINTS_all_28102020.csv")
gi$lineSex <- paste0(gi$line, gi$sex)
head(gi)
table(gi$Gen)

ggplot(data=gi[gi$strategy == "GenGen",], aes(x=Gen, y=genInt, group=lineSex, colour=lineSex)) + geom_line() +
  facet_grid(cols=c(vars(Group), vars(trait)), rows=vars(scenario))

ggplot(data=gi[gi$strategy == "ClassGen",], aes(x=Gen, y=genInt, group=lineSex, colour=lineSex)) + geom_line() +
  facet_grid(cols=c(vars(Group), vars(trait)), rows=vars(scenario))



#######################################
## Partitioned trends
#######################################
part <- read.csv("Results/Partition_Import_28102020.csv", sep= " ")
head(part)
#part <- read.csv("Partition_Import.csv", sep= " ")
part$Import <- plyr::revalue(part$Import, c("0_0" = "0", "10_10" = "10", "50_50" = "50", "100_100" = "100", "0_100" = "100BD"))
part$Import <- factor(part$Import, c("0", "10", "50", "100", "100BD"))
plotTGV$PlotGroup <- factor(plotTGV$PlotGroup, c("home_0", "home_10", "home_50", "home_100", "home_100BD", "import"))

partL <- part %>% pivot_longer(cols = c(home, import, Sum), names_to = "Partition" )

ggplot(data = partL[partL$Population == "home" & partL$Scenario == "GenGen",], aes(x=Generation, y = value, group = Partition, colour = Partition)) + 
  geom_line(size=1) + 
  facet_grid(cols = vars(Import), rows = vars(Cor)) + theme_bw(base_size=16) + 
  scale_colour_manual("", values = viridis::plasma(8)[c(1,4,6)])

ggplot(data = partL[partL$Population == "home" & partL$Scenario == "ClassGen",], aes(x=Generation, y = value, group = Partition, colour = Partition)) + 
  geom_line(size=1) + 
  facet_grid(cols = vars(Import), rows = vars(Cor)) + theme_bw(base_size=16) + 
  scale_colour_manual("", values = viridis::plasma(8)[c(6,4,1)])


#######################################
## Check fathers TGVs
#######################################
allMales <- read.csv("Results/AllBulls_TGVs_raw.csv")
homeMales <- read.csv("Results/HomeFathers_TGVs_raw.csv")

# All males selecte
allMales_m <- allMales %>% pivot_longer(c(gvNormUnres1, gvNormUnres3))
# Check the by-trait comparison
ggplot(data = filter(allMales_m, strategy == "GenGen"), aes(x = Generation, y = value, group = Group, colour = Group)) + 
  geom_point() + 
  facet_grid(cols = vars(scenario), rows = vars(name))

# Check only for the one observed trait
allMales_oneTrait <- allMales_m %>% filter((name == "gvNormUnres1" & Group == "home") | (name == "gvNormUnres3" & Group == "import"))
head(allMales_oneTrait)
table(allMales_oneTrait$Group)
ggplot(data = filter(allMales_oneTrait, strategy == "GenGen"), aes(x = Generation, y = value, group = Group, colour = Group)) + 
  geom_point() + 
  facet_grid(cols = vars(scenario))


# Check the same only for the males used in home population
homeMales <- read.csv("Results/HomeFathers_TGVs_raw.csv")

# Males used in home population
homeMales_m <- homeMales %>% pivot_longer(c(gvNormUnres1, gvNormUnres3))
# Check the by-trait comparison
ggplot(data = filter(homeMales_m, strategy == "GenGen"), aes(x = Generation, y = value, group = Group, colour = Group)) + 
  geom_point() + 
  facet_grid(cols = vars(scenario), rows = vars(name))

# Check only for the one observed trait
homeMales_oneTrait <- homeMales_m %>% filter((name == "gvNormUnres1" & Group == "home") | (name == "gvNormUnres3" & Group == "import"))
head(homeMales_oneTrait)
table(homeMales_oneTrait$Group)
ggplot(data = filter(homeMales_oneTrait, strategy == "GenGen"), aes(x = Generation, y = value, group = Group, colour = Group)) + 
  geom_point() + 
  facet_grid(cols = vars(scenario))



######################################
# Reference size
#####################################
ref <- read.csv("Results/ReferenceSize_Import_15122020.csv")[,-1]
head(ref)
refA <- ref %>% group_by(Group, scenario, strategy, trait) %>% summarise(meanSize = mean(V1))
refA %>%  ggplot(aes(x = scenario, y = meanSize, group = Group,fill = Group)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_grid(rows = vars(strategy), cols = vars(trait))
