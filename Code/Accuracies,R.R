setwd("~/EddieDir/10K/SU55_import/GenGen0_0_012/")
p <- read.table("SimulatedData/PedigreeAndGeneticValues_cat.txt", header=T, sep=" ")
ps <- read.csv("PopulationSplit.txt")
colnames(ps) <- c("Group", "Indiv")
p <- merge(p, ps, by="Indiv")
table(p$Generation[p$cat == "pb"])
table(p$Generation[p$Indiv %in% p$Father[p$cat == "nr"]])

head(p)
library(reshape)
library(dplyr)
library(ggplot2)
pA <- p %>% group_by(Generation, Group) %>%  summarise(mean1 = mean(gvNormUnres1), mean3 = mean(gvNormUnres3)) 
pAm <- melt(as.data.frame(pA), id.vars = c("Generation", "Group"))
head(pAm)
pAm_one <- rbind(pAm[pAm$Group =="home" & pAm$variable == "mean1",], pAm[pAm$Group =="import" & pAm$variable == "mean3",])
ggplot(data=pAm_one, aes(x=Generation, y=value, group=Group, colour=Group)) + geom_line()


acc <- data.frame()

for (rep in c(0, 1, 2, 4)) {
  for (scenario in c("0_0", "10_10", "50_50", "100_100", "0_100")) {

    p <- read.table(paste0("GenGen", rep, "_", scenario, "13/SimulatedData/PedigreeAndGeneticValues_cat.txt"), header=T, sep=" ")
    ps <- read.csv(paste0("GenGen", rep, "_", scenario, "13/PopulationSplit.txt"))
    colnames(ps) <- c("Group", "Indiv")
    p <- merge(p, ps, by="Indiv")
    p$CatAge <- paste0(p$cat, p$age)
    
    ri <- read.table("renumbered_Solutions_import_60")
    colnames(ri) <- c("renID", "Indiv", "EBVi")
    p <- merge(p, ri, by="Indiv", all.x=T)
    rh <- read.table("renumbered_Solutions_home_60")
    colnames(rh) <- c("renID", "Indiv", "EBVh")
    p <- merge(p, rh, by="Indiv", all.x=T)
    
    
    
    PED <- p[p$CatAge %in% c("genTest1", "telF1",  "k3", "k4", "k5", "k6", "pBM4", "pBM5", "pBM6", "bm7", "gpb2", "gpb3", "gpb4", "gpb5", "gpb6", "pb6", "pb7", "pb8", "pbv9", "vhlevljeni1", "cak5"),]
    
    
    PED$CatAge <- plyr::revalue(PED$CatAge , c("genTest1" = "Male candidates", "telF1" = "Female candidates", 
                                                "k3" = "Mothers", 
                                                "k4" = "Mothers", 
                                                "k5" = "Mothers",
                                                "k6" = "Mothers", 
                                                "pBM4" = "Mothers", 
                                                "pBM5" = "Mothers", 
                                                "pBM6" = "Mothers", 
                                                "bm7" = "Mothers", 
                                                "gpb2" = "Fathers", 
                                                "gpb3" = "Fathers", 
                                                "gpb4" = "Fathers", 
                                                "gpb5" = "Fathers", 
                                                "gpb6" = "Fathers", 
                                                "pb6" = "Fathers", 
                                                "pb7" = "Fathers", 
                                                "pb8" = "Fathers", 
                                                "pb9" = "Fathers", 
                                                "vhlevljeni1" = "Male candidates1",
                                                "cak5" = "Male candidates"))
    
    PED0 <- PED[PED$Generation > 39,]
    accTmp <- PED0 %>% group_by(CatAge, Group) %>% dplyr::summarise(accEBVh = cor(gvNormUnres1, EBVh), accEBVi = cor(gvNormUnres3, EBVi))
    accTmp$Scenario <- scenario
    accTmp$rep <- rep
    acc <- rbind(acc, accTmp)
  }
}

write.csv(acc, "Acc_LastGen.csv", quote=F, row.names=F)

