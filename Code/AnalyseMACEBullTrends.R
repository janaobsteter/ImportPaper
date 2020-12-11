# MACE Analysis
pvs <- read.csv("/home/jana/Documents/PhD/Projects/inProgress/GenomicStrategies_Import/MACE_Slo_Tuji_avgust2020.csv")
pvs$DAT_ROJSTVO1 <- as.Date(pvs$DAT_ROJSTVO, format="%d-%m-%Y")
head(pvs)
pvs$BY <- format(pvs$DAT_ROJSTVO1, "%Y")
pvs <- pvs[pvs$PASMA %in% 1:3,]
pvs <- pvs[pvs$DRZ_ORIG_ZIVAL %in% c("AT", "CH", "DE", "FR", "SI", "US", "FR"),]
pvs$PASMA <- as.character(pvs$PASMA)
pvs$PASMA <- plyr::revalue(pvs$PASMA, c("1" = "BSW", "2" = "SIM", "3" = "HOL"))

library(dplyr)
library(ggplot2)
pvsM <- pvs %>% group_by(BY, DRZ_ORIG_ZIVAL, PASMA) %>% summarize(MeanBV = mean(VREDNOST_12_PV), SdBV = sd(VREDNOST_12_PV))
pvsM$PASMA <- as.factor(pvsM$PASMA)
pvsM$BY <- as.numeric(pvsM$BY)
pvsM <- pvsM[pvsM$BY > 1980,]
pvsM$Group <- ifelse(pvsM$DRZ_ORIG_ZIVAL == "SI", "Slo", "Tuji")
#po državah
ggplot(data = pvsM, aes(x=BY, y=MeanBV, colour=DRZ_ORIG_ZIVAL, group=DRZ_ORIG_ZIVAL, fill=DRZ_ORIG_ZIVAL)) + 
  geom_ribbon(aes(ymin=(MeanBV - SdBV), ymax=(MeanBV + SdBV), fill=DRZ_ORIG_ZIVAL), alpha=0.1, colour=NA) +
  geom_line() + facet_grid(rows = vars(PASMA)) +
  theme_bw(base_size = 16)


pvs$Group <- ifelse(pvs$DRZ_ORIG_ZIVAL == "SI", "Slo", "Tuji")
pvsMG <- pvs %>% group_by(BY, Group, PASMA) %>% summarize(MeanBV = mean(VREDNOST_12_PV), SdBV = sd(VREDNOST_12_PV))
#Slo-tuji
ggplot(data = pvsMG, aes(x=BY, y=MeanBV, colour=Group, group=Group, fill=Group)) + 
  geom_ribbon(aes(ymin=(MeanBV - SdBV), ymax=(MeanBV + SdBV), fill=Group), alpha=0.1, colour=NA) +
  geom_line() + facet_grid(rows = vars(PASMA)) +
  theme_bw(base_size = 16)

ggplot(data = pvsM, aes(x=BY, y=MeanBV, colour=DRZ_ORIG_ZIVAL, group=DRZ_ORIG_ZIVAL, fill=DRZ_ORIG_ZIVAL)) + 
  geom_ribbon(aes(ymin=(MeanBV - SdBV), ymax=(MeanBV + SdBV), fill=DRZ_ORIG_ZIVAL), alpha=0.1, colour=NA) +
  geom_bar(stat="identity", position="dodge") + facet_grid(rows = vars(PASMA)) +
  theme_bw(base_size = 16)

head(pvsM)
as.data.frame(tail(pvsM, 50))

diff <- data.frame(Pasma = NA, Pop = NA, BY = NA, Diff = NA)
for (pasma in c("HOL", "BSW", "SIM")) {
  for (by in 1980:2016) {
      tmp <- pvsM[pvsM$PASMA == pasma & pvsM$BY == by,]
      for (pop in unique(pvsM$DRZ_ORIG_ZIVAL)) {
        diff <- rbind(diff, c(pasma, pop, by, 
                              ifelse(length(tmp$MeanBV[tmp$DRZ_ORIG_ZIVAL == pop] - tmp$MeanBV[tmp$DRZ_ORIG_ZIVAL == "SI"]) != 0, (tmp$MeanBV[tmp$DRZ_ORIG_ZIVAL == pop] - tmp$MeanBV[tmp$DRZ_ORIG_ZIVAL == "SI"]), 0)))
      
    }
    
  }
}

diff <- diff[!is.na(diff$Pasma),]
diff$Diff <- as.numeric(diff$Diff)
ggplot(data = diff, aes(x=BY, y=Diff, group=Pop, colour=Pop, fill=Pop)) + geom_bar(stat="identity", position="dodge") + 
  geom_hline(yintercept=12) + facet_grid(rows = vars(Pasma))

