## this is a script to combine the heritability partitions
######################################################
#EBV
######################################################
library(AlphaPart)
library(mltools)


for (rep in 0:5) {
  for (strategy in c("GenGen", "ClassGen")) {
    for (scenario in c("0_0", "10_10", "50_50", "100_100", "0_100")) { 
      # Read in the pedigree
      Ped <- read.table(paste0(strategy, str(rep), "_", scenario, "13/SimulatedData/PedigreeAndGeneticValues_cat.txt"), header=TRUE, sep=" ")
      ps <- read.csv("PopulationSplit.txt")
      colnames(ps) <- c("Group", "Indiv")
      Ped <- merge(Ped, ps, by="Indiv")
      Ped <- Ped[Ped$Generation %in% 20:60, c("Generation", "IId", "FId", "MId", "sex", "Group", "gvNormUnres1", "gvNormUnre3")]

      ########################################################
      #TGV
      #######################################################
    
      PedEval1$TbvT1_s <- (PedEval1$TbvT1 -  mean(PedEval1$TbvT1[PedEval1$Generation == 20])) / sd(PedEval1$TbvT1[PedEval1$Generation == 20])
      PedEval1$TbvT2_s <- (PedEval1$TbvT2 -  mean(PedEval1$TbvT2[PedEval1$Generation == 20])) / sd(PedEval1$TbvT2[PedEval1$Generation == 20])
      PedEval1$TbvI_s <- 0.5 * (PedEval1$TbvT1_s + PedEval1$TbvT2_s)
      

      Ped$PopulationGender = paste(Ped$Group, Ped$sex, sep = "_")

            print("Partition genetic values")
      
      Part = AlphaPart(x = as.data.frame(Ped), sort = FALSE,
                        colId = "IId", colFid = "FId", colMid = "MId",
                        colPath = "PopulationGender", colAGV = c("gvNormUnres1", "gvNormUnres3"))

      
      # Cut the "home" animals outs
      # TGVs
      
      PartHome <- Part
      PartHome$gvNormUnres1 <- PartHome$gvNormUnres1[PartHome$gvNormUnres1$Group == "home",]
      PartHome$gvNormUnres3 <- PartHome$gvNormUnres3[PartHome$gvNormUnres3$Group == "home",]
      
      PartImport <- Part
      PartImport$gvNormUnres1 <- PartImport$gvNormUnres1[PartImport$gvNormUnres1$Group == "import",]
      PartImport$gvNormUnres3 <- PartImport$gvNormUnres3[PartImport$gvNormUnres3$Group == "import",]
      

      # Summarise partitions by Group (home \ import)
      PartHomeSummary = summary(object = PartHome, by = "Generation")
      PartImportSummary = summary(object = PartImport, by = "Generation")

      PartitionHomeDF <- data.frame()
      PartitionImportDF <- data.frame()
      for (trait in c("gvNormUnres1", "gvNormUnres3")) {
        # Home
        t1 <- PartHomeSummary[[trait]]$abs
        t1$Population = "home"
        t1$Trait <- trait
        t1$value <- "Tbv"
        PartitionHomeDF <- rbind(PartitionHomeDF, t1)
        
        # Import
        t1 <- PartImportSummary[[trait]]$abs
        t1$Population = "import"
        t1$Trait <- trait
        t1$value <- "Tbv"
        PartImportSummary <- rbind(PartImportSummary, t1)
        

        
      }
    }
    write.table(PartHomeSummary, paste0("PartitionHome_", rep, ".csv"), quote=FALSE, row.names=FALSE)
    write.table(PartImportSummary, paste0("PartitionImport_", rep, ".csv"), quote=FALSE, row.names=FALSE)
  }
    

  
write.table(PartHomeSummary, paste0("PartitionHome_.csv"), quote=FALSE, row.names=FALSE)
write.table(PartImportSummary, paste0("PartitionImport.csv"), quote=FALSE, row.names=FALSE)
