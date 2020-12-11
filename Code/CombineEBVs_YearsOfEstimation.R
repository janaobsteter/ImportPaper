###################################
# This is a script to combine the EBVs from multiple years of estimation
###################################

library(dplyr)

# Create an empty tibble
meanEBVs <- tibble()

# For each generation read in the EBVs, add generation and year and bind them to the meanEBVs tibble
for (gen in 21:40) {
  # Read in solutions for gome and import populations
  rI <- read.table(paste0("renumbered_Solutions_import_", gen))[,-1]
  rH <- read.table(paste0("renumbered_Solutions_home_", gen))[,-1]
  
  # Read in the dataframes with ID and Generation Information
  gH <- read.csv("GenPed_EBVhome.txt")
  gI <- read.csv("GenPed_EBVimport.txt")
  
  # Rename the columns
  colnames(rH) <- c("Indiv", "EBVh")
  colnames(rI) <- c("Indiv", "EBVi")
  
  # Add the group to the EBVs, either home of import
  rH$Group <- ifelse(rH$Indiv %in% gH$Indiv, "home", "import")
  rI$Group <- ifelse(rH$Indiv %in% gH$Indiv, "home", "import")
  
  # Merge the dataframes, add the Generation of the individuals and year of the estimation
  r <- merge(rH, rI, by=c("Indiv", "Group"))
  r$Generation <- rep(1:gen, each=8640*2)
  r$YearEBV <- gen
  
  # Bind to meanEBVs
  meanEBVs <- bind_rows(meanEBVs,
                    r %>% group_by(Generation, Group, YearEBV) %>% summarise(EBVi = mean(EBVi), EBVh = mean(EBVh)))
  
  
}

write.csv(as.data.frame(meanEBVs), "MeanEBVs.txt", quote=F, row.names=F)