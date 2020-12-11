import pandas as pd

scenarios = ["10_10", "50_50", "100_100", "0_100"]
rep = [4]

catDict = {
  "genTest1" : "Male candidates",
  "telF1" : "Female candidates",
  "k3" : "Mothers",
  "k4" : "Mothers",
  "k5" : "Mothers",
  "k6" : "Mothers",
  "pBM4" : "Mothers",
  "pBM5" : "Mothers",
  "pBM6" : "Mothers",
  "bm7" : "Mothers",
  "gpb2" : "Fathers",
  "gpb3" : "Fathers",
  "gpb4" : "Fathers",
  "gpb5" : "Fathers",
  "gpb6" : "Fathers",
  "pb6" : "Fathers",
  "pb7" : "Fathers",
  "pb8" : "Fathers",
  "pb9" : "Fathers",
  "vhlevljeni1" : "Male candidates1",
  "cak5" : "Male candidates"}

# Initiate the dataframe for TGVs for all bulls and bulls used in the home population
acc = pd.DataFrame()



# Check the results for each rep and each scenario within
for rep in rep:
    for scenario in scenarios:
        for group in ['home', 'import']:
            accRep = pd.read_csv("GenGen" + str(rep) + "_" + scenario + "13/Accuracies_CatAgehome.csv")
            accRep = accRep[accRep.AgeCat.isin(catDict.keys())]
            accRep.loc[:, "Cat"] = [catDict[x] for x in list(accRep.AgeCat)]
            accMean = accRep.groupby(["Cat", "Gen"])[['Cor']].mean()
            accRep.loc[:, "Group"] =  group
            accRep.loc[:, "Scenario"] = scenario
            accRep.loc[:, "Rep"] = rep
            # Concatenate the dataframes
            acc = acc.append(accRep)


# Write the files
acc.groupby(['Scenario', 'Group', 'Gen', 'Cat'])[['Cor']].mean().to_csv("Accuracies.csv")
