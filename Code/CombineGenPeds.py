import pandas as pd

scenarios = ["10_10", "50_50", "100_100", "0_100"]
rep = [4]

# Initiate the dataframe for TGVs for all bulls and bulls used in the home population
EBVs = pd.DataFrame()

# Check the results for each rep and each scenario within
for rep in rep:
    for scenario in scenarios:
        pedI = pd.read_csv("GenGen" + str(rep) + "_" + scenario + "13/SimulatedData/GenPed_EBVhome.txt")
        pedH = pd.read_csv("GenGen" + str(rep) + "_" + scenario + "13/SimulatedData/GenPed_EBVimport.txt")
        pedH.loc[:, "Group"] = "home"
        pedI.loc[:, "Group"] = "import"

        ped = pedH.append(pedI)
        pedMean = ped.groupby(['Generation', 'Group'])[['EBV', 'gvNormUnres1', 'gvNormUnres3']].mean()
        pedMean.loc[:, "scenario"] = scenario
        pedMean.loc[:, "rep"] = rep
        EBVs = EBVs.append(pedMean)

EBVs.groupby(['scenario', 'Generation', 'Group'])[['EBV', 'gvNormUnres1', 'gvNormUnres3']].mean().to_csv("EBVs.csv")