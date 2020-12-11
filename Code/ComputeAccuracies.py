import pandas as pd

scenarios = ["10_10", "50_50", "100_100", "0_100"]
rep = [0,1,2,4]


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
  "cak5" : "Male candidates",
"izl59": "izl"}

accMean = pd.DataFrame()
# Check the results for each rep and each scenario within
for rep in rep:
    for scenario in scenarios:
        PED = pd.read_csv("GenGen" + str(rep) + "_" + scenario + "13/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")
        ps = pd.read_csv("GenGen" + str(rep) + "_" + scenario + "13/PopulationSplit.txt", names = ['Group', 'Indiv'], header=1, low_memory=False)
        # Merge the files
        PED = PED.merge(ps, on="Indiv")
        for group in ['home', 'import']:
            
            for gen in range(21, 61):
                ped = PED[PED.Group == group]
                print("Generation: " + str(gen))
                sol = pd.read_table("GenGen" + str(rep) + "_" + scenario + "13/renumbered_Solutions_" + group + "_" + str(gen), names=['renID', 'Indiv', 'EBV'+str(gen)], sep=" ")
                cat = pd.read_csv("GenGen" + str(rep) + "_" + scenario + "13/Categories_gen" + str(gen) + "DF_" + (str(0) if group == "home" else str(1)) + ".csv")
                catm = pd.melt(cat)
                catm = catm.dropna()
                catm.value = catm.value.astype(int)
                catm.columns = ['Cat' + str(gen), 'Indiv']

                ped = ped.merge(sol[['Indiv', 'EBV' + str(gen)]], on="Indiv", how="left")
                ped = ped.merge(catm, on="Indiv", how="left")

                ped.loc[:, "Age" + str(gen)] = gen - ped.Generation
                
                ped.loc[:, "CatAge" + str(gen)] = ped.loc[:, "Cat" + str(gen)] + ped.loc[:, "Age" + str(gen)].astype(str)
                
                ped = ped[ped['CatAge'+str(gen)].isin(catDict.keys())]
                ped.loc[:, "Category" + str(gen)] = [catDict[x] for x in list(ped.loc[:, 'CatAge' + str(gen)])]
                
                #Compute
                if group == "home":
                    pedMean = pd.DataFrame(ped.groupby(["Category" + str(gen)])[['EBV' + str(gen), 'gvNormUnres1']].corr().iloc[0::2,-1])
                    pedMean = pedMean.reset_index()
                    pedMean.columns = ["Category", "EBVGen", "Cor"]
                elif group == "import":
                    pedMean = pd.DataFrame(ped.groupby(["Category" + str(gen)])[['EBV' + str(gen), 'gvNormUnres3']].corr().iloc[0::2,-1])
                    pedMean = pedMean.reset_index()
                    pedMean.columns = ["Category", "EBVGen", "Cor"]
                pedMean.loc[:, "Group"] = group
                pedMean.loc[:, "Scenario"] = scenario
                pedMean.loc[:, "Rep"] = rep
                pedMean.loc[:, "Generation"] = gen
                
                # Concatenate the dataframes
                accMean = accMean.append(pedMean)


accMean.groupby(['Group', 'Scenario', 'Generation', 'Category'])[['Cor']].mean().to_csv("AccuraciesPost_Rep0.csv")

import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sns.set_theme(style="darkgrid")
g = sns.FacetGrid(accMean, col="Scenario", row="Group", hue="Category")
g.map(sns.lineplot, "Generation", "Cor")
plt.show()