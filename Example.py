import gurobipy as gp
import numpy as np
import pickle

from src import *
from src.CreateEColiComm import *
from src.CreateCommModel import CreateCommModel
from src.parseXML import parseXML
from src.FindMaxGrowthRate import FindMaxGrowthRate
from src.RemoveBlocked import RemoveBlocked
from src.FluxCouplingAnalysis import FluxCouplingAnalysis

MetabolicMdl = parseXML('iAF1260.xml')
#Adding reaction 'met__L_c' + 'h_c' --> 'met__L_p' + 'h_p'
MetabolicMdl[0] = np.hstack((MetabolicMdl[0],np.zeros( (MetabolicMdl[0].shape[0],1) )))
MetabolicMdl[0][MetabolicMdl[3].index('met__L_c')][MetabolicMdl[0].shape[1] - 1] = MetabolicMdl[0][MetabolicMdl[3].index('h_c')][MetabolicMdl[0].shape[1] - 1] = -1
MetabolicMdl[0][MetabolicMdl[3].index('met__L_p')][MetabolicMdl[0].shape[1] - 1] = MetabolicMdl[0][MetabolicMdl[3].index('h_p')][MetabolicMdl[0].shape[1] - 1] = 1

#Replacing the ±999999 with ±1000 lets this code work, probably because of getting rid of large coefficients, which may cause numerical instability
for Idx in range(len(MetabolicMdl[1])):
    if MetabolicMdl[1][Idx] < -1000:
        MetabolicMdl[1][Idx] = -1000
    if MetabolicMdl[2][Idx] > 1000:
        MetabolicMdl[2][Idx] = 1000

MetabolicMdl[4].append('METt3pp')#Add newly added reactions' name
MetabolicMdl[5].append('METt3pp')
#Add new reactions' bounds
MetabolicMdl[1] = np.append(MetabolicMdl[1],0)
MetabolicMdl[2] = np.append(MetabolicMdl[2],1000)

TheListOfMdls = TwoComm(MetabolicMdl) #Load E. Coli consortia
CommMdl = CreateCommModel(TheListOfMdls,"EX")
print(CommMdl[0].shape,*[len(CommMdl[Num]) for Num in range(1,6)])
#For this example, set the community uptake/export reaction bounds to equal that of the original model's exchange reactions'
for Idx,Rxn in enumerate(MetabolicMdl[4]):
    if Rxn.startswith("EX"):
        TheIdx = CommMdl[4].index(Rxn + "_u")
        CommMdl[1][TheIdx] = MetabolicMdl[1][Idx]
        CommMdl[2][TheIdx] = MetabolicMdl[2][Idx]
GrowthRate,TheMdl,BioIndices = FindMaxGrowthRate(CommMdl,"BIOMASS")

TrimmedLP,TrimmedComm,TheKVT = RemoveBlocked(TheMdl,CommMdl,GrowthRate*.99,BioIndices)#Using the exact value of GrowthRate may lead to an infeasible model. Use one that's slightly less
# TrimmedLP.write("TrimmedLP.mps")
with open("TrimmedData.pkl",'wb') as File:
    pickle.dump((TrimmedComm,TheKVT),File)
with open("FCAResults.pkl","wb") as File:
    pickle.dump(FluxCouplingAnalysis(TrimmedLP,TheKVT),File)