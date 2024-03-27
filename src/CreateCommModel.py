import numpy as np
def CreateCommModel(TheModels, ExchangeRxnID):
    BoundsCommunityUptakeExport = {}
    for Mdl in TheModels:
        for Idx,Rxn in enumerate(Mdl[4]):
            if Rxn.split('_')[0] == ExchangeRxnID:
                try: #(Exchanged metabolite,[Lower bound, upper bound, number of occurences of rxn])
                    BoundsCommunityUptakeExport[Rxn][1][0]+=Mdl[1][Idx]
                    BoundsCommunityUptakeExport[Rxn][1][1]+=Mdl[2][Idx]
                    BoundsCommunityUptakeExport[Rxn][1][2]+=1
                except KeyError:
                    BoundsCommunityUptakeExport[Rxn] = (Mdl[3][np.nonzero(Mdl[0][:,Idx])[0][0]],[Mdl[1][Idx],Mdl[2][Idx],1])
    # print(BoundsCommunityUptakeExport)
    ExchangedMets = [Val[0] for Val in BoundsCommunityUptakeExport.values()]
    BoundsCommunityUptakeExport = dict((Key,(Val[1][0]/Val[1][2],Val[1][1]/Val[1][2])) for Key,Val in BoundsCommunityUptakeExport.items()) #Taking the average of the bounds

    NumOfExchangeRxns = len(BoundsCommunityUptakeExport)
    TotalNumOfRxns = sum(Mdl[0].shape[1] for Mdl in TheModels) + NumOfExchangeRxns

    MatrixForModels = [np.hstack((TheModels[0][0],np.zeros((TheModels[0][0].shape[0], TotalNumOfRxns - TheModels[0][0].shape[1]))))] + ([None] * len(TheModels)) #Number of entries is len(TheModels) + 1. + 1 is here for community uptake and export reactions

    #Merging metabolic models
    RxnsIncluded = TheModels[0][0].shape[1]
    for Idx in range(1,len(TheModels)):
        MatrixForModels[Idx] = np.hstack((np.zeros((TheModels[Idx][0].shape[0], RxnsIncluded) ),TheModels[Idx][0], np.zeros((TheModels[Idx][0].shape[0], TotalNumOfRxns - RxnsIncluded - TheModels[Idx][0].shape[1]))  ))
        RxnsIncluded+=TheModels[Idx][0].shape[1] #Number of rxns the community stoichiometric matrix contains from the community so far
    #First entry is [First model's S_Matrix, zeros], second entry is [Zeros, second model's s_matrix, zeros],...[Last model's S_Matrix, zeros]
    
    #Creating entires for uptake and export, then merging them
    MatrixForModels[-1] = np.hstack((np.zeros((NumOfExchangeRxns, RxnsIncluded)),-np.identity(NumOfExchangeRxns)))
    ColIdx = 0
    for Mdl in TheModels:
        RowIdx = 0
        for Idx in range(len(Mdl[4])):
            if Mdl[4][Idx].startswith("EX_"):
                MatrixForModels[-1][RowIdx][Idx + ColIdx] = 1 #The original exchange reactions now interact with the community uptake and export, turning them from "met[e] <=>" to "met[e] <=> met[u]"
                RowIdx = RowIdx + 1
        ColIdx = ColIdx + Mdl[0].shape[1]

    return [np.vstack(tuple(MatrixForModels)),\
        [LB for Mdl in TheModels for LB in Mdl[1]] + ([-float('inf')] * len(BoundsCommunityUptakeExport)),\
        [UB for Mdl in TheModels for UB in Mdl[2]] + ([float('inf')] * len(BoundsCommunityUptakeExport)),\
        [Met + "_c" + str(Idx) for Idx,Mdl in enumerate(TheModels) for Met in Mdl[3]] + [Met + "_u" for Met in ExchangedMets],\
        [Rxn + "_c" + str(Idx) for Idx,Mdl in enumerate(TheModels) for Rxn in Mdl[4]] + [Key + "_u" for Key in BoundsCommunityUptakeExport],\
        [Rxn + "_c" + str(Idx) for Idx,Mdl in enumerate(TheModels) for Rxn in Mdl[5]] + [Key + "_u" for Key in BoundsCommunityUptakeExport]]