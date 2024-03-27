import numpy as np
def FluxCouplingAnalysis(TheMdl,FluxVectorTable):
    print("FluxCouplingAnalysis")
    TheVars = TheMdl.getVars()
    RxnIndicesDict = dict((Var.varname,Idx) for Idx,Var in enumerate(TheVars))
    
    FullyCoupled = [None] * len(TheVars)#E.g. if FCTbl[2] = 0, then rxns 0 and 2 are fully coupled
    FCTbl = np.identity(len(TheVars),np.uint8)
    
    #len(TheVars) should be equal to len(FluxVectorTable)
    TheMdl.setObjective(0)
    for Idx,Var in enumerate(TheVars): #RxnA
        if FullyCoupled[Idx] != None: #Earlier it was found that the Idx'th reaction was coupled to the n'th reaction, where case n < Idx. No need to find Idx'th reaction's couplings since they match n'th.
            continue
        #Each element of FluxVectorTable is (Idx'th nonzero flux,tuple of other element's (Rxn name, Rxn nonzero flux))
        for Elem in FluxVectorTable[Idx][1]: #RxnB
            if Var.varname == Elem[0]:
                continue
            TheConstr = TheMdl.addLConstr(Var == FluxVectorTable[Idx][0])
            OtherConstr = TheMdl.addLConstr(TheVars[RxnIndicesDict[Elem[0]]] == 0)
            TheMdl.optimize()
            TheMdl.remove((TheConstr,OtherConstr))
            
            AtoB = TheMdl.status != 2#If AtoB is true, for RxnA to be active, RxnB must be active
            
            TheConstr = TheMdl.addLConstr(TheVars[RxnIndicesDict[Elem[0]]] == Elem[1])
            OtherConstr = TheMdl.addLConstr(Var == 0)
            TheMdl.optimize()
            TheMdl.remove((TheConstr,OtherConstr))
            
            BToA = TheMdl.status != 2#If BToA is true, for RxnB to be active, RxnA must be active
            
            if AtoB:
                if BToA:
                    TheConstr = TheMdl.addLConstr(Var == FluxVectorTable[Idx][0])
                    OtherConstr = TheMdl.addLConstr(TheVars[RxnIndicesDict[Elem[0]]] <= Elem[1] - 1e-6) #After fixing RxnA's flux, can RxnB's flux be less than what was observed previously?
                    TheMdl.optimize()
                    TheMdl.remove(OtherConstr)
                    if TheMdl.status != 2:
                        OtherConstr = TheMdl.addLConstr(TheVars[RxnIndicesDict[Elem[0]]] >= Elem[1] + 1e-6) #After fixing RxnA's flux, can RxnB's flux be more than what was observed previously?
                        TheMdl.optimize()
                        TheMdl.remove(OtherConstr)
                        if TheMdl.status != 2:
                            FCTbl[Idx][RxnIndicesDict[Elem[0]]] = FCTbl[RxnIndicesDict[Elem[0]]][Idx] = 1
                            FullyCoupled[RxnIndicesDict[Elem[0]]] = Idx
                        else:
                            FCTbl[Idx][RxnIndicesDict[Elem[0]]] = FCTbl[RxnIndicesDict[Elem[0]]][Idx] = 2
                    else:
                        FCTbl[Idx][RxnIndicesDict[Elem[0]]] = FCTbl[RxnIndicesDict[Elem[0]]][Idx] = 2
                    TheMdl.remove(TheConstr)
                else:
                    FCTbl[Idx][RxnIndicesDict[Elem[0]]] = 3
                    FCTbl[RxnIndicesDict[Elem[0]]][Idx] = 4
            elif BToA:
                FCTbl[Idx][RxnIndicesDict[Elem[0]]] = 4
                FCTbl[RxnIndicesDict[Elem[0]]][Idx] = 3
        # raise Exception(FCTbl[0,:].tolist())
    for Idx,Num in enumerate(FullyCoupled): #Fully coupled reactions are coupled to the same reactions
        if Num != None:
            FCTbl[Idx,:] = FCTbl[Num,:]
            # FCTbl[:,Idx] = FCTbl[:,Num]
    return FCTbl