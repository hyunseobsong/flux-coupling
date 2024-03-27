def RemoveBlocked(LPMdlParam,CommMdlParam,GrowthRateParam,BiomassIndices):
    # BioRxnNames = (CommMdlParam[4][Num] for Num in BiomassIndices)
    
    TheVars = LPMdlParam.getVars()
    LPMdlParam.addLConstr(sum(TheVars[Num] for Num in BiomassIndices) == GrowthRateParam)#∑X_k = X_0
    
    for Idx,Var in enumerate(TheVars):#LB*X_k ≤ V_k ≤ UB*X_k
        if CommMdlParam[4][Idx][-1] != 'u':
            LPMdlParam.addLConstr((TheVars[BiomassIndices[int(Var.varname[-1])]] * CommMdlParam[1][Idx] / GrowthRateParam) <= Var)
            LPMdlParam.addLConstr(Var <= (TheVars[BiomassIndices[int(Var.varname[-1])]] * CommMdlParam[2][Idx] / GrowthRateParam))
    LPMdlParam.setObjective(0)
    LPMdlParam.optimize()
    
    EachVarVal = [(Var.varname,Var.X) for Var in TheVars]
    DefinitelyUnblocked = [abs(Var[1]) > 1e-7 for Var in EachVarVal]
    
    Blocked = [False] * len(TheVars)
    FluxVectorTable = []
    # raise Exception("Sum", sum(DefinitelyUnblocked), min(abs(Var.X) for Var in TheVars if Var.X))
    for Idx,Var in enumerate(TheVars):
        if DefinitelyUnblocked[Idx]:
            FluxVectorTable.append((EachVarVal[Idx][1],tuple(Var for Var in EachVarVal if Var[1])))
            continue
        # print(Var.varname)
        TheConstr = LPMdlParam.addLConstr(Var >= 1e-7)
        LPMdlParam.optimize()
        LPMdlParam.remove(TheConstr)
        if LPMdlParam.status != 2:
            if CommMdlParam[1][Idx] < 0:
                TheConstr = LPMdlParam.addLConstr(Var <= -1e-7)
                LPMdlParam.optimize()
                LPMdlParam.remove(TheConstr)
                if LPMdlParam.status != 2:
                    Blocked[Idx] = True
                else:
                    FluxVectorTable.append((Var.X,tuple((InnerVar.varname,InnerVar.X) for InnerVar in TheVars if InnerVar.X)))
            else:
                Blocked[Idx] = True
        else:
            FluxVectorTable.append((Var.X,tuple((InnerVar.varname,InnerVar.X) for InnerVar in TheVars if InnerVar.X)))
    LPMdlParam.remove([Var for Idx,Var in enumerate(TheVars) if Blocked[Idx]])
    LPMdlParam.update()
    LPMdlParam.remove([Constr for Constr in LPMdlParam.getConstrs() if not LPMdlParam.getRow(Constr).size()]) #Remove constraints like "0.0 = nonzero value", which happens if the variables in the original constraint were all removed
    LPMdlParam.optimize()
    print("Removing variables",*[Var for Var in LPMdlParam.getVars() if "BIOMASS" in Var.varname],len(LPMdlParam.getVars()),set(Constr.rhs for Constr in LPMdlParam.getConstrs()))
    #Model.getRow(constraint) 
    # quit()
    
    
    TrimmedComm = [None] * len(CommMdlParam)
    TrimmedComm[0] = CommMdlParam[0][:,[not Bool for Bool in Blocked]]
    
    RetainedMets = TrimmedComm[0].any(axis=1)
    TrimmedComm[0] = TrimmedComm[0][RetainedMets,:]
    
    TrimmedComm[3] = [CommMdlParam[3][Idx] for Idx,Bool in enumerate(RetainedMets) if Bool]
    for Num in 1,2,4,5:
        TrimmedComm[Num] = [CommMdlParam[Num][Idx] for Idx,Bool in enumerate(Blocked) if not Bool]
    
    # assert all(Rxn in TrimmedComm[4] for Rxn in DB_Rxns)
    
    # BiomassIndices = [TrimmedComm[4].index(Name) for Name in BioRxnNames]
    # print("New BiomassIndices",BiomassIndices)
    #Removing the variables from the model isn't enough. E.g., you have a constraint (Var_C2 == 1), then you remove Var_C2. You're left with (0.0 == 1).
    # TheMdl = gp.Model()
    
    # #Adding back SteadyCom constraints
    # TheVars = TheMdl.addVars(len(TrimmedComm[1]),lb=(LB if TrimmedComm[4][Idx][-1] == 'u' else -float('inf') for Idx,LB in enumerate(TrimmedComm[1])),ub=(UB if TrimmedComm[4][Idx][-1] == 'u' else float('inf') for Idx,UB in enumerate(TrimmedComm[2])),name=TrimmedComm[4]).values()
    
    # TheMdl.addMConstr(TrimmedComm[0],TheVars,'=',[0]*TrimmedComm[0].shape[0])#Steady state constraint
    # for Idx,Var in enumerate(TheVars):#LB*X_k ≤ V_k ≤ UB*X_k
    #     if TrimmedComm[4][Idx][-1] != 'u':
    #         TheMdl.addLConstr((TheVars[BiomassIndices[int(TrimmedComm[4][Idx][-1])]] * TrimmedComm[1][Idx] / GrowthRateParam) <= Var)
    #         TheMdl.addLConstr(Var <= (TheVars[BiomassIndices[int(TrimmedComm[4][Idx][-1])]] * TrimmedComm[2][Idx] / GrowthRateParam))
    # TheMdl.addLConstr(sum(TheVars[Num] for Num in BiomassIndices) == GrowthRateParam)#∑X_k = X_0
    # TheMdl.write("FromTrimmed.mps")
    # TheMdl.write("FromTrimmed.lp")
    # print("TheMdl.status",TheMdl.status,TheVars[-1].varname)
    
    return LPMdlParam,TrimmedComm,FluxVectorTable