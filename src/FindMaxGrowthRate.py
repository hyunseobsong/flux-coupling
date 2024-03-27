import gurobipy as gp
def FindMaxGrowthRate(TheMetMdl,BioString,µ = .5): #µ:"Initial guess". Constraints will need to be removed and re-added once µ changes
    TheMdl = gp.Model()
    TheVars = TheMdl.addVars(len(TheMetMdl[1]),lb=(LB if TheMetMdl[4][Idx][-1] == 'u' else -float('inf') for Idx,LB in enumerate(TheMetMdl[1])),ub=(UB if TheMetMdl[4][Idx][-1] == 'u' else float('inf') for Idx,UB in enumerate(TheMetMdl[2])),name=TheMetMdl[4]).values()
    TheMdl.addMConstr(TheMetMdl[0],TheVars,'=',[0] * TheMetMdl[0].shape[0])#Steady state constraint

    BiomassIndices = [Idx for Idx in range(len(TheMetMdl[4])) if BioString in TheMetMdl[4][Idx]]
    # print(*BiomassIndices)
    TotalCommunityBiomass = 1

    RangeOfµ = [None] * 2
    while True:
        #LB*X_k ≤ V_k ≤ UB*X_k
        ConstraintsForEqNine =\
        [TheMdl.addLConstr((TheVars[BiomassIndices[int(TheMetMdl[4][Idx][-1])]] * TheMetMdl[1][Idx] / µ) <= TheVars[Idx]) for Idx in range(len(TheVars)) if TheMetMdl[4][Idx][-1] != 'u'] +\
        [TheMdl.addLConstr(TheVars[Idx] <= (TheVars[BiomassIndices[int(TheMetMdl[4][Idx][-1])]] * TheMetMdl[2][Idx] / µ)) for Idx in range(len(TheVars)) if TheMetMdl[4][Idx][-1] != 'u']
        
        TheMdl.setObjective(sum(TheVars[Num] for Num in BiomassIndices) / µ, gp.GRB.MAXIMIZE)
        TheMdl.optimize()
        TheMdl.remove(ConstraintsForEqNine)
        # print("Glucose uptake",[(Var.lb,Var.ub) for Var in TheVars if 'EX_glc__D_e' in Var.varname])
        # print(µ,TheMdl.ObjVal,RangeOfµ,[TheVars[Num].X for Num in BiomassIndices],[Var.X for Var in TheVars if 'EX_glc__D_e' in Var.varname])
        # print("Highest flux",max(abs(Var.X) for Var in TheVars))
        if abs(TheMdl.ObjVal - TotalCommunityBiomass) < 1e-8:
            print("Final µ",µ,TheMdl.ObjVal,[TheVars[Num].X / µ for Num in BiomassIndices])
            TheMdl.update()
            return µ,TheMdl,BiomassIndices
        elif RangeOfµ[0] and RangeOfµ[1]:
            # Oldµ = µ
            if TheMdl.ObjVal >= TotalCommunityBiomass:
                RangeOfµ[0] = µ
            else:
                RangeOfµ[1] = µ
            µ = (RangeOfµ[0] + RangeOfµ[1]) / 2
        elif TheMdl.ObjVal >= TotalCommunityBiomass and not RangeOfµ[0]:
            RangeOfµ[0] = µ
            µ*=1.5
        elif not RangeOfµ[1]:
            RangeOfµ[1] = µ
            µ*=0.5