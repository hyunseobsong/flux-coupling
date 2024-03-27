from copy import deepcopy

def TwoComm(MetMdl):
    EachEC = [deepcopy(MetMdl) for _ in range(2)]
    #KO for auxotrophy
    argH = MetMdl[4].index('ARGSL')
    lysA = MetMdl[4].index('DAPDC')
    
    #KO for inability to export
    argO = MetMdl[4].index('ARGt3pp')
    lysO = MetMdl[4].index('LYSt3pp')
    
    EachEC[0][1][argH] = EachEC[0][2][argH] = EachEC[1][1][lysA] = EachEC[1][2][lysA] = EachEC[0][1][argO] = EachEC[0][2][argO] = EachEC[1][1][lysO] = EachEC[1][2][lysO] = 0
    
    ArgEx = MetMdl[4].index("EX_arg__L_e")
    LysEx = MetMdl[4].index("EX_lys__L_e")
    
    #Preventing uptake of metabolites it can synthesize
    EachEC[0][1][LysEx] = EachEC[1][1][ArgEx] = 0
    
    #Allow uptake of metabolites it's auxotrophic in
    EachEC[0][1][ArgEx] = EachEC[1][1][LysEx] = -1
    
    return EachEC

def ThreeComm(MetMdl):
    EachEC = [deepcopy(MetMdl) for _ in range(3)]
    #KO for auxotrophy
    argH = MetMdl[4].index('ARGSL')
    lysA = MetMdl[4].index('DAPDC')
    metA = MetMdl[4].index('HSST')
    
    #KO for inability to export
    argO = MetMdl[4].index('ARGt3pp')
    lysO = MetMdl[4].index('LYSt3pp')
    yjeH = MetMdl[4].index('METt3pp')
    
    EachEC[0][1][argH] = EachEC[0][2][argH] = EachEC[0][1][metA] = EachEC[0][2][metA] = \
    EachEC[1][1][lysA] = EachEC[1][2][lysA] = EachEC[1][1][metA] = EachEC[1][2][metA] = \
    EachEC[2][1][argH] = EachEC[2][2][argH] = EachEC[2][1][lysA] = EachEC[2][2][lysA] = \
    EachEC[0][1][argO] = EachEC[0][2][argO] = EachEC[0][1][yjeH] = EachEC[0][2][yjeH] = \
    EachEC[1][1][lysO] = EachEC[1][2][lysO] = EachEC[1][1][yjeH] = EachEC[1][2][yjeH] = \
    EachEC[2][1][argO] = EachEC[2][2][argO] = EachEC[1][1][lysO] = EachEC[1][2][lysO] = 0
    
    ArgEx = MetMdl[4].index("EX_arg__L_e")
    LysEx = MetMdl[4].index("EX_lys__L_e")
    MetEx = MetMdl[4].index("EX_met__L_e")
    
    #Preventing uptake of metabolites it can synthesize
    EachEC[0][1][LysEx] = EachEC[1][1][ArgEx] = EachEC[2][1][MetEx] = 0
    
    #Allow uptake of metabolites it's auxotrophic in
    EachEC[0][1][ArgEx] = EachEC[0][1][MetEx] = EachEC[1][1][LysEx] = EachEC[1][1][MetEx] = EachEC[2][1][ArgEx] = EachEC[2][1][LysEx] = -1
    return EachEC
def FourComm(MetMdl):
    FourECs = [deepcopy(MetMdl) for _ in range(4)]
    argH = MetMdl[4].index('ARGSL')
    lysA = MetMdl[4].index('DAPDC')
    metA = MetMdl[4].index('HSST')
    ilvE = MetMdl[4].index('PPNDH')

    argO = MetMdl[4].index('ARGt3pp')
    lysO = MetMdl[4].index('LYSt3pp')
    yjeH = MetMdl[4].index('METt3pp')
    yddG = MetMdl[4].index('PHEt2rpp')

    #"all relevant AA exchange reactions"
    ArgEx = MetMdl[4].index("EX_arg__L_e")
    MetEx = MetMdl[4].index("EX_met__L_e")
    LysEx = MetMdl[4].index("EX_lys__L_e")
    PheEx = MetMdl[4].index("EX_phe__L_e")

    #Blocking various reactions
    FourECs[0][1][lysA] = FourECs[0][2][lysA] = FourECs[0][1][metA] = FourECs[0][2][metA] = FourECs[0][1][yddG] = FourECs[0][2][yddG] =\
    FourECs[1][1][argH] = FourECs[1][2][argH] = FourECs[1][1][yjeH] = FourECs[1][2][yjeH] = FourECs[1][1][ilvE] = FourECs[1][2][ilvE] =\
    FourECs[2][1][argH] = FourECs[2][2][argH] = FourECs[2][1][lysO] = FourECs[2][2][lysO] = FourECs[2][1][ilvE] = FourECs[2][2][ilvE] =\
    FourECs[3][1][argO] = FourECs[3][2][argO] = FourECs[3][1][lysA] = FourECs[3][2][lysA] = FourECs[3][1][metA] = FourECs[3][2][metA] = 0

    #organism-specific AA uptake rates
    FourECs[0][1][ArgEx] = FourECs[0][2][PheEx] = FourECs[1][1][LysEx] = FourECs[1][2][MetEx] = FourECs[2][2][LysEx] = FourECs[2][1][MetEx] = FourECs[3][2][ArgEx] = FourECs[3][1][PheEx] = 0

    FourECs[0][1][LysEx] = FourECs[0][1][MetEx] = FourECs[1][1][ArgEx] = FourECs[1][1][PheEx] = FourECs[2][1][ArgEx] = FourECs[2][1][PheEx] = FourECs[3][1][LysEx] = FourECs[3][1][MetEx] = -1#maximum uptake rate for cross feeding AAs
    return FourECs