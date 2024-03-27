import xmltodict
import numpy as np

def parseXML(filename):
    with open(filename) as InputFile:
        Parsed = xmltodict.parse(InputFile.read())
    # Getting the values used to represent bounds    
    MyDict = dict((Elem['@id'],float(Elem['@value'])) for Elem in Parsed['sbml']['model']['listOfParameters']['parameter'])
    EveryRxn = Parsed['sbml']['model']['listOfReactions']['reaction']
    # print(len(EveryRxn)) #With e_coli_core.xml, 95 reactions as expected

    # print("Number of metabolites: ",len(Parsed['sbml']['model']['listOfSpecies']['species']))
    MetIDs = [child['@id'].lstrip('M_')
        for child in Parsed['sbml']['model']['listOfSpecies']['species']]
    TheMatrix = np.zeros((len(MetIDs), len(EveryRxn)))
    TheLBs = np.zeros(len(EveryRxn))
    TheUBs = np.zeros(len(EveryRxn))
    RxnIDs = [None] * len(EveryRxn)
    RxnNames = [None] * len(EveryRxn)

    for Idx,TheRxn in enumerate(EveryRxn):
        RxnIDs[Idx] = TheRxn['@id'][2:]
        try:
            RxnNames[Idx] = TheRxn['@name']
        except KeyError:
            RxnNames[Idx] = "NoName"
        TheLBs[Idx] = MyDict[TheRxn['@fbc:lowerFluxBound']]
        TheUBs[Idx] = MyDict[TheRxn['@fbc:upperFluxBound']]
        
        if 'listOfReactants' in TheRxn:
            TheReactants = TheRxn['listOfReactants']['speciesReference']
            if isinstance(TheReactants, dict):
                TheReactants = [TheReactants]
            for Metabolite in TheReactants:
                TheMatrix[MetIDs.index(Metabolite['@species'].lstrip('M_'))][Idx] = - \
                                       float(Metabolite['@stoichiometry'])
        if 'listOfProducts' in TheRxn:
            TheProducts = TheRxn['listOfProducts']['speciesReference']
            if isinstance(TheProducts, dict):
                TheProducts = [TheProducts]
            for Species in TheProducts:
                TheMatrix[MetIDs.index(Species['@species'].lstrip('M_'))
                                       ][Idx] = float(Species['@stoichiometry'])
    return [TheMatrix, TheLBs, TheUBs, MetIDs, RxnIDs, RxnNames]