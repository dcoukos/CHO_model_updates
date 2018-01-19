'''
    This file contains functions pertaining to the treatment of the data that 
    is returned from BRENDA, including functions that find the KEGG Ids of that
    data, as well as functions that organize that data and simplify it in order
    to keep only the necessary information for improving the model. 
'''

#TODO: WRITE TESTS FOR THIS CODE

import os
import json
import requests
import copy
import pickle
import cobra
from fuzzywuzzy import fuzz
from enum import Enum, auto  #Enum breaks spyder autocomplete

'''
def auto():
    ''TODO: delete class. Is here only so auto-complete will work in spyder while
    working on project with ENUM class. '   
    return None
'''
def initializeModelUpdate():
    '''Creates deep copy of BiGG Model representation, to which updated 
        information will be added. '''
    global MODEL_UPDATE
    global BIGG_MODEL
    MODEL_UPDATE = copy.deepcopy(BIGG_MODEL)

def treatTurnoverData():
    '''Filters and eliminates unnecessary data, choosing best turnover data.
    
    This function is meant specifically for data of the DataType.turnover Enum. 
    All of the selection processes (filtering by organism, direction, 
    wild-type, and magnitude) occur in functions called by 
    treatTurnoverData().
    '''
    
    global MODEL_UPDATE
    potential_updates = {} 
    
    storeBiggRepresentation()
    treated_output = openJson('JSONs/treated_BRENDA_output.json')
    if os.stat('JSONs/BRENDA_KEGG_IDs.json') == 0:
        brendaToKegg(treated_output) 
    brenda_keggs = correctJson('JSONs/BRENDA_KEGG_IDs.json')
    brenda_no_kegg = openJson('JSONs/BRENDA_no_KEGG.json')
    
    
    eliminateUnmatchable(potential_updates,treated_output, brenda_keggs, 
                         brenda_no_kegg)
    selectBestData(potential_updates, DataType.turnover)
    applyNewData(potential_updates, DataType.turnover)
    
    
def eliminateUnmatchable(potential_updates,treated_brenda_output, 
                         brenda_keggs, brenda_no_kegg):
    '''Eliminates entries from BRENDA which cannot be matched to model data. 
    
    The return from BRENDA is quite large, and only useful if pertinent to 
    metabolites in the model. The data must therefore be pruned. When correct
    matches are found they are added to a new data structure, which is also 
    organized for forward and backwards reactions. 
    
    Args:
        potential_updates: empty dict which will contain Enzyme entries. This
            dict is where the incoming data will be processed and selected. 
        treated_brenda_output: This is a dict containing data from BRENDA. Keys are 
            reaction Ids. First order of keys are BiGG reaction Ids, and the second order 
            keys are the names of metabolites in BRENDA. Thereafter there is a 
            list of outputs, obtained programatically from BRENDA.
        brenda_keggs: dict containing KEGG ids and metabolite names in BRENDA. 
            First order of keys are BiGG reaction Ids, and the second order 
            keys are the names of metabolites in BRENDA.
        brenda_no_kegg: dict containing metabolites for which no KEGG id was 
            found. 
    
    Note:
        Data passed to this function as treated_brenda_output, brenda_keggs, 
        and brenda_no_kegg must at least have this format:
            { reaction1:{ metabolite1: ..., 
                          metabolite2: ...,
                          metabolite3: ...,
                         },
                ...     
            }
    
    '''
    #TODO: determine whether matching by name is even worth it. 

    unmatched = matchById(potential_updates, brenda_keggs, 
                          treated_brenda_output, DataType.turnover)
    matchByName(potential_updates, unmatched, treated_brenda_output, 
                DataType.turnover)
    matchByName(potential_updates, brenda_no_kegg, treated_brenda_output, 
                DataType.turnover)

    return potential_updates    

def matchById(potential_updates, brenda_keggs, treated_brenda_output, 
              data_type):
    '''Tries to match metabolites by KEGG Id. 
    
    Args:
        potential_updates: empty dict which will contain Enzyme entries. This
            dict is where the incoming data will be processed and selected. Is 
            updated by this function. 
        brenda_keggs: dict of KEGG codes (keys) and corresponding names 
            (values)
        treated_brenda_output: This is a dict containing data from BRENDA. 
            Keys are reaction Ids. First order of keys are BiGG reaction Ids, 
            and the second order keys are the names of metabolites in BRENDA. 
            Thereafter there is a list of outputs, obtained programatically 
            from BRENDA.
        data_type: DataType which tells matchById how to add new data to the 
            potential_updates dict.
    
    Return:  
        unmatched: dict of BRENDA metabolites that could not be matched by name. 
    '''
    unmatched = {}
    global BIGG_MODEL

    for reaction in brenda_keggs:  #Checking for matches for each reaction.
        unmatched[reaction] = []
        if reaction not in potential_updates:
            potential_updates[reaction] = Enzyme(reaction, getEcNumber(
                                                treated_brenda_output[
                                                        reaction]))             
        for ID, name in brenda_keggs[reaction].items():
            try: 
                bigg_name = BIGG_MODEL[reaction].with_kegg(ID) # = checked.
                if bigg_name in BIGG_MODEL[reaction].forward: 
                    data = getData(treated_brenda_output[reaction], name, 
                               'forward', data_type)                    
                    if bigg_name not in potential_updates[reaction].forward:
                        potential_updates[reaction].forward[bigg_name] = []
                    for entry in data:
                        potential_updates[reaction].forward[bigg_name].append(
                            MetaboliteCandidate(bigg_name, kegg=ID,
                                                data=entry))
                elif bigg_name in BIGG_MODEL[reaction].backward:
                    data = getData(treated_brenda_output[reaction], name, 
                               'forward', data_type)
                    if bigg_name not in potential_updates[reaction].backward:
                        potential_updates[reaction].backward[bigg_name] = []
                    for entry in data:
                        potential_updates[reaction].backward[bigg_name].append(
                            MetaboliteCandidate(bigg_name, kegg=ID, 
                                                data=entry))
                else:
                    raise DataMissingError()
                
            except KeyError:
                unmatched[reaction].append(name)
            except DataMissingError:
                print('Metabolite not added correctly to BIGG_MODEL.')              
    return unmatched
  
    
def matchByName(potential_updates, unmatched, treated_brenda_output, 
                data_type):
    '''Tries to fuzzy match metabolite names that cannot be matched via KEGG.
        
    Args:
        potential_updates: dict which will contain Enzyme entries. This
            dict is where the incoming data will be processed and selected. 
            Metabolite entries which are matched are added here.
        unmatched: dict of the metabolite entries from BRENDA which have not 
            yet been matched to BiGG metabolite by previous function calls.
        treated_brenda_output: This is a dict containing data from BRENDA. 
            Keys are reaction Ids. First order of keys are BiGG reaction Ids, 
            and the second order keys are the names of metabolites in BRENDA. 
            Thereafter there is a list of outputs, obtained programatically 
            from BRENDA.
        data_type: DataType which tells matchById how to add new data to the 
            potential_updates dict.
    
    '''

    global BIGG_MODEL
    
    for reaction in unmatched:
        for metabolite in unmatched[reaction]:
            #REACTANT IS A KEY TO MATCH WITH METABOLITE
            for reactant in BIGG_MODEL[reaction].forward:
                if fuzz.token_set_ratio(metabolite, reactant.name) > 0.85:
                    data = getData(treated_brenda_output[reaction], metabolite, 
                               'forward', data_type)
                    for entry in data:
                        potential_updates[reaction].forward[
                            reactant.name].append(MetaboliteCandidate(
                            reactant.name, data=entry))
            for product in BIGG_MODEL[reaction].backward:
                if fuzz.token_set_ratio(metabolite, product.name) > 0.85:
                    data = getData(treated_brenda_output[reaction], metabolite, 
                               'forward', data_type)
                    for entry in data:
                        potential_updates[reaction].backward[
                            product.name].append(MetaboliteCandidate(
                            product.name, data=entry))

def selectBestData(model_updates, data_type):
    '''Selects best data to add to model based on selection criteria
    
    This function selects the best data based on Organism, wild-type, and the 
    magnitude of the data. 
    
    Args:
        model_updates: contains new data for the model, is filtered by this 
            function. The final candidates are returned through this arg.
        data_type: instance of DataType Enum, describing where new data is to 
            be placed.
            
    '''
    #first the data must be filtered by organism. 
    
    filtered_by_organism = {}    
    for reaction in model_updates:
        filtered_by_organism[reaction] = Enzyme.copyOnlyEnzyme(model_updates
                            [reaction])
        for metabolite_name in model_updates[reaction].forward:
            filtered_by_organism[reaction].forward[metabolite_name] = []
            closest_organism, indices = findClosestOrganism(model_updates[
                    reaction].forward[metabolite_name])
            best_found = False
            for organism in Organism:
                while not best_found:
                    for index in indices[organism]:
                        filtered_by_organism[reaction].forward[
                                metabolite_name].append(indices[
                                closest_organism][index])
                        if organism is closest_organism:
                            best_found = True   
        
        #SAME as above but for reverse reaction                    
        for metabolite_name in model_updates[reaction].backward:
            filtered_by_organism[reaction].backward[metabolite_name] = []
            closest_organism, indices = findClosestOrganism(model_updates[
                    reaction].backward[metabolite_name])
            best_found = False
            for organism in Organism:
                while not best_found:
                    for index in indices[organism]:
                        filtered_by_organism[reaction].backward[
                                metabolite_name].append(indices[
                                closest_organism][index])
                        if organism is closest_organism:
                            best_found = True    
    
    filtered_by_wild_type = {}
    for reaction in filtered_by_organism:
        filtered_by_wild_type[reaction] = Enzyme.copyOnlyEnzyme(
                filtered_by_organism[reaction])
        for metabolite_name in filtered_by_organism[reaction].forward:
            wild_type = []
            for entry in model_updates[reaction].forward[metabolite_name]:
                if entry.wild_type is not None:
                    wild_type.append(entry)
            if wild_type != []:
                filtered_by_wild_type[reaction].forward[
                        metabolite_name] = wild_type
            else:
                filtered_by_wild_type[reaction].forward[
                        metabolite_name] = copy.copy(filtered_by_organism
                        [reaction].forward[metabolite_name])
            
        #SAME as above but for reverse reactions
        for metabolite_name in filtered_by_organism[reaction].backward:
            wild_type = []
            for entry in model_updates[reaction].backward[metabolite_name]:
                if entry.wild_type is not None:
                    wild_type.append(entry)
            if wild_type != []:
                filtered_by_wild_type[reaction].backward[metabolite_name] = \
                wild_type
            else:
                filtered_by_wild_type[reaction].backward[
                        metabolite_name] = copy.copy(filtered_by_organism
                        [reaction].backward[metabolite_name])
    
    #-----Choose data according to magnitude preference.-----------------
    #work on filtered_by_wild_type)
    
    if data_type is DataType.turnover:
        for reaction in filtered_by_wild_type:
            for metabolite in filtered_by_wild_type[reaction].forward:
                model_updates[reaction].forward = chooseHighestTurnover(
                        filtered_by_wild_type[reaction].forward[metabolite])
            for metabolite in filtered_by_wild_type[reaction].backward:
                model_updates[reaction].backward = chooseHighestTurnover(
                        filtered_by_wild_type[reaction].backward[metabolite])
  
def chooseHighestTurnover(metabolite_entries):
    '''Selects entry with highest turnover
    
    Args:
        metabolite_entries: list of selected entries for a metabolite.
    
    Returns:
        entry with the highest turnover
    '''
    highest_index = 0
    highest_turnover = 0
    for index, entry in enumerate(metabolite_entries):
        if entry.turnover > highest_turnover:
            highest_index = index
    
    return metabolite_entries[highest_index]          
        
def findClosestOrganism(metabolite_entries):
    '''finds which organism closest to hamster is in list of MetaboliteCandidate
    
    Args:
        metabolite_entries: list of MetaboliteCandidate
        
    Returns: 
        Organism: enum of closest organism closest to hamster
        dict: dictionary of Organism and indices matching to the organism
    '''
    has_hamster = False
    has_mouse = False
    has_rat = False
    has_human = False
    has_fly = False
    has_yeast = False
    has_coli = False
    indices = {
        Organism.hamster : [],
        Organism.mouse: [],
        Organism.rat: [],
        Organism.human: [],
        Organism.fly: [],
        Organism.yeast: [],
        Organism.coli: []}

    for index, metabolite in enumerate(metabolite_entries):
        if metabolite.organism == 'Cricetulus griseus':
            indices[Organism.hamster].append(index)
            has_hamster = True
        elif metabolite.organism == 'Mus musculus':
            has_mouse = True
            indices[Organism.mouse].append(index)
        elif metabolite.organism == 'Rattus norvegicus':
            has_rat = True
            indices[Organism.rat].append(index)
        elif metabolite.organism == 'Homo sapiens':
            has_human = True
            indices[Organism.human].append(index)
        elif metabolite.organism == 'Drosophila melanogaster':
            has_fly = True
            indices[Organism.fly].append(index)
        elif metabolite.organism == 'Saccharomyces cerevisiae':
            has_yeast = True
            indices[Organism.yeast].append(index)
        elif metabolite.organism == 'Escherichia coli':
            has_coli = True
            indices[Organism.coli].append(index)
    
    if has_hamster:
        return Organism.hamster, indices
    elif has_mouse:
        return Organism.mouse, indices
    elif has_rat:
        return Organism.rat, indices
    elif has_human:
        return Organism.human, indices
    elif has_fly:
        return Organism.fly, indices
    elif has_yeast:
        return Organism.yeast, indices
    elif has_coli:
        return Organism.coli, indices
    else:
        return None

def brendaToKegg(data): 
    '''Tries to match BRENDA metabolite names to KEGG ids. 
        
        BRENDA metabolites that can be matched to a KEGG id are stored in 
        JSONs/BRENDA_KEGG_IDs.json and those that cannot are stored in 
        JSONs/BRENDA_no_KEGG.json. This is to counteract the long time it takes 
        to collect data. 
    
        
    Args: 
        data: This is a dict containing data from BRENDA. Keys are reaction Ids. 
    
    Data passed to this function must at least have this format:
        { reaction1:{ metabolite1: ..., 
                      metabolite2: ...,
                      metabolite3: ...,
                     },
            ...     
        }
    
    Note: I have already run this from Cossio's computer, and the JSON files
            are already populated. 
            
     TODO: ARE NOT POPULATED CORRECTLY
        
    '''
    
    
    for reaction in data:
        [brenda_kegg_ids, brenda_no_kegg] = getKeggIds(data[reaction].keys())
        kegg_data = {}
        no_kegg_data = []
        try:
            kegg_data_r = open('JSONs/BRENDA_KEGG_IDs.json')
        except ValueError:
            kegg_data_r = {}
        try:
            no_kegg_data_r = open('JSONs/BRENDA_no_KEGG.json')
        except ValueError:
            no_kegg_data_r = {}
        for metabolite in brenda_kegg_ids:
            kegg_data.update({metabolite:brenda_kegg_ids[metabolite]})
        
        for metabolite in brenda_no_kegg:
            no_kegg_data.append(metabolite)
    
        kegg_data_r.update({reaction:kegg_data})
        no_kegg_data_r.update({reaction:no_kegg_data})
        write('JSONs/BRENDA_KEGG_IDs.json', kegg_data_r)
        write('JSONs/BRENDA_no_KEGG.json', no_kegg_data_r)

def correctJson(path):
    ''' Ensures that KEGG codes are the keys and metabolite names are the values. 
    
    When running brendaToKegg() in Jupyter Notebook I noticed a funny error:
    some, but not all, of the key:value entries were inverted. This appeared to
    occur 'at random'. This function ensures that the KEGG code is the value.
    
    Args:
        path: filepath to KEGG ids and BRENDA metabolite names. 
        
    Return:
        corrected_json: JSON object with KEGG codes as the dict keys.   
    '''
    
    brenda_keggs = openJson(path)
    corrected_json = {}
    
    for reaction in brenda_keggs.keys():
        for supposed_code in brenda_keggs[reaction].keys():
            if is_number(supposed_code[2:]):
                corrected_json[reaction][supposed_code] = brenda_keggs \
                                                            [reaction] \
                                                            [supposed_code]
            else:
                corrected_json[reaction][brenda_keggs[reaction]
                                            [supposed_code]] = supposed_code
    
    return corrected_json
            
def openJson(path):
    '''Shortened call to open JSON files.'''
    with open(path) as r: 
        return json.load(r)

def write(path, data): 
    '''Shortened call to JSON dumps with indent = 4'''
    with open(path, 'w') as wr:
            json.dump(data, wr, indent = 4)
            
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass

    return False


def getData(reaction, metabolite, directionality, data_type):
    '''Chooses data to give to matchers. 
    
    This function is intended to improve code reusability so that the parent 
    functions can treat the data without knowing what type of data it is. 
    
    Args:
        reaction: treated brenda output for one EC number. This is a dict
            where the keys are metabolites, and the values are lists of entries
            pertaining to that metabolite name. 
        metabolite: metabolite name
        directionality: string indicating whether .forward or .backward 
            attribute is being assessed. 
        data_type: instance of DataType Enum. Tells the function which data is 
            pertinent. 
        
    Returns:
        metabolite_data: a list containing data which is used to construct
            MetaboliteCandidate instances. 
    '''
    
    metabolite_data = []    
    
    for entry in getattr(reaction, directionality):
        metabolite_entry  = {}
        if data_type is DataType.turnover:
            metabolite_entry['turnover'] =  reaction[metabolite]['turnoverNumber']
        elif data_type is DataType.specific_activity:
            metabolite_entry['specific activity'] = reaction[metabolite] \
                                                       ['specificActivity']
        else:
            metabolite_entry['molecular weight'] = reaction[metabolite] \
                                                       ['molecularWeight']
        try:
            metabolite_entry['wild-type'] = reaction[metabolite]['wild-type']
        except:
            pass
        metabolite_entry['organism'] = reaction[metabolite]['organism']
        metabolite_data.append(metabolite_entry)
    return metabolite_data

def applyNewData(updates, data_type):
    '''Adds updated data to global MODEL_UPDATES variable
    
    Args: 
        updates: updated model information. dict of Enzymes.
        data_type: type of data being added to model.   
    '''
    global MODEL_UPDATE
    for reaction in updates:
        for metabolite in updates[reaction].forward:
            try:
                MODEL_UPDATE[reaction].forward[metabolite] = 

def getEcNumber(reaction):
    '''Gets EC number to make new Enzyme object. "Flat is better than nested"
    '''    
    for metabolite in reaction:
        for entry in reaction[metabolite]:
            try: 
                return entry['ecNumber']
            except KeyError:
                continue
            except:
                raise
    
def cleanMetaboliteNames(brenda_metabolites):
    '''
        TODO: fill out function
    '''
    pass
    
def getKeggIds(reactants):
    '''
    
    TODO: fill out documentation. 
    '''

    reactant_to_KEGG = {}
    reactant_no_KEGG = []
    reactant_counter = 0

    for reactant in reactants[reactant_counter:]:
        try:
            if " - reduced" in reactant:
                reactant = "reduced " + reactant[:-10]
            cts_output = requests.get("http://cts.fiehnlab.ucdavis.edu/"
                                      "service/convert/Chemical%20Name/KEGG/"
                                      + reactant)
            reactant_to_KEGG[str(json.loads(cts_output.text)[0]
                                 ['result'][0])] = reactant
            reactant_counter = reactant_counter + 1
            print('Reactant ' + reactant + ' added to KEGG')
        except:
            # Sometimes the request glitches, so we try again.
            try:
                cts_output = requests.get("http://cts.fiehnlab.ucdavis.edu/"
                                          "service/convert/Chemical%20Name/"
                                          "KEGG/"+reactant)
                if json.loads(cts_output.text)[0]['result']:
                    reactant_to_KEGG[reactant] = str(json.loads(cts_output.
                                                     text)[0]['result'][0])

                else:
                    reactant_no_KEGG.append(reactant)
                    print('Reactant ' + reactant + ' not added to KEGG')
                reactant_counter = reactant_counter + 1
            except:
                continue
    return [reactant_to_KEGG, reactant_no_KEGG]


def storeBiggRepresentation():
    '''Creates and stores a representation of the BiGG Model. 
    
        This creates and pickles a representation of the BiGG Model which 
        should be a little bit easier to work with. 
    '''
    global BIGG_MODEL
    bigg_model = cobra.io.load_json_model('JSONs/BIGG_master_modle.json')
     
    for reaction in bigg_model.reactions:
        BIGG_MODEL[reaction.id] = Enzyme(reaction.id)
        
        for reactant in reaction.reactants: 
            BIGG_MODEL[reaction.id].forward[reactant.name] = Metabolite(
                                                            reactant.name, 
                                                            bigg = reactant.id)
        for reactant in reaction.products: 
            BIGG_MODEL[reaction.id].backward[reactant.name] = Metabolite(
                                                            reactant.name, 
                                                            bigg = reactant.id)
    pickle.dump(BIGG_MODEL, 'bigg_model_object.pickle', protocol = -1)
    
def loadGlobalBigg():
    global BIGG_MODEL
    BIGG_MODEL = pickle.load('bigg_model_object.pickle')

def downloadModelKeggs():
    '''Downloads kegg IDs of metabolites in BiGG model.

    This function is responsible for calling addKeggToMetabolites, which will 
    try to query the chemical translation service for the kegg ids. It handles 
    forward and reverse reactions separately, and if a kegg is found, it 
    updates the with_kegg attribute of the enzyme. 
    
    '''
    global BIGG_MODEL
    
    
    for reaction in BIGG_MODEL:
        for reactant in BIGG_MODEL[reaction].forward:
            kegg_found, kegg_id = addKeggToMetabolites(reactant)[0] 
            if kegg_found:    
                BIGG_MODEL[reaction].with_kegg[kegg_id] = reactant.name
        for product in BIGG_MODEL[reaction].backward:
            kegg_found, kegg_id = addKeggToMetabolites(product)[0] 
            if kegg_found:
                BIGG_MODEL[reaction].with_kegg[kegg_id] = product.name
    
    saveBiggModel()

def addKeggToMetabolites(metabolite):
    '''Adds the kegg id to the local representation of the BiGG Model. 
    
    Uses requests to query the chemical translation service for the kegg id of 
    the metabolite, which is then added to the kegg attribute of the Metabolite.
    
    Returns: 
        bool: True if metabolite kegg was found, false if it was not found. 
            This allows the calling function to add the metabolite to the dict
            of metabolites with a corresponding KEGG ID. 
        kegg: This is the kegg code of the 
    '''
    request_counter = 0
    
    while request_counter < 3:
        try:
            metabolite_name = metabolite.name
            if " - reduced" in metabolite_name:
                metabolite_name = "reduced " + metabolite_name[:-10]
            cts_output = requests.get("http://cts.fiehnlab.ucdavis.edu/"
                                      "service/convert/Chemical%20Name/"
                                      "KEGG/"+metabolite_name)
            kegg = [str(json.loads(cts_output.text)[0]['result']
                                    [0])]
            metabolite.kegg = kegg
            print('Metabolite ' + metabolite_name + 'kegg found.')
            return True, kegg
        except:
            print('EXCEPTED')
            request_counter = request_counter + 1
            continue
    
    return False, None
                
def saveBiggModel():
    '''Stores model using pickle module.
    
    '''
    global BIGG_MODEL
    pickle_out = open('BiGG_model.pickle', 'wb')
    pickle.dump(BIGG_MODEL, pickle_out)
    pickle_out.close()
    
    
BIGG_MODEL = {} #Probably not the best way to do this...
MODEL_UPDATE = {}
    

class DataType(Enum):
    turnover = auto()
    specific_activity = auto()
    molecular_weight = auto()
  
class Organism(Enum):
    hamster = auto()
    mouse = auto()
    rat = auto()
    human = auto()
    fly = auto()
    yeast = auto()
    coli = auto()
    
class Enzyme():
    '''Class containing reactions and metabolites pertaining to an enzyme. 
    
    Attributes:
        ID: the BiGG idenitifier of the enzyme
        EC: Enzyme commision number
        metabolites: a dict of metabolites
        forward: a dict of metabolites participating in the forward reaction
        backward: a dict of metabolites participating in the backward reaction
        with_kegg: a dict of metabolites with a KEGG identifier, where the KEGG
            is the key, and the metabolite is the valye. 
    '''
    def __init__(self, bigg, ec = None):
        self.bigg = bigg
        if ec is not None:
            self.EC = ec
        self.forward = {}
        self.backward = {}
        self.with_kegg = {}
        self.forward_turnover = None
        self.backward_turnover = None
        self.has_f_wild_type = None
        self.has_b_wild_type = None
    
    def returnWithKegg():
        '''Returns a dict of metabolites by their KEGG ids.'''
        
    def updateKegg():
        '''iterates through all metabolites. If they have a Kegg ID, they are 
        added to '''
        pass
    
    def copyOnlyEnzyme(enzyme):
        enz = Enzyme(enzyme.bigg, enzyme.EC)
        return enz
    
    def chooseBiggestTurnovers():
        '''Chooses best value from forward metabolites.
        
            This function sets the value of forward_tunover and 
            backward_turnover. 
        '''
    
    
    
class Metabolite():
    '''Class containing all information of interest pertaining to a metabolite. 
        
    
    Attributes:
        name: the name of the metabolite
        bigg: BiGG identifier
        kegg: kegg identifier
        turn: turnover number 
        spec: specific activity (turnover number divided by the 
              molecular weight)
        molw: molecular weight 
    '''

    def __init__(self, name, bigg = None, kegg = None, **kwdata):
        '''TODO: Some info on this constructor'''
        self.turnover = None
        self.molecular_weight = None
        self.specific_activity = None
        self.name = name
        if bigg == None:
            self.bigg = ''
        else: 
            self.bigg = bigg
        if kegg == None:
            self.kegg = ''
        else:
            self.kegg = kegg
        if kwdata is not None:
            for identifier, data in kwdata.items(): 
                if identifier == 'turnover':
                    self.turn = data
                if identifier == 'specific activity':
                    self.specific_activity = data
                if identifier == 'molecular weight':
                    self.molecular_weight = data
   
class MetaboliteCandidate(Metabolite):
    '''Metabolite class intended to be used for data selection.
    
    On top of the actual values, selection of the potential metabolites data
        is based on the organism that data is provenant from, and whether the 
        molecule tested was wild-type or not
    ''' 
    
    def __init__(self, name, bigg=None, kegg=None, **kwdata):
        self.organism == None
        self.wild_type == None
        super().__init__(self, name, bigg, kegg, **kwdata)
        if kwdata is not None:
            for identifier, value in kwdata.items():
                if identifier == 'organism':
                    self.organism = value
                if identifier == 'wild-type':
                    self.wild_type = value
    
class DataMissingError(Exception):
    pass