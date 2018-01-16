'''
    This file contains functions pertaining to the treatment of the data that 
    is returned from BRENDA, including functions that find the KEGG Ids of that
    data, as well as functions that organize that data and simplify it in order
    to keep only the necessary information for improving the model. 
'''


import cobra
import json
import requests


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
    def _init_(self, bigg, EC = None, metabolites = {}, forward = {}, backward = {}):
        self.bigg = bigg
        self.EC = EC
        self.metabolites = metabolites
        self.forward = forward
        self.backwards = backwards
        self.with_kegg = {}
    
    def returnWithKegg():
        '''Returns a dict of metabolites by their KEGG ids.'''
        
    def updateKegg():
        '''iterates through all metabolites. If they have a Kegg ID, they are 
        added to '''
        pass
    
    
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

    def _init_(self, name, bigg = None, kegg = None, turnover = None, 
               specific_activity = None, molecular_weight = None):
        self.name = name
        self.bigg = bigg
        self.kegg = kegg
        self.turn = turnover
        self.spec = specific_activity
        self.molw = molecular_weight

def TreatTurnoverData():
    
    eliminateUnmatchable()
    selectBestData()
    
  
    
def organizeForwardsBackwards(Data):
    '''
    '''

def eliminateUnmatchable(brenda_data):
    '''Eliminates entries from BRENDA which cannot be matched to model data. 
    
    The return from BRENDA is quite large, and only useful if pertinent to 
    metabolites in the model. The data must therefore be pruned. When correct
    matches are found they are added to a new data structure, which is also 
    organized for forward and backwards reactions. 
    
    Args:
        brenda_data: This is a dict containing data from BRENDA. Keys are 
                        reaction Ids.
        
    Returns: 
        pruned_data: Data from brenda with unmatchable entries removed. 
    
    Data passed to this function must at least have this format:
        { reaction1:{ metabolite1: ..., 
                      metabolite2: ...,
                      metabolite3: ...,
                     },
            ...     
        }
    
    '''
    
    treated_output = openJson('JSONs/treated_BRENDA_output.json')
    brendaToKegg(treated_output)
    brenda_keggs = correctJson('JSONs/treated_BRENDA_output.json')
    bigg_model = cobra.io.load_json_model('JSONs/BIGG_master_modle.json')
    [directional_model, unmatched] = matchById(brenda_keggs, bigg_model)
    
    
    
def fuzzyMatchNames():
    '''Tries to match metabolite names that cannot be matched via KEGG.'''
   

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

def matchById(brenda_keggs, bigg_model):
    '''Tries to match metabolites by KEGG Id. 
    
    Args:
        brenda_keggs: dict of KEGG codes (keys) and corresponding names (values)
        bigg_model: cobra model of CHO downloaded from BiGG database. 
    
    Return: 
        new_model: dict of model with matched metabolites. 
        no_match: dict of BRENDA metabolites that could not be matched by name. 
    
    TODO: add new_model structure to this docstring
    '''
    bigg_representation = {}
    for reaction in bigg_model.reactions:
        bigg_representation[reaction.name] = {}
        for bigg_model.reactions:
            pass
        '''TODO: finish restructuring of classes: Enzyme and Metabolite, and then 
        restructure the new data structures accordingly. '''
    for reaction in brenda_keggs:
        for metabolite in brenda_keggs[reaction]:
            
        
        
    
def cleanMetaboliteNames(brenda_metabolites):
    '''
        TODO: fill out function
    '''
    return 1
    
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