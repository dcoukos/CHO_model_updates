'''
    This file contains functions pertaining to the treatment of the data that 
    is returned from BRENDA, including functions that find the KEGG Ids of that
    data, as well as functions that organize that data and simplify it in order
    to keep only the necessary information for improving the model. 
'''


import cobra
import json
from BRENDA_extract import getKeggIds




def BrendaToKegg(): 
    '''
    This is simple script which writes two json files based on information
    retrieved from the chemical translation service (http://cts.fiehnlab.ucdavis.edu/).
    Using metabolite names retrieved from BRENDA for the turnover rate of enzymes, the
    translation service is queried to see if an idenitifiable KEGG ID is available for
    that enzyme. This aids in the idenitification of metabolites present in our CHO
    model from BiGG (http://bigg.ucsd.edu/models/iCHOv1). This requires a lot of time,
    as querying is a slow process. Therefore, this script autoupdates two json files.
    One contains the output of the translation service for metabolites with a KEGG Id,
    and the second contains the name of metabolites for which an Id was not found.
'''
    bigg_model = cobra.io.load_json_model('JSON/BIGG_master_modle.json')
    
    def getModelKegg():
        with open('JSON/Model_KEGG_IDs.json', 'r') as infile:
                return json.load(infile)
    
    with open('JSON/treated_BRENDA_output.json', 'r') as json_file:
        treated_output =  json.load(json_file)
    
    no_data = []
    selected_entries = {}
    model_KEGG_ids = getModelKegg()
    
    for reaction in treated_output:
        [brenda_kegg_ids, brenda_no_kegg] = getKeggIds(treated_output[reaction].keys())
        kegg_data = {}
        no_kegg_data = []
        try:
            with open('JSON/BRENDA_KEGG_IDs.json') as kegg:
                kegg_data_r = json.load(kegg)
        except ValueError:
            kegg_data_r = {}
        try:
            with open('JSON/BRENDA_no_KEGG.json') as no_kegg:
                no_kegg_data_r = json.load(no_kegg)
        except ValueError:
            no_kegg_data_r = {}
        for metabolite in brenda_kegg_ids:
            kegg_data.update({metabolite:brenda_kegg_ids[metabolite]})
    
    
        for metabolite in brenda_no_kegg:
            no_kegg_data.append(metabolite)
    
    
        kegg_data_r.update({reaction:kegg_data})
        no_kegg_data_r.update({reaction:no_kegg_data})
        with open('JSON/BRENDA_KEGG_IDs.json', 'w') as kegg:
            json.dump(kegg_data_r, kegg, indent = 4)
    
        with open('JSON/BRENDA_no_KEGG.json', 'w') as no_kegg:
            json.dump(no_kegg_data_r, no_kegg, indent = 4)
