'''
This script collects the KEGG IDs of metabolites present within the iCHOv1
model downloaded from BiGG.uscd.edu. The purpose of this is so that the
metabolites can be compared with reactions available from BRENDA, in order
to obtain the correct information for a given EC number.


'''
import requests
import json
import cobra
from BeautifulSoup import BeautifulSoup as Soup


def main():
    writeBrendaParameters()
    writeKeggIds()
'''
Loading the xml into BeautifulSoup allows me to extract data present in the
xml (namely the EC numbers for the reactions) which is lost when loaded as a
a cobra model.
'''

'''This is more useful than loading the iCHOv1 as a cobra model, because we are
   only interested in the reactions that have an accompanying EC number.'''

def writeBrendaParameters():
    '''
    Obtains BRENDA parameters from the iCHOv1_xml file. Only entries that have an
    identifiable EC number and BiGG ID are treated and saved.
    The EC Number, enzyme name, and metabolite names are saved under a BiGG Id to
    a JSON file: BRENDA_parameters.json
    '''
    bigg_model = cobra.io.load_json_model('BIGG_master_modle.json')
    BRENDA_parameters = getBrendaParametersAndReactants(bigg_model)[0]
    writeToJson('BRENDA_parameters.json', BRENDA_parameters)

def writeKeggIds():
    '''
    Obtains KEGG ids for the reactant names given in the bigg model for reactions
    in the iCHO_v1.xml.
    Metabolite names and the KEGG Ids are saved to a json file: Model_KEGG_IDs.json
    '''
    bigg_model = cobra.io.load_json_model('BIGG_master_modle.json')
    reactants = getBrendaParametersAndReactants(bigg_model)[1]
    model_KEGG_IDs = getKeggIds(reactants)
    writeToJson('Model_KEGG_IDs.json', model_KEGG_IDs)

def getKeggIds(reactants):

    reactant_to_KEGG = {}
    reactant_no_KEGG = []
    reactant_counter = 0

    for reactant in reactants[reactant_counter:]:
        try:
            if " - reduced" in reactant:
                reactant = "reduced " + reactant[:-10]
            cts_output = requests.get("http://cts.fiehnlab.ucdavis.edu/service/convert/Chemical%20Name/KEGG/"+reactant)
            reactant_to_KEGG[str(json.loads(cts_output.text)[0]['result'][0])] = reactant
            reactant_counter = reactant_counter + 1
            print('Reactant '+ reactant + ' added to KEGG')
        except:
            #Sometimes the request glitches, so we try again.
            try:
                cts_output = requests.get("http://cts.fiehnlab.ucdavis.edu/service/convert/Chemical%20Name/KEGG/"+reactant)
                if json.loads(cts_output.text)[0]['result']:
                    reactant_to_KEGG[reactant] = str(json.loads(cts_output.text)
                                                    [0]['result'][0])

                else:
                    reactant_no_KEGG.append(reactant)
                    print('Reactant '+ reactant + ' not added to KEGG')
                reactant_counter = reactant_counter + 1
            except:
                continue
    return [reactant_to_KEGG, reactant_no_KEGG]


#def queryBrenda():
  #ERROR: list has no attribute items.
def writeToJson(filename, data):
    if filename == 'Model_KEGG_IDs.json':
        inv_data = {v: k for k, v in data.items()}
        with open(filename, 'w') as outfile:
            json.dump(inv_data, outfile, indent=4)
    else:
        with open(filename, 'w') as outfile:
            json.dump(data, outfile, indent=4)
def getBrendaParametersAndReactants(bigg_model):
    handler = open('iCHOv1.xml').read()
    soup = Soup(handler, 'xml')
    reactions = soup.find_all('reaction')
    BRENDA_parameters = {}
    reactant_dict = {}

    for reaction in reactions:
        links = reaction.find_all("li")
        BiGG_ID = ""
        EC_number = ""
        reactants = []
        name = ""


        for link in links:
            if "identifiers.org/ec-code" in link['resource']:
                EC_number = link['resource'][31:]
                name = reaction['name']
            elif "identifiers.org/bigg.reaction" in link['resource']:
                BiGG_ID = reaction['id']
                if 'R_' in BiGG_ID:
                    BiGG_ID = BiGG_ID[2:]

        if EC_number and BiGG_ID:
            try:   #Links in xml downloaded from BiGG that do not correspond to actual addresses in their database.
                for reactant in bigg_model.reactions.get_by_id(BiGG_ID).reactants:
                    reactants.append(reactant.name)
            except:
                pass
            BRENDA_parameters[BiGG_ID]= [EC_number, name, reactants]
            for react in reactants:
                reactant_dict[react] = ""

    reactants = list(reactant_dict.keys())
    return [BRENDA_parameters, reactants]
