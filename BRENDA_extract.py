'''
This script collects the KEGG IDs of metabolites present within the iCHOv1
model downloaded from BiGG.uscd.edu. The purpose of this is so that the 
metabolites can be compared with reactions available from BRENDA, in order 
to obtain the correct information for a given EC number. 


'''
import requests
import json
import cobra
from bs4 import BeautifulSoup as Soup


'''
Loading the xml into BeautifulSoup allows me to extract data present in the 
xml (namely the EC numbers for the reactions) which is lost when loaded as a 
a cobra model. 
'''
handler = open("iCHOv1.xml").read() 
bigg_model = cobra.io.load_json_model("BIGG_master_modle.json")
soup = Soup(handler, 'xml')
reactions = soup.find_all('reaction')

BiGG_to_EC = {}
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
        BiGG_to_EC[BiGG_ID]= [EC_number, name, reactants]
        for react in reactants:
            reactant_dict[react] = ""
            
reactants = list(reactant_dict.keys())
reactant_to_KEGG = {}
reactants_no_KEGG = []
reactant_counter = 0

for reactant in reactants[reactant_counter:]:
    try:
        if " - reduced" in reactant:
            reactant = "reduced " + reactant[:-10]
        cts_output = requests.get("http://cts.fiehnlab.ucdavis.edu/service/convert/Chemical%20Name/KEGG/"+reactant)
        reactant_to_KEGG[reactant] = str(json.loads(cts_output.text)[0]['result'][0])
        
        reactant_counter = reactant_counter + 1
        print("reactant " + str(reactant_counter) + " : " + reactant)
    except: 
        try:
            cts_output = requests.get("http://cts.fiehnlab.ucdavis.edu/service/convert/Chemical%20Name/KEGG/"+reactant)
            if json.loads(cts_output.text)[0]['result']:
                reactant_to_KEGG[reactant] = str(json.loads(cts_output.text)[0]['result'][0])
            else: 
                print("Reactant " + reactant +' has no KEGG')
                reactants_no_KEGG.append(reactant)
            reactant_counter = reactant_counter + 1
            print("reactant " + str(reactant_counter) + " : " + reactant)
        except:
            print("Incorrect data retrieval for reactant: " + reactant)
            print("Last correctly saved index: " + str(reactant_counter-1))
            continue