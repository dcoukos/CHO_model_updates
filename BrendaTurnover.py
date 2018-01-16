'''
This script allows us to download and organize enzyme
turnover entries for a given map of BiGG IDs and EC
numbers. This data is taken from a JSON file.

Python 2 because SOAPpy is outdated, and recommended
for communicating with BRENDA.
'''
import os
import sys
import string
import hashlib
import json
import SOAPpy
from SOAPpy import SOAPProxy
from BRENDA_extract import getKeggIds

def main():
    BRENDA_parameters = importBrendaParameters()
    BRENDA_output = importBrendaEntries(BRENDA_parameters)
    saveBrendaOutput(BRENDA_output)
    treated_output = treatBrendaOutput(BRENDA_output)
    saveTreatedEntries()


def importBrendaParameters():
    with open('JSONs/BRENDA_parameters.json', 'r') as json_file:
        return json.load(json_file)
    

def importBrendaTurnovers(BRENDA_parameters):
    output = {}
    for ID in BRENDA_parameters:
        EC_number = BRENDA_parameters[ID][0]
        reactants = BRENDA_parameters[ID][2]
        try:
            # in some intranets an issue: how to use a web proxy for WS. Here
            # we assume a set environment variable 'http_proxy'.·
            # This is common in unix environments. SOAPpy does not like
            # a leading 'http://'
            if os.environ.has_key("http_proxy"):
                my_http_proxy=os.environ["http_proxy"].replace("http://","")
            else:
                 my_http_proxy=None

            from SOAPpy import SOAPProxy ## for usage without WSDL file

            endpointURL = "http://www.brenda-enzymes.org/soap/brenda_server.php"
            proxy = SOAPProxy(endpointURL, http_proxy=my_http_proxy)
            #proxy = SOAPProxy(endpointURL)
            password = hashlib.sha256("Feynman").hexdigest()
            parameters = 'cossio@cim.sld.cu,'+password+',ecNumber*'+EC_number
            new_EC = proxy.getEcNumber(parameters)
            if('transferred to' in new_EC):
                new_EC = new_EC.rsplit(' ', 1)[1]
                new_EC = new_EC[:-1]
                EC_number = new_EC
                BRENDA_parameters[ID][0] = new_EC
            output[ID] = proxy.getTurnoverNumber(parameters)
            return output
        except:
            raise
    with open('JSONs/BRENDA_parameters_v2.json', 'w') as outfile:
        json.dump(output, outfile, indent=4)

def saveBrendaOutput(output):
    with open('JSONs/BRENDA_output.json', 'w') as outfile:
        json.dump(output, outfile, indent=4)

def treatBrendaOutput(output):
    '''
    Removes unnecessary parameters from entries and
    checks to see if enzymes characterized were
    wild-type or mutant.
    '''
    treated_output = {}
    treated_output = {ID: [{item.split('*')[0]: item.split('*')[1] 
        for item in entry.split('#') 
        if len(item.split('*')) > 1} 
        for entry in output[ID].split('!')] for ID in output.keys()}

    treated_output_by_substrate = treated_output
    for ID in treated_output:
        new_entry = {}
        for entry in treated_output[ID]:
            for data_point, description in entry.iteritems():
                if data_point == 'substrate':
                    if description in new_entry:
                        new_entry[description].append(entry)
                    else:
                        new_entry[description] = []
        treated_output[ID] = new_entry

    no_data = []

    for ID in treated_output: 
        if output[ID] == '':
            no_turnover_data.append(ID)
        else:
            empty = bool
            for substrate in treated_output[ID]:        
                commentary_treated = False
                wild_type = False
                for entry in treated_output[ID][substrate]:
                    if entry == [] :
                        continue
                    else:
                        for key,value in entry.iteritems():
                            if (key == 'commentary') and 'wild' in value:
                                wild_type = True
                                commentary_treated = True
                            elif (key == 'commentary') and 'mutant' in value:
                                wild_type = False
                                commentary_treated = True
                        print(ID)
                        entry.pop('literature')
                        entry.pop('substrate')
                        entry.pop('ligandStructureId')
                        entry.pop('turnoverNumberMaximum')
                        entry.pop('commentary', 'No comment')
                        if wild_type:
                            entry['wild-type'] = True
                        elif not wild_type and commentary_treated:
                            entry['wild-type'] = False


def saveTreatedEntries():
    with open('JSONs/treated_BRENDA_output.json', 'w') as outfile:
        json.dump(treated_output, outfile, indent=4)


def importKeggIds():
    with open('JSONs/Model_KEGG_IDs.json', 'r') as json_file:
        return json.load(json_file)
    

selectEntries()
                
main()
