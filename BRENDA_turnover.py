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

def main():
    BRENDA_parameters = import_BRENDA_parameters()
    BRENDA_output = import_BRENDA_entries(BRENDA_parameters)
    save_BRENDA_output()
    treated_BRENDA_output = treat_BRENDA_output(BRENDA_output)
    save_treated_entries()


def import_BRENDA_parameters():
    with open('BRENDA_parameters.json', 'r') as json_file:
        return json.load(json_file)

def import_BRENDA_entries(BRENDA_parameters):
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

            from SOAPpy import WSDL ## for extracting the URL of the endpoint (server script) from the WSDL file
            from SOAPpy import SOAPProxy ## for usage without WSDL file

            endpointURL = "http://www.brenda-enzymes.org/soap/brenda_server.php"
            proxy = SOAPProxy(endpointURL, http_proxy=my_http_proxy)
            #proxy = SOAPProxy(endpointURL)
            password = hashlib.sha256("Feynman").hexdigest()
            parameters = 'cossio@cim.sld.cu,'+password+',ecNumber*'+EC_number
            output[ID] = proxy.getTurnoverNumber(parameters)
        except:
            raise

save_BRENDA_output():
    with open('BRENDA_output.json', 'w') as outfile:
        json.dump(output, outfile, indent=4)

treat_BRENDA_output(BRENDA_output):
'''
    Removes unnecessary parameters from entries and
    checks to see if enzymes characterized were
    wild-type or mutant.
'''
    treated_output = {}
    treated_output = {ID: [{item.split('*')[0]: item.split('*')[1] for item in entry.split('#')
         if len(item.split('*')) > 1}
         for entry in output[ID].split('!')] for ID in output.keys()}

         no_turnover_data = []

    for ID in treated_output:
        if output[ID] == '':
            no_turnover_data.append(ID)
        else:
            for entry in treated_output[ID]:
                commentary_treated = False
                wild_type = False
                for data_point, description in entry.iteritems():
                    if (data_point == 'commentary') and 'wild' in description:
                        wild_type = True
                        commentary_treated = True
                    elif (data_point == 'commentary') and 'mutant' in description:
                        wild_type = False
                        commentary_treated = True
                print ID
                entry.pop('literature')
                entry.pop('ligandStructureId')
                entry.pop('turnoverNumberMaximum')
                entry.pop('commentary', 'No comment')
                if wild_type:
                    entry['wild-type'] = True
                elif not wild_type and commentary_treated:
                    entry['wild-type'] = False
def save_treated_entries():
    with open('treated_BRENDA_output.json', 'w') as outfile:
        json.dump(treated_output, outfile, indent=4)

main()
