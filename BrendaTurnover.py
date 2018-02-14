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
import cobra_services as CS
from SOAPpy import SOAPProxy
from BRENDA_extract import getKeggIds


def main():
    BRENDA_parameters = importBrendaParameters()
    treated_output = importBrendaEntries(BRENDA_parameters)
    simplified_output = simplifyBrendaOutput(treated_output)
    saveTreatedEntries(simplified_output)


def importBrendaParameters():
    with open('JSONs/BRENDA_parameters.json', 'r') as json_file:
        return json.load(json_file)


def importBrendaTurnovers(BRENDA_parameters):
    output = {}
    for ID in BRENDA_parameters:
        EC_number = BRENDA_parameters[ID][0]
        # TODO : include this as a parameter somewhere to remove acc. info
        raw_output = CS.brenda(EC_number, 'cossio@cim.sld.cu', 'Feynman')
        treated_output = CS.parse_brenda_raw_output(raw_output)

    with open('JSONs/BRENDA_parameters_v2.json', 'w') as outfile:
        json.dump(output, outfile, indent=4)

def saveBrendaOutput(output):
    with open('JSONs/BRENDA_output.json', 'w') as outfile:
        json.dump(output, outfile, indent=4)

def simplifyBrendaOutput(output):
    '''
    Removes unnecessary parameters from entries and
    checks to see if enzymes characterized were
    wild-type or mutant.
    '''
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


def saveTreatedEntries(output):
    with open('JSONs/treated_BRENDA_output.json', 'w') as outfile:
        json.dump(output, outfile, indent=4)


if __name__ == '__main__':
    main()
