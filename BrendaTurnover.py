'''
This script allows us to download and organize enzyme
turnover entries for a given map of BiGG IDs and EC
numbers. This data is taken from a JSON file.

Python 2 because SOAPpy is outdated, and recommended
for communicating with BRENDA.
'''
import json
import cobra_services as CS
from BeautifulSoup import BeautifulSoup as Soup


def main():
    ec_numbers = getEcNumbers('iCHOv1.xml')
    treated_output = importBrendaTurnovers(ec_numbers)
    simplified_output = simplifyBrendaOutput(treated_output)
    saveTreatedEntries(simplified_output)


def importBrendaParameters():
    with open('JSONs/BRENDA_parameters.json', 'r') as json_file:
        return json.load(json_file)


def importBrendaTurnovers(ec_numbers):
    output = {}
    for ID in ec_numbers:
        # TODO : include this as a parameter somewhere to remove acc. info
        raw_output = CS.brenda(ec_numbers[ID], 'cossio@cim.sld.cu', 'Feynman')
        output[ID] = CS.parse_brenda_raw_output(raw_output)

    return output


def saveBrendaOutput(output):
    with open('JSONs/BRENDA_output.json', 'w') as outfile:
        json.dump(output, outfile, indent=4)


def simplifyBrendaOutput(output):
    '''
    Removes unnecessary parameters from entries and
    checks to see if enzymes characterized were
    wild-type or mutant.
    '''
    new_output = {}
    for ID in output:
        new_entry = {}
        for entry in output[ID]:
            for data_point, description in entry.iteritems():
                if data_point == 'substrate':
                    if description in new_entry:
                        new_entry[description].append(entry)
                    else:
                        new_entry[description] = []
        output[ID] = new_entry

        if output[ID] != '':
            for substrate in output[ID]:
                commentary_treated = False
                wild_type = False
                for entry in output[ID][substrate]:
                    if entry == []:
                        continue
                    else:
                        for key, value in entry.iteritems():
                            if (key == 'commentary') and 'wild' in value:
                                wild_type = True
                                commentary_treated = True
                            elif (key == 'commentary') and 'wt' in value:
                                wild_type = True
                                commentary_treated = True
                            elif (key == 'commentary') and 'muta' in value:
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
            new_output[ID] = output[ID]
    return new_output


def saveTreatedEntries(output):
    with open('JSONs/treated_BRENDA_output.json', 'w') as outfile:
        json.dump(output, outfile, indent=4)


def getEcNumbers(file):
    handler = open(file).read()
    soup = Soup(handler, 'xml')
    reactions = soup.find_all('reaction')
    ec_numbers = {}

    for reaction in reactions:
        links = reaction.find_all("li")
        BiGG_ID = ""
        EC_number = ""

        for link in links:
            if "identifiers.org/ec-code" in link['resource']:
                EC_number = link['resource'][31:]
            elif "identifiers.org/bigg.reaction" in link['resource']:
                BiGG_ID = reaction['id']
                if 'R_' in BiGG_ID:
                    BiGG_ID = BiGG_ID[2:]

        if EC_number and BiGG_ID:
            # CHANGED: now we just save the EC number, since the names of the
                # enzymes and metabolites were never used for querying Brenda.
            ec_numbers[BiGG_ID] = EC_number

    return ec_numbers


if __name__ == '__main__':
    main()
