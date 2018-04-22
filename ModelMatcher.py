'''
    This file contains a script which matches the new bigg reaction names in
    the iCHO model downloaded from BiGG to the older names in the iCHOv1_K1
    model from the paper.

    This allows us match data form BRENDA by reaction to the iCHOv1_K1 model.
    This is because the xml for this model did not contain EC codes, while the
    xml from BiGG did.
'''
import sys
import json
import cobra
import argparse
from statistics import median
from multiprocessing.pool import Pool
from bs4 import BeautifulSoup as Soup
from progress.bar import Bar


def openJson(path):
    '''Shortened call to open JSON files.'''
    with open(path) as r:
        return json.load(r)


def write(path, data):
    '''Shortened call to JSON dumps with indent = 4'''
    with open(path, 'w') as wr:
        json.dump(data, wr, indent=4)


def MatchReaction(bid, reactants, products, reactions):
    ''' PSEUDOCODE:

    If for every metabolite there is a matching link, return true. Else
    return false.'''

    r_names = []
    for reactant in reactants:
        r_names.append('M_' + reactant.id)
    p_names = []
    for product in products:
        p_names.append('M_' + product.id)
    for reaction in reactions:
        try:
            reactant_list = reaction.find('listOfReactants')
            products_list = reaction.find('listOfProducts')
            k1_reactants = reactant_list.find_all('speciesReference')
            k1_products = products_list.find_all('speciesReference')
            rids = []
            pids = []
            for reactant in k1_reactants:
                rids.append(reactant['species'])
            for product in k1_products:
                pids.append(product['species'])
            if rids == r_names and pids == p_names:
                return reaction['id']
        except AttributeError:
            pass
    return None


def mainSubProcess(reactions, process):
    k1_model = cobra.io.read_sbml_model('iCHOv1_K1_final.xml')
    handler = open('iCHOv1.xml').read()
    handler = open('iCHOv1_K1_final.xml').read()
    k1_xml = Soup(handler, 'xml')
    total = len(reactions)
    if process == 1:
        bar = Bar('Matching models ', max=total)
    xml_k1_reactions = k1_xml.find_all('reaction')
    v1_to_k1 = {}
    for reaction in reactions:
        if reaction.id not in k1_model.reactions:
            reactants_bigg = reaction.reactants
            products_bigg = reaction.products
            v1_to_k1[reaction.id] = MatchReaction(reaction.id, reactants_bigg,
                                                  products_bigg,
                                                  xml_k1_reactions)
        if process == 1:
            bar.next()
    return v1_to_k1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='This is a script to match kegg codes and turnovers'
        'between the CHO and CHO K1 models. ')
    parser.add_argument('-s', '--search', help='Match bigg ids',
                        action='store_true')
    parser.add_argument('-m', '--match', help='Match turnovers',
                        action='store_true')
    args = parser.parse_args()

    if args.search:
        print('Loading models....')
        bigg_model = cobra.io.read_sbml_model('iCHOv1.xml')

        v1_to_k1 = {}
        reactions1 = []
        reactions2 = []
        reactions3 = []
        reactions4 = []
        counter = 0
        for reaction in bigg_model.reactions:
            if counter % 4 == 0:
                reactions1.append(reaction)
            if counter % 4 == 1:
                reactions2.append(reaction)
            if counter % 4 == 2:
                reactions3.append(reaction)
            if counter % 4 == 3:
                reactions4.append(reaction)
            counter = counter + 1
        with Pool(processes=4) as pool:
            lm_1 = pool.apply_async(mainSubProcess, (reactions1, 1,))
            print('Process on core 1 started')
            lm_2 = pool.apply_async(mainSubProcess, (reactions2, 2,))
            print('Process on core 2 started')
            lm_3 = pool.apply_async(mainSubProcess, (reactions3, 3,))
            print('Process on core 3 started')
            lm_4 = pool.apply_async(mainSubProcess, (reactions4, 4,))
            print('Process on core 4 started')
            print('Porting k1 xml file to processes...')
            pool.close()
            pool.join()
        v1_to_k1.update(lm_1.get())
        v1_to_k1.update(lm_2.get())
        v1_to_k1.update(lm_3.get())
        v1_to_k1.update(lm_4.get())

        write('JSONs/v1_to_k1.json', v1_to_k1)
    elif args.match:
        conv_dict = openJson('JSONs/v1_to_k1.json')
        k1_to_v1 = {}
        updates = openJson('JSONS/model_updates.json')
        k1_model = cobra.io.read_sbml_model('iCHOv1_K1_final.xml')
        k1_updates = {}
        turn = []
        for met in updates:
            turn.extend(updates[met].values())
        median = median(turn)
        for v1 in conv_dict:
            k1_id = conv_dict[v1]
            if k1_id is not None:
                k1_to_v1[k1_id[2:]] = v1

        for rxn in k1_model.reactions:
            rxn_id = rxn.id
            if rxn_id in updates:
                k1_updates[rxn_id] = updates[rxn_id]
            elif rxn_id in k1_to_v1:
                k1_updates[rxn_id] = updates[k1_to_v1[rxn_id]]
            else:
                k1_updates[rxn_id] = {}
                k1_updates[rxn_id]['forward'] = median
                k1_updates[rxn_id]['backward'] = median
        write('JSONs/k1_updates_ave.json', k1_updates)
