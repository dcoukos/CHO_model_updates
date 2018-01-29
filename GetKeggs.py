import json
import sys
import requests
import cobra
from multiprocessing import Pool
from bs4 import BeautifulSoup as Soup
from progress.bar import Bar
from DataTreatment import openJson, write


def extractBiggKeggs(reactions, cm_param, p_number):
    """Extracts Kegg metabolite IDs from iCHOv1.xml.

        In order to match turnover data between BRENDA and the BiGG model, we
        need a code that can be applied to both sets of metabolite data.
        Currently the only available piece of data that we can use is the kegg
        identifier. While the chemical translation service has to be used for
        the metabolites from BRENDA, the iCHOv1.xml file contains many kegg
        ids.

        Parameters
        ----------
        reactions : list
            list of cobra.reaction
        cm_param : string
            command line parameter, telling function whether it's running as
                a test.
        p_number: int
            Process number which is running this function
        Returns
        -------
        dict
            dict containing kegg codes extracted by subprocess for metabolites
    """
    if cm_param == 'model':
        handler = open('iCHOv1.xml').read()
    else:
        handler = open('Unit Tests/sample_xml.xml').read()
    soup = Soup(handler, 'xml')
    local_model = {}
    total = len(reactions)
    if p_number == 1:
        bar = Bar('Processing: ' + cm_param, max=total)
    for reaction in reactions:
        react_id = reaction.id
        local_model[react_id] = {}
        local_model[react_id]['reactants'] = {}
        local_model[react_id]['products'] = {}
        if p_number == 1:
            bar.next()
        for reactant in reaction.reactants:
            rid = 'M_' + reactant.id
            local_model[react_id]['reactants'][rid] = []
            species = soup.find('species', id=rid)
            if species:
                links = species.find_all('li')
                no_kegg = True
                for link in links:
                    if 'identifiers.org/kegg.compound' in link['resource']:
                        local_model[react_id]['reactants'][rid].append(
                            link['resource'][37:])
                    elif 'identifiers.org/kegg.drug' in link['resource']:
                        local_model[react_id]['reactants'][rid].append(
                            link['resource'][33:])
                    elif 'identifiers.org/kegg.glycan' in link['resource']:
                        local_model[react_id]['reactants'][rid].append(
                            link['resource'][35:])
                        no_kegg = False
                if no_kegg:
                    local_model[react_id]['reactants'][rid] = None
        for product in reaction.products:
            pid = 'M_' + product.id
            local_model[react_id]['products'][pid] = []
            species = soup.find('species', id=pid)
            if species:
                links = species.find_all('li')
                no_kegg = True
                for link in links:
                    if 'identifiers.org/kegg.compound' in link['resource']:
                        local_model[react_id]['products'][pid].append(
                            link['resource'][37:])
                    elif 'identifiers.org/kegg.drug' in link['resource']:
                        local_model[react_id]['products'][pid].append(
                            link['resource'][33:])
                    elif 'identifiers.org/kegg.glycan' in link['resource']:
                        local_model[react_id]['products'][pid].append(
                            link['resource'][35:])
                        no_kegg = False
                if no_kegg:
                    local_model[react_id]['products'][pid] = None
    return local_model


def getBrendaKeggs(reactions, process):
    '''

    TODO: Make sure this returns the correct data structure for the
    metabolite_no_kegg.
    '''
    # TODO: restructure function to work with Brenda return.
    #       Source of data loss?
    # TODO: make sure the function knows how to deal with empty data.
    # TODO: Rewrite program to call this function with multiprocessing
    total = len(reactions)
    if process == 1:
        bar = Bar('Retrieving kegg codes for brenda ouput: ', max=total)
    metabolite_to_KEGG = {}
    metabolite_no_KEGG = {}
    for bigg_id in reactions:
        bar.next()
        metabolite_to_KEGG[bigg_id] = {}
        metabolite_no_KEGG[bigg_id] = []
        for metabolite in reactions[bigg_id]:
            request_counter = 0
            while request_counter < 3:
                try:
                    if ' - reduced' in metabolite:
                        metabolite = 'reduced ' + metabolite[:-10]
                    cts_output = requests.get(
                                    'http://cts.fiehnlab.ucdavis.edu/'
                                    'service/convert/Chemical%20Name/'
                                    'KEGG/' + metabolite)
                    kegg_id = str(json.loads(cts_output.text)[0]['result'][0])
                    metabolite_to_KEGG[bigg_id][kegg_id] = metabolite
                    request_counter = 3
            # TODO: What type of error is expected to be thrown here?
            # This may be a source of missing entries.
                except IndexError:  # Just a guess!
                    if request_counter == 2:
                        metabolite_no_KEGG[bigg_id].append(metabolite)
                    request_counter = request_counter + 1

    return metabolite_to_KEGG, metabolite_no_KEGG


if __name__ == '__main__' and len(sys.argv) > 1:
    if sys.argv[1] == 'model' or sys.argv[1] == 'test-model':
        print('Opening model...')
        model = cobra.io.read_sbml_model('iCHOv1.xml')
        reactions1 = []
        reactions2 = []
        reactions3 = []
        reactions4 = []
        local_model = {}
        counter = 0
        for reaction in model.reactions:
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
            lm_1 = pool.apply_async(extractBiggKeggs, (reactions1,
                                    sys.argv[1], 1,))
            print('Process on core 1 started')
            lm_2 = pool.apply_async(extractBiggKeggs, (reactions2,
                                    sys.argv[1], 2,))
            print('Process on core 2 started')
            lm_3 = pool.apply_async(extractBiggKeggs, (reactions3,
                                    sys.argv[1], 3,))
            print('Process on core 3 started')
            lm_4 = pool.apply_async(extractBiggKeggs, (reactions4,
                                    sys.argv[1], 4,))
            print('Process on core 4 started')
            print('Porting model xml file to processes...')
            pool.close()
            pool.join()
        local_model.update(lm_1.get())
        local_model.update(lm_2.get())
        local_model.update(lm_3.get())
        local_model.update(lm_4.get())
        if sys.argv[1] == 'model':
            write('JSONs/iCHOv1_keggs.json', local_model)
        else:
            write('Unit Tests/iCHOv1_keggs_test.json', local_model)
    elif sys.argv[1] == 'brenda-keggs':
        print('Opening brenda output...')
        treated_brenda_output = openJson('JSONs/treated_BRENDA_output.json')
        brenda_keggs = {}
        no_keggs = {}
        reactions1 = {}
        reactions2 = {}
        reactions3 = {}
        reactions4 = {}
        counter = 0
        for reaction in treated_brenda_output:
            if counter % 4 == 0:
                reactions1[reaction] = treated_brenda_output[reaction]
            if counter % 4 == 1:
                reactions2[reaction] = treated_brenda_output[reaction]
            if counter % 4 == 2:
                reactions3[reaction] = treated_brenda_output[reaction]
            if counter % 4 == 3:
                reactions4[reaction] = treated_brenda_output[reaction]
            counter = counter + 1
        with Pool(processes=4) as pool:
            r_set1 = pool.apply_async(getBrendaKeggs, (reactions1, 1,))
            print('Process on core 1 started')
            r_set2 = pool.apply_async(getBrendaKeggs, (reactions2, 2,))
            print('Process on core 2 started')
            r_set3 = pool.apply_async(getBrendaKeggs, (reactions3, 3,))
            print('Process on core 3 started')
            r_set4 = pool.apply_async(getBrendaKeggs, (reactions4, 4,))
            print('Process on core 4 started')
            pool.close()
            pool.join()
        brenda_keggs.update(r_set1[0].get())
        brenda_keggs.update(r_set2[0].get())
        brenda_keggs.update(r_set3[0].get())
        brenda_keggs.update(r_set4[0].get())

        no_keggs.update(r_set1[1].get())
        no_keggs.update(r_set2[1].get())
        no_keggs.update(r_set3[1].get())
        no_keggs.update(r_set4[1].get())

        write('JSONs/brenda_keggs.json', brenda_keggs)
        write('JSONs/brenda_no_keggs.json', no_keggs)
