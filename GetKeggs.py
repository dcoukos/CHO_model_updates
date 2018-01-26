import sys
import json
import cobra
from multiprocessing import Pool
from bs4 import BeautifulSoup as Soup
from progress.bar import Bar


def write(path, data):
    '''Shortened call to JSON dumps with indent = 4'''
    with open(path, 'w') as wr:
        json.dump(data, wr, indent=4)


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
                    if 'identifiers.org/kegg' in link['resource']:
                        local_model[react_id]['reactants'][rid].append(
                            link['resource'][37:])
                        no_kegg = False
                if no_kegg:
                    local_model[react_id]['reactants'][rid] = None
        for product in reaction.products:
            pid = 'M_' + product.id
            species = soup.find('species', id=pid)
            if species:
                links = species.find_all('li')
                no_kegg = True
                for link in links:
                    if 'identifiers.org/kegg' in link['resource']:
                        local_model[react_id]['products'][pid] = \
                            link['resource'][37:]
                        no_kegg = False
                if no_kegg:
                    local_model[react_id]['products'][pid] = None
    return local_model


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
            print('Process 1 started')
            lm_2 = pool.apply_async(extractBiggKeggs, (reactions2,
                                    sys.argv[1], 2,))
            print('Process 2 started')
            lm_3 = pool.apply_async(extractBiggKeggs, (reactions3,
                                    sys.argv[1], 3,))
            print('Process 3 started')
            lm_4 = pool.apply_async(extractBiggKeggs, (reactions4,
                                    sys.argv[1], 4,))
            print('Process 4 started')
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
