'''
  This file contains the necessary suite of tests to evaluate the functioning
  of DataTreatment.py, under normal circumstances.

'''
import unittest
import cobra
from DataTreatment import getData, openJson, DataType, Organism, \
    MetaboliteCandidate, Enzyme, Metabolite, matchById, writeEnzymes


class TestGetData(unittest.TestCase):
    '''contains functions for testing getData function.'''

    def test_getData(self):
        '''Ensure getData can extract and return correct data as a list of dict
        '''
        sample_reaction = openJson('Unit Tests/sample_brenda_reaction.json')
        expected_return = openJson('Unit Tests/expected_getData_return.json')

        returned_data = getData(sample_reaction, 'thioredoxin',
                                DataType.turnover)
        for index, entry in enumerate(returned_data):
            self.assertDictEqual(returned_data[index],
                                 expected_return[index])


class TestReturnAttributes(unittest.TestCase):

    def test_returnAttributes(self):
        '''Tests that Metabolite.returnAttributes returns correct values. '''
        data = {
            'wild-type': True,
            'organism': Organism.mouse,
            'turnover': 7
        }
        correct_return = {
            'turnover': 7,
            'name': 'metabolite',
            'bigg': 'biggly',
            'kegg': None,
            'specific activity': None,
            'molecular weight': None
        }
        candidate = MetaboliteCandidate('metabolite', 'biggly', **data)
        candictionary = candidate.returnAttributes()
        self.assertDictEqual(candictionary, correct_return)


class TestApplyBestData(unittest.TestCase):
    '''Methods for testing applyHighestTurnover'''

    def test_applyBestData(self):
        '''
            Expected Behavior: expects the
        '''


class MatchById(unittest.TestCase):
    # TODO: make sure this test is written correctly.
    '''Contains functions to test Match By Id
    Note:
        openJson already checked in TestFileIO.py.
    '''

    def loadMadeUpModel(path):
        '''Load a cobra model using a structure based on Enzymes and Metabolite
        Candidates

        Args:
            path: filepath to made up mode.
        Note:
            This function, based on DataTreatment.storeBiGGRepresentation is
            for testing purposes, and allows to load a made up model.
        '''
        made_up_model = cobra.io.load_json_model(path)
        local_repr = {}
        for reaction in made_up_model.reactions:
            if reaction.id == 'CSND' or reaction.id == 'DHPM1':
                local_repr[reaction.id] = Enzyme(reaction.id)

                for reactant in reaction.reactants:
                    local_repr[reaction.id].forward[reactant.name] = \
                        Metabolite(reactant.name, bigg=reactant.id)
                for product in reaction.products:
                    local_repr[reaction.id].backward[product.name] = \
                        Metabolite(product.name, bigg=product.id)
        return local_repr

    brenda_keggs = openJson('Unit Tests/sample_brenda_keggs.json')
    treated_brenda_output = openJson(
            'Unit Tests/sample_simple_brenda_output.json')
    data_type = DataType.turnover
    potential_updates_dict = openJson(
            'Unit Tests/correct_potential_updates.json')
    simple_test_model = loadMadeUpModel('Unit Tests/simple_test_model.json')
    simple_test_model['CSND'].with_kegg['C00380'] = 'cyt'
    simple_test_model['CSND'].with_kegg['D00323'] = '5-fluorocyt'
    simple_test_model['DHPM1'].with_kegg['C00148'] = 'DL-p'
    correct_potential_updates = {}
    for reaction in potential_updates_dict:
        correct_potential_updates.update({reaction: Enzyme(reaction)})
        if reaction in brenda_keggs:
            for kegg in brenda_keggs[reaction]:
                brenda_name = brenda_keggs[reaction][kegg]
                if kegg in simple_test_model[reaction].with_kegg:
                    if simple_test_model[reaction].with_kegg[kegg] in\
                            simple_test_model[reaction].forward:
                        if treated_brenda_output[reaction][
                                brenda_keggs[reaction][kegg]] != []:
                            name = simple_test_model[reaction].with_kegg[kegg]
                            correct_potential_updates[reaction].forward[
                                name] = []
                            for entry in treated_brenda_output[reaction][
                                    brenda_name]:
                                data = {
                                        'organism': entry['organism'],
                                        'wild-type': entry['wild-type'],
                                        'turnover': entry['turnoverNumber']
                                        }
                                correct_potential_updates[reaction].forward[
                                    name].append(MetaboliteCandidate(
                                        brenda_name, data))

                    elif simple_test_model[reaction].with_kegg[kegg] in\
                            simple_test_model[reaction].backward:
                            name = simple_test_model[reaction].with_kegg[kegg]
                    correct_potential_updates[reaction].backward[name] = []
                    for entry in treated_brenda_output[reaction][brenda_name]:
                        data = {
                                'organism': entry['organism'],
                                'wild-type': entry['wild-type'],
                                'turnover': entry['turnoverNumber']
                                }
                        correct_potential_updates[reaction].backward[
                            name].append(MetaboliteCandidate(
                                brenda_name, data))

    correct_unmatched = openJson('Unit Tests/correct_unmatched.json')

    def test_matchById_potential_updates(self):
        '''matchByName should match BiGG metabolites with BRENDA metabolites
        given a file containing their respective KEGG Ids.
        '''
        potential_updates = {}
        matchById(potential_updates, MatchById.brenda_keggs,
                  MatchById.treated_brenda_output, MatchById.data_type,
                  MatchById.simple_test_model)
        writeEnzymes('Unit Tests/return_matchById_potential_updates.json',
                     potential_updates)
        potential_updates_as_dict = {}
        correct_potential_updates_as_dict = {}
        for enzyme in potential_updates:
            potential_updates_as_dict[enzyme] = \
                potential_updates[enzyme].getDict()
            for enzyme in MatchById.correct_potential_updates:
                correct_potential_updates_as_dict[enzyme] = \
                    MatchById.correct_potential_updates[enzyme].getDict()

        self.assertDictEqual(potential_updates_as_dict,
                             correct_potential_updates_as_dict,
                             msg='Potential updates incorrect.')


if __name__ == '__main__':
    unittest.main()
