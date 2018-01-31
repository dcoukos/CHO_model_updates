'''
  This file contains the necessary suite of tests to evaluate the functioning
  of DataTreatment.py, under normal circumstances.

'''
import unittest

from DataTreatment import getData, openJson, DataType, Organism, \
    MetaboliteCandidate


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


if __name__ == '__main__':
    unittest.main()
