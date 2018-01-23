#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 14:26:03 2018

@author: dimitricoukos
"""

import unittest
import cobra
from DataTreatment import DataType, openJson, matchById, write,\
    Enzyme, Metabolite, MetaboliteCandidate, writeEnzymes

class MatchById(unittest.TestCase):
    '''Contains functions to test Match By Id
       
    Note:
        openJson already checked in TestFileIO.py
    ''' 
    def loadMadeUpModel(path):
        '''Load a cobra model using a structure based on Enzymes and Metabolite
        Candidates
        
        Args: 
            path: filepath to made up mode. 
        Note:
            This function, based on DataTreatment.storeBiGGRepresentation is for 
            testing purposes, and allows to load a made up model.
        '''
        made_up_model = cobra.io.load_json_model(path)
        local_repr = {}
        for reaction in made_up_model.reactions:
            if reaction.id == 'CSND' or reaction.id == 'DHPM1':
                local_repr[reaction.id] = Enzyme(reaction.id)
    
                for reactant in reaction.reactants:
                    local_repr[reaction.id].forward[reactant.name] = Metabolite(
                        reactant.name, bigg=reactant.id)
                for product in reaction.products:
                    local_repr[reaction.id].backward[product.name] = Metabolite(
                        product.name, bigg=product.id)
        return local_repr 
   
    brenda_keggs = openJson('Unit Tests/sample_brenda_keggs.json')
    treated_brenda_output = openJson(
            'Unit Tests/sample_simple_brenda_output.json')
    data_type = DataType.turnover
    potential_updates_dict = openJson(
            'Unit Tests/correct_potential_updates.json')
    simple_test_model = loadMadeUpModel('Unit Tests/simple_test_model.json')
    simple_test_model['CSND'].with_kegg['C00380'] = 'cyt' #reactant, forward
    simple_test_model['CSND'].with_kegg['D00323'] = '5-fluorocyt'# product, backward
    simple_test_model['DHPM1'].with_kegg['C00148'] = 'DL-p'  #product, backward  
    correct_potential_updates = {}
    for reaction in potential_updates_dict:
        correct_potential_updates.update({reaction:Enzyme(reaction)})
        if reaction in brenda_keggs:
            for kegg in brenda_keggs[reaction]: 
                brenda_name = brenda_keggs[reaction][kegg]                
                if kegg in simple_test_model[reaction].with_kegg:
                    if simple_test_model[reaction].with_kegg[kegg] in\
                    simple_test_model[reaction].forward: 
                        if treated_brenda_output[reaction][brenda_keggs[reaction][kegg]]!=[]:
                            name = simple_test_model[reaction].with_kegg[kegg]
                            correct_potential_updates[reaction].forward[name] = []
                            for entry in treated_brenda_output[reaction][brenda_name]:
                                data = {
                                        'organism': entry['organism'], 
                                        'wild-type': entry['wild-type'],
                                        'turnover': entry['turnoverNumber']
                                        }
                                correct_potential_updates[reaction].forward[name].append(
                                        MetaboliteCandidate(brenda_name, data))
                
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
                        correct_potential_updates[reaction].backward[name].append(
                                MetaboliteCandidate(brenda_name, data))
  
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
   
    def test_matchById_unmatched(self):
        '''Tests that the unmatched dict returns correctly. '''
        potential_updates = {}
        unmatched = matchById(potential_updates, 
                              MatchById.brenda_keggs,
                              MatchById.treated_brenda_output, 
                              MatchById.data_type, 
                              MatchById.simple_test_model)
        write('Unit Tests/return_matchById_unmatched.json', unmatched)
   
        #unmatched only needs to return the BRENDA names of the unmatched 
        #metabolites. matchByName will look for the data. (Maintains code 
        #similarity.)
        for reaction in MatchById.correct_unmatched:
            self.assertListEqual(sorted(MatchById.correct_unmatched[reaction]),
                                 sorted(unmatched[reaction]), msg='Unmatched '
                                 'return incorrect.')
       
       
if __name__ == '__main__':
    unittest.main()