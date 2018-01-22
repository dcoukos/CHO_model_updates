#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 17:20:19 2018

@author: dimitricoukos
"""

import unittest
#import json
from DataTreatment import  Enzyme, MetaboliteCandidate
    
class CorrectConstruction(unittest.TestCase):
    '''Class contains methods to verify correct construction of classes from 
    DataTreatment.py'''
    
    def test_buildExpectedEnzyme(self):
        new_enzyme = Enzyme('biggly', '2.2.2.2b')
        self.assertEqual('biggly', new_enzyme.bigg)
        self.assertEqual('2.2.2.2b', new_enzyme.EC)
        self.assertEqual({}, new_enzyme.forward)        
        self.assertEqual({}, new_enzyme.backward)
        self.assertEqual({}, new_enzyme.with_kegg)
        self.assertEqual(None, new_enzyme.forward_turnover)
        self.assertEqual(None, new_enzyme.backward_turnover)
        self.assertEqual(None, new_enzyme.has_f_wild_type)
        self.assertEqual(None, new_enzyme.has_b_wild_type)
        
    def test_buildExpectedMetaboliteCandidate(self):
        metabolite_data_1 = {
                'turnover': 345,
                'specific activity': 2.1,
                'molecular weight': 234,
                'organism': 'Homo sapiens',
                'wild-type': False
                }
        #test order
        metabolite_data_2 = {
                'turnover': 345,
                'wild-type': False,
                'molecular weight': 234,
                'organism': 'Homo sapiens',              
                'specific activity': 2.1
                }
        #test missing data
        metabolite_data_3 = {
            'turnover': 345,
            'organism': 'Homo sapiens',
            'wild-type': False
            }
         #test incorrect key
        metabolite_data_4 = {
            'turnover': 345,
            'spasmatic activity': 2.1,
            'molecular party': 234,
            'organism': 'Homo sapiens',
            'wild-one': True
            }        
    
        new_mc_1 = MetaboliteCandidate('a nongeneric name', metabolite_data_1)
        new_mc_2 = MetaboliteCandidate('a nongeneric name', metabolite_data_2)
        new_mc_3 = MetaboliteCandidate('a nongeneric name', 'bigg!', 
                                       '2.3.4', **metabolite_data_3)
        new_mc_4 = MetaboliteCandidate('a nongeneric name', metabolite_data_4)
        
        mc_list_1 = [new_mc_1.name, 
                     new_mc_1.bigg,
                     new_mc_1.kegg,
                     new_mc_1.molecular_weight,
                     new_mc_1.turnover,
                     new_mc_1.organism,
                     new_mc_1.specific_activity,
                     new_mc_1.wild_type]
        mc_list_2 = [new_mc_2.name, 
                     new_mc_2.bigg,
                     new_mc_2.kegg,
                     new_mc_2.molecular_weight,
                     new_mc_2.turnover,
                     new_mc_2.organism,
                     new_mc_2.specific_activity,
                     new_mc_2.wild_type]
        
        self.assertListEqual(mc_list_1,mc_list_2)
        self.assertEqual('a nongeneric name', new_mc_3.name)
        self.assertEqual('bigg!', new_mc_3.bigg)
        self.assertEqual('2.3.4', new_mc_3.kegg)
        self.assertEqual(None, new_mc_3.molecular_weight)
        self.assertEqual(None, new_mc_3.specific_activity)
        self.assertEqual(None, new_mc_4.molecular_weight)
        self.assertEqual(None, new_mc_4.specific_activity)
        self.assertEqual(None, new_mc_4.wild_type)
        
        
if __name__ == '__main__':
    unittest.main()
    