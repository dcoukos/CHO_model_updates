#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 23:45:56 2018

@author: dimitricoukos
"""
import unittest
import json
import DataTreatment 

class sampleData(unittest.TestCase):
    
 
    def test_file_correction(self):
        '''correct JSON should be able to reverse KEGG codes and BRENDA names
            that have been reversed.
        '''
        brenda_keggs = DataTreatment.correctJSON('tests/incorrect_json.json')
        with open('tests/correct_json.json') as infile:
            correct_file = json.load(infile, indent=4)
        self.assertEqual(brenda_keggs, correct_file)
        
if __name__ == '__main__':
    unittest.main()