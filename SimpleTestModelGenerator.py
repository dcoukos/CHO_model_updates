#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 20:51:19 2018

@author: dimitricoukos
"""
import cobra
import cobra.test
from cobra import Reaction, Metabolite
model = cobra.test.create_test_model('textbook')

reaction = Reaction('CSND')
reaction2 = Reaction('DHPM1')
reaction.name = 'Test reaction 1'
reaction.subsystem = 'cytoplasm'

cyt = Metabolite(
    'cyt',
    name='cyt',
    compartment='c')
fluorocyt = Metabolite(
    "5-fluorocyt",
    name='5-fluorocyt',
    compartment='c'
)
unmatchable = Metabolite(
    'Unmatchable_metabolite',
    name='Unmatchable_metabolite',
    compartment='c'
)

DL_p = Metabolite(
    'DL-p'
)

ethylhydantoin = Metabolite(
    'ethylhydantoin'
)


reaction.add_metabolites({
    cyt: -1.0,
    fluorocyt: 1.0,
    unmatchable: 1.0
})

reaction2.add_metabolites({
    DL_p: 1.0,
    ethylhydantoin: -1.0
})

model.add_reaction(reaction)
model.add_reaction(reaction2)

cobra.io.save_json_model(model, 'Unit Tests/simple_test_model.json')
