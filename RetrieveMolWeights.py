#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 16:03:28 2018

@author: dimitricoukos

Test: in command line:
    python RetrieveUniProt.py 'Unit Tests/sample_brenda_parameters.json'
"""
# Sandbox

import cobra
import argparse
from statistics import mean
from Bio.KEGG import REST
from Bio.SeqUtils import molecular_weight
import cobra_services as CS
from multiprocessing import Pool
from sympy.logic.boolalg import simplify_logic
from urllib.error import HTTPError
from DataTreatment import write
from progress.bar import Bar


def weightsFromModel(ids, model, process):
    '''
        Retrieve molecular weight based on gene reaction rules from model.
        Once the logic is parsed, kegg can be queried almost directly.

        ids: dictionary with Bigg IDs as keys, with None values
        model: cobra model.
    '''
    add_not_found = []
    rules = {}
    sequences = {}
    mol_weights = {}
    for id in ids:
        rules[id] = model.reactions.get_by_id(id).gene_reaction_rule
    char_to_int = {'a': '0', 'b': '1', 'c': '2', 'd': '3', 'e': '4',
                   'f': '5', 'g': '6', 'h': '7', 'i': '8', 'j': '9'}
    if process == 1:
        pass
        bar = Bar('Retrieving Weights: ', max=len(ids))
    for id in rules:
        if process == 1:
            bar.next()
        if rules[id] == '':
            continue
        expr = rules[id]
        weights = []
        split_rules = []
        expr = expr.replace('and', '&')
        expr = expr.replace('or', '|')
        let = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
        # print('*', end='')
        try:
            # print('!', end='')
            if 's' in expr:
                expr = expr.replace('s', '')
            for char in expr:
                if char not in ['(', ')', '&', '|', ' ']:
                    expr = expr.replace(char, let[int(char)])
            if '(' in expr or '&' in expr:
                ors = len(expr.split('|'))
                ands = len(expr.split('&'))
                if ors + ands > 10:
                    continue
                else:
                    expr = simplify_logic(expr, form='dnf')
        except:
            # print('@', end='')
            print('Expression throwing error:',  expr, ' || ID: ', id)
            if expr == '':
                continue
            else:
                raise
        clause_weights = []
        for clause in str(expr).replace('(', '').replace(')',
                                                         '').split(' | '):
            comb_weight = 0
            split_rules = clause.split(' & ')
            for entry in split_rules:
                # print('#', end='')

                for char in entry:
                    if char not in ['(', ')', '&', '|', ' ']:
                        entry = entry.replace(char, char_to_int[char])
                try:
                    # print('$', end='')
                    if entry in sequences:
                        sequence = sequences[entry]
                    else:
                        sequence = kegggene_to_sequence('cge', entry)
                        sequences[entry] = sequence
                        # TODO: do this sort of bug check with all 4 processes.
                        # if process == 1:
                        # print('\t|| ', id)
                        # print('request')
                except HTTPError:
                    # print('%', end='')
                    mess = 'Id: ' + id + ' || Address: cge:' + entry
                    add_not_found.append(mess)
                    continue
                comb_weight += sequence_weight(sequence)
                # print('^', end='')
            clause_weights.append(comb_weight)
        mol_weights[id] = mean(clause_weights)
    print('Process ', process, 'done finding minimum weights')
    return mol_weights


def sequence_weight(sequence):
    "Return weight in Daltons"
    ambigous_count = sequence.count('X')
    mod_seqeunce = sequence.replace('X', '')
    weight = molecular_weight(mod_seqeunce, seq_type='protein')
    weight = weight + 110*ambigous_count  # Estimate
    return weight


def minimum(weights):
    mini = 0
    for weight in weights:
        if weight > mini:
            mini = weight
    return mini

# TODO: fix cobra_services so that these work from the package


def kegggene_to_sequence(organism, kegggene):
    """
    Get the sequence of a gene
    """

    text = REST.kegg_get(organism.lower() + ':' + kegggene).read()

    start_index = text.index('AASEQ')
    end_index = text.index('NTSEQ')

    raw_code = text[start_index: end_index].split('\n', 1)[1]
    code = raw_code.split('\n')
    sequence = ''
    for piece in code:
        sequence = sequence + piece.strip()

    return sequence


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='This is a script to model'
        'CHO cell growth using updated specific activities and cobra.')
    parser.add_argument('-s', '--single', help='Run script as single process',
                        action='store_true')
    args = parser.parse_args()
    mol_weights = {}
    if args.single:
        model = cobra.io.read_sbml_model('iCHOv1_K1_final.xml')
        ids = {}
        for reaction in model.reactions:
            ids[reaction.id] = None
        mol_weights = weightsFromModel(ids, model, 1)
    else:
        sub_dict_1 = {}
        sub_dict_2 = {}
        sub_dict_3 = {}
        sub_dict_4 = {}

        counter = 0
        model = cobra.io.read_sbml_model('iCHOv1_K1_final.xml')

        for reaction in model.reactions:
            if counter % 4 == 0:
                sub_dict_1[reaction.id] = None
            if counter % 4 == 1:
                sub_dict_2[reaction.id] = None
            if counter % 4 == 2:
                sub_dict_3[reaction.id] = None
            if counter % 4 == 3:
                sub_dict_4[reaction.id] = None
            counter = counter + 1

        try:
            with Pool(processes=4) as pool:
                mw1 = pool.apply_async(weightsFromModel, (sub_dict_1,
                                                          model, 1,))
                mw2 = pool.apply_async(weightsFromModel, (sub_dict_2,
                                                          model, 2,))
                mw3 = pool.apply_async(weightsFromModel, (sub_dict_3,
                                                          model, 3,))
                mw4 = pool.apply_async(weightsFromModel, (sub_dict_4,
                                                          model, 4,))
                pool.close()
                pool.join()
        except:
            raise
        mol_weights.update(mw1.get())
        mol_weights.update(mw2.get())
        mol_weights.update(mw3.get())
        mol_weights.update(mw4.get())

    write('JSONs/molecular_weights_ave.json', mol_weights)
