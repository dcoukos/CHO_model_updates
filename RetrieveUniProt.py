#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 16:03:28 2018

@author: dimitricoukos
"""
import re
from urllib.error import HTTPError
from Bio.KEGG import REST
from Bio.SeqUtils import molecular_weight
from DataTreatment import openJson, write

mammals = ['HSA', 'PTR', 'PPS', 'GGO', 'PON', 'NLE', 'MCC', 'MCF', 'RRO',
           'RBB', 'CJC', 'SBQ', 'MMU', 'RNO', 'CGE', 'NGI', 'HGL', 'OCU',
           'TUP', 'CFA', 'AML', 'UMR', 'ORO', 'FCA', 'PTG', 'AJU', 'BTA',
           'BOM', 'BIU', 'PHD', 'CHX', 'OAS', 'SSC', 'CFR', 'CDK', 'LVE',
           'OOR', 'ECB', 'EPZ', 'EAI', 'MYB', 'MYD', 'HAI', 'RSS', 'LAV',
           'TMU', 'MDO', 'SHR', 'OAA']

animals = ['HSA', 'PTR', 'PPS', 'GGO', 'PON', 'NLE', 'MCC', 'MCF', 'RRO',
           'RBB', 'CJC', 'SBQ', 'MMU', 'RNO', 'CGE', 'NGI', 'HGL', 'OCU',
           'TUP', 'CFA', 'AML', 'UMR', 'ORO', 'FCA', 'PTG', 'AJU', 'BTA',
           'BOM', 'BIU', 'PHD', 'CHX', 'OAS', 'SSC', 'CFR', 'CDK', 'LVE',
           'OOR', 'ECB', 'EPZ', 'EAI', 'MYB', 'MYD', 'HAI', 'RSS', 'LAV',
           'TMU', 'MDO', 'SHR', 'OAA', 'GGA', 'MGP', 'CJO', 'TGU', 'GFR',
           'FAB', 'PHI', 'CCW', 'FPG', 'FCH', 'CLV', 'EGZ', 'AAM', 'ASN',
           'AMJ', 'PSS', 'CMY', 'SEA', 'ACS', 'PVT', 'PBI', 'GJA', 'XLA',
           'XTR', 'NPR', 'DRE', 'SRX', 'SGH', 'IPU', 'TRU', 'TNG', 'LCO',
           'NCC', 'MZE', 'OLA', 'XMA', 'NFU', 'LCF', 'HCQ', 'ELS', 'SFM',
           'LCM', 'CMK']


def returnBestAddress(gene_list, loop):
    '''

    Function in split into blocks to prevent unnecessary function calls.'''
    if loop == 'best':
        gene_dict = {}
        for entry in gene_list:
            split_seq = list(filter(None, re.split("[, \-:()]+", entry)))
            f_seq = []
            if split_seq == []:
                continue
            for entry in split_seq:
                f_seq.append(entry.strip())
            gene_dict[f_seq[0]] = list(filter(lambda x: x.isdigit(), f_seq))
        if 'CGE' in gene_dict:
            for index, entry in enumerate(gene_dict['CGE']):
                gene_dict['CGE'][index] = 'cge: ' + entry
            return gene_dict['CGE']
        elif 'MMU' in gene_dict:
            for index, entry in enumerate(gene_dict['CGE']):
                gene_dict['MMU'][index] = 'mmu: ' + entry
            return gene_dict['MMU']
        elif 'RNO' in gene_dict:
            for index, entry in enumerate(gene_dict['RNO']):
                gene_dict['RNO'][index] = 'rno: ' + entry
            return gene_dict['RNO']
        elif 'HSA' in gene_dict:
            for index, entry in enumerate(gene_dict['HSA']):
                gene_dict['HSA'][index] = 'hsa: ' + entry
            return gene_dict['HSA']
        else:
            loop = 'mammals'
    if loop == 'mammals':
        mammal_match = set(gene_dict.keys()).intersection(mammals)
        if bool(mammal_match):
            return mammal_match
        else:
            loop = 'vertebrates'
    if loop == 'vertebrates':
        animal_match = set(gene_dict.keys()).intersection(animals)
        if bool(animal_match):
            return animal_match
        else:
            loop = 'csm'  # Stands for "common simple models"
    if loop == 'csm':
        if 'DME' in gene_dict:
            for index, entry in enumerate(gene_dict['CGE']):
                gene_dict['DME'][index] = 'dme: ' + entry
            return gene_dict['DME']
        elif 'SCE' in gene_dict:
            for index, entry in enumerate(gene_dict['SCE']):
                gene_dict['SCE'][index] = 'sce: ' + entry
            return gene_dict['SCE']
        elif 'ECO' in gene_dict:
            for index, entry in enumerate(gene_dict['ECO']):
                gene_dict['ECO'][index] = 'eco: ' + entry
            return gene_dict['ECO']


def loopHandler(mol_weights, bigg_id, genes, loop):
    searching = True
    while searching:
        best = returnBestAddress(genes, loop)
        if not best:
            if loop == 'best':
                loop = 'mammals'
                break
            if loop == 'mammals':
                loop = 'vertebrates'
                break
            if loop == 'vertebrates':
                loop = 'csm'
                break
            if loop == 'csm':
                searching = False
                return None
        searching = False
    mol_weights[bigg_id]['genes'] = best
    mol_weights[bigg_id]['sequences'] = []
    mol_weights[bigg_id]['molecular_weights'] = []
    mol_weights[bigg_id]['uniprot_ids'] = []
    if loop == 'best' or loop == 'csm':
        for address in best:
            try:
                fillData(mol_weights, bigg_id, address, genes, loop)
            except HTTPError as err:
                if err.code == 404:
                    pass
    else:
        for gene in best:
            for address in best[gene]:
                try:
                    fillData(mol_weights, bigg_id, address, genes, loop)
                except HTTPError as err:
                    if err.code == 404:
                        pass


def fillData(mol_weights, bigg_id, address, genes, loop):
    text = REST.kegg_get(address).read()
    start_index = text.index('AASEQ')
    end_index = text.index('NTSEQ')
    raw_code = text[start_index:end_index].split('\n', 1)[1]
    code = raw_code.split('\n')
    sequence = ''
    for piece in code:
        sequence = sequence + piece.strip()
    mol_weights[bigg_id]['sequences'].append(sequence)
    ambigous_count = sequence.count('X')
    mod_seqeunce = sequence.replace('X', '')
    weight = molecular_weight(mod_seqeunce, seq_type='protein')
    weight = weight + 110*ambigous_count  # Estimate
    mol_weights[bigg_id]['molecular_weights'].append(weight)
    uni_index = text.find('UniProt')
    if uni_index != -1:
        mol_weights[bigg_id]['uniprot_ids'].append(text[uni_index+9:
                                                        uni_index+15])


if __name__ == '__main__':
    mol_weights = {}
    brenda_parameters = openJson('JSONs/brenda_parameters.json')
    for bigg_id in brenda_parameters:
        # TODO: restructure area to feed loop arguments.
        # TODO: implement dict data structure.
        # TODO: implement control methods for incomplete data.
        mol_weights[bigg_id] = {}
        print('Currently processing BiGG id: ' + bigg_id)
        try:
            ec_number = brenda_parameters[bigg_id][0]
            mol_weights[bigg_id]['ec_number'] = ec_number
            print(ec_number)
            text = REST.kegg_get('ec:' + ec_number).read()
            try:
                start_index = text.index('GENES') + 5
                end_index = text.index('DBLINKS')
            except ValueError:
                continue
            gene_string = text[start_index: end_index]
            genes = gene_string.split('\n')
            for index, gene in enumerate(genes):
                genes[index] = gene.strip()
        # CHANGED: no need to select closest organism here. Done in
        # returnBestAddress.
        except HTTPError as err:
            if err.code == 404:
                print('Excepted: Error 1')
                continue
            else:
                raise
        except ValueError:
            start_index = text.index('Now EC ')
            print('New EC')
            new_ec = text[start_index+6: start_index+20].split(',')[0]
            mol_weights[bigg_id]['ec_number'] = new_ec
            try:
                text = REST.kegg_get('ec:' + new_ec).read()
                try:
                    start_index = text.index('GENES') + 5
                    end_index = text.index('DBLINKS')
                except ValueError:
                    continue
                gene_string = text[start_index: end_index]
                genes = gene_string.split('\n')
                for index, gene in enumerate(genes):
                    genes[index] = gene.strip()
            except:
                print(new_ec)
                raise
        loop = 'best'
        searching = True
        while searching:
            try:
                loopHandler(mol_weights, bigg_id, genes, loop)
                searching = False
            except HTTPError as err:
                if err.code == 404 and loop == 'csm':
                    searching = False
            except TypeError as err:
                if loop == 'best':
                    loop = 'mammals'
                if loop == 'mammals':
                    loop = 'vertebrates'
                if loop == 'vertebrates':
                    loop = 'csm'
                if loop == 'csm':
                    searching = False

write('JSONs/molecular_weights.json', mol_weights)
