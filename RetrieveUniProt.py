#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 16:03:28 2018

@author: dimitricoukos
"""

from urllib.error import HTTPError
from Bio.KEGG import REST, Enzyme
from Bio.SeqUtils import molecular_weight
from DataTreatment import openJson, write

def returnBestAddress(gene_list):
    gene_dict = {}
    for entry in gene_list:
        org, gene = entry.split(':')
        gene_dict[org] = gene.split('(')[0].strip()
    if 'CGE' in gene_dict:
        return 'cge: ' + gene_dict['CGE']
    elif 'MMU' in gene_dict:
        return 'mmu: '+ gene_dict['MMU']
    elif 'RNO' in gene_dict:
        return 'rno: '+ gene_dict['RNO']
    elif 'HSA' in gene_dict:
        return 'hsa: '+ gene_dict['HSA']
    mammal_match = set(gene_list).intersection(mammals)
    if bool(mammal_match):
        return mammal_match.pop()
    elif 'DME' in gene_dict:
        return 'dme: '+ gene_dict['DME']
    animal_match = set(gene_list).intersection(animals)
    if bool(animal_match):
        return animal_match.pop()
    elif 'SCE' in gene_dict:
        return 'sce: '+ gene_dict['SCE']
    elif 'ECO' in gene_dict:
        return 'eco: '+ gene_dict['ECO']

def returnAddress(gene):
    org, address = gene.split(':')
    return org.lower() + ': ' + address.split('(')[0].strip()

mammals = ['HSA', 'PTR', 'PPS', 'GGO', 'PON', 'NLE', 'MCC', 'MCF', 'RRO',
    'RBB', 'CJC', 'SBQ', 'MMU', 'RNO', 'CGE', 'NGI', 'HGL', 'OCU', 'TUP',
    'CFA', 'AML', 'UMR', 'ORO', 'FCA', 'PTG', 'AJU', 'BTA', 'BOM', 'BIU',
    'PHD', 'CHX', 'OAS', 'SSC', 'CFR', 'CDK', 'LVE', 'OOR', 'ECB', 'EPZ',
    'EAI', 'MYB', 'MYD', 'HAI', 'RSS', 'LAV', 'TMU', 'MDO', 'SHR', 'OAA']

animals = ['HSA', 'PTR', 'PPS', 'GGO', 'PON', 'NLE', 'MCC', 'MCF', 'RRO',
    'RBB', 'CJC', 'SBQ', 'MMU', 'RNO', 'CGE', 'NGI', 'HGL', 'OCU', 'TUP',
    'CFA', 'AML', 'UMR', 'ORO', 'FCA', 'PTG', 'AJU', 'BTA', 'BOM', 'BIU',
    'PHD', 'CHX', 'OAS', 'SSC', 'CFR', 'CDK', 'LVE', 'OOR', 'ECB', 'EPZ',
    'EAI', 'MYB', 'MYD', 'HAI', 'RSS', 'LAV', 'TMU', 'MDO', 'SHR', 'OAA',
    'GGA', 'MGP', 'CJO', 'TGU', 'GFR', 'FAB', 'PHI', 'CCW', 'FPG', 'FCH',
    'CLV', 'EGZ', 'AAM', 'ASN', 'AMJ', 'PSS', 'CMY', 'SEA', 'ACS', 'PVT',
    'PBI', 'GJA', 'XLA', 'XTR', 'NPR', 'DRE', 'SRX', 'SGH', 'IPU', 'TRU',
    'TNG', 'LCO', 'NCC', 'MZE', 'OLA', 'XMA', 'NFU', 'LCF', 'HCQ', 'ELS',
    'SFM', 'LCM', 'CMK']

if __name__ == '__main__':
    mol_weights = {}
    brenda_parameters = openJson('JSONs/brenda_parameters.json')
    for bigg_id in brenda_parameters:
        #IN PROGRESS:
        print('Currently processing BiGG id: ' + bigg_id)
        try:
            ec_number = brenda_parameters[bigg_id][0]
            print(ec_number)
            text = REST.kegg_get('ec:' + ec_number).read()
            start_index = text.index('GENES') + 5
            end_index = text.index('DBLINKS')
            gene_string = text[start_index: end_index]
            genes = gene_string.split('\n')
            for index, gene in enumerate(genes):
                genes[index] = gene.strip()
            close_genes = list(filter(lambda x: x[0:3] in (animals + ['SCE',
                'ECO']) and x[3] == ':', genes ))
            print(close_genes)
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
            try:
                text = REST.kegg_get('ec:' + new_ec).read()
                start_index = text.index('GENES') + 5
                end_index = text.index('DBLINKS')
                gene_string = text[start_index: end_index]
                genes = gene_string.split('\n')
                for gene in genes:
                    gene = gene.strip()
                close_genes = list(filter(lambda x: x[0:3] in (animals + ['SCE',
                    'ECO']) and x[3] == ':', genes ))
                print(close_genes)
            except:
                print(new_ec)
                raise
        try:
            text = REST.kegg_get(returnBestAddress(close_genes)).read()
            start_index = text.index('AASEQ')
            end_index = text.index('NTSEQ')
            raw_code = text[start_index:end_index].split('\n',1)[1]
            code = raw_code.split('\n')
            sequence = ''
            for piece in code:
                sequence = sequence + piece.strip()
            mol_weights[bigg_id] = molecular_weight(sequence, seq_type='protein')
            print('Molecular weight for ' + bigg_id + ' added to uniprot_ids')
        except HTTPError or ValueError as err:
            if err.code == 404:
                print('Excepted: Error 2')
                for gene in close_genes:
                    try:
                        text = REST.kegg_get(returnAddress(gene)).read()
                        start_index = text.index('AASEQ')
                        end_index = text.index('NTSEQ')
                        raw_code = text[start_index:end_index].split('\n',1)[1]
                        code = raw_code.split('\n')
                        sequence = ''
                        for piece in code:
                            sequence = sequence + piece.strip()
                        mol_weights[bigg_id] = molecular_weight(sequence,
                            seq_type='protein')
                        print('Molecular weight for ' + bigg_id + ' added to uniprot_ids')
                        break
                    except HTTPError as err:
                        print('Excepted: Error 3')
                    except ValueError:
                        print('Excepted: Error 4')
        except ValueError:
            print(text)
            raise

    write('JSONs/molecular_weights.json', mol_weights)
