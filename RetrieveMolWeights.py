#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 16:03:28 2018

@author: dimitricoukos

Test: in command line:
    python RetrieveUniProt.py 'Unit Tests/sample_brenda_parameters.json'
"""
import sys
import cobra_services as CS
from multiprocessing import Pool
from urllib.error import HTTPError
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


def returnBestAddress(genes, loop):
    """Searches for available genes matching kegg enzyme entry.

    This function searches 'sequentially'. It returns the best available model
    organism genes. Organisms phylogenetically closer to Cricetulus griseus are
    preferred, but they are chosen by approximation. A detailed study of the
    phylogenetic tree has not been done for this project. Hopefully going
    sequentially increases both readability and efficiency.

    Parameters
    ----------
    genes : dict
        key: value pair is organism: address
    loop : string
        Indicates the highest potential group of matching organisms to search
        in.

    Returns
    -------
    dict
        key: kegg organism code. value: gene addresses for enzyme and organism

    """
    if loop == 'best':
        if 'CGE' in genes:
            return genes['CGE']
        elif 'MMU' in genes:
            return genes['MMU']
        elif 'RNO' in genes:
            return genes['RNO']
        elif 'HSA' in genes:
            return genes['HSA']
        else:
            loop = 'mammals'
    if loop == 'mammals':
        mammal_match = set(genes.keys()).intersection(mammals)
        if bool(mammal_match):
            return mammal_match
        else:
            loop = 'vertebrates'
    if loop == 'vertebrates':
        animal_match = set(genes.keys()).intersection(animals)
        if bool(animal_match):
            return animal_match
        else:
            loop = 'csm'  # Stands for "common simple models"
    if loop == 'csm':
        if 'DME' in genes:
            return genes['DME']
        elif 'SCE' in genes:
            return genes['SCE']
        elif 'ECO' in genes:
            return genes['ECO']


def loopHandler(mol_weights, ec_number, genes, loop):
    """Calls the correct loop of returnBestAddress based on best potential genes
        matches.

    Parameters
    ----------
    mol_weights : list
        empty list. will contain estimated molecular weights of enzymes.
    ec_number : string
    genes : list
        Addresses of genes corresponding to ec number.
    loop : string
    """
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
    mol_weights[ec_number]['weights'] = []
    mol_weights[ec_number]['uniprot_ids'] = []
    if loop == 'best' or loop == 'csm':
        for address in best:
            organism = best          # for readability
            try:
                fillData(mol_weights, ec_number, organism, address)
            except HTTPError as err:
                if err.code == 404:
                    pass
    else:
        for gene in best:
            for address in best[gene]:
                organism = best[gene]          # for readability
                try:
                    fillData(mol_weights, ec_number, organism, address)
                except HTTPError as err:
                    if err.code == 404:
                        pass


def fillData(mol_weights, ec_number, organism, address):
    """Searches kegg for enzyme uniprot id and AA sequence.

    Parameters
    ----------
    mol_weights : dict
        object containing all information collected by program.
    ec_number : string
        enzyme classification number used to organize data.
    address : string
        gene address for sequence lookup.
    """
    mol_weights[ec_number]['genes'].append(organism.lower() + ':' + address)
    sequence = CS.kegggene_to_sequence(organism, address)
    weight = CS.sequence_weight(sequence)
    mol_weights[ec_number]['weights'].append(weight)
    uniprot = CS.kegggene_to_uniprotid(organism, address)
    if uniprot:
        mol_weights[ec_number]['uniprot_ids'].uniprot


def mainSubprocess(bigg_ids, del_ec):
    """Main function called by each multiprocessing.process.

    Parameters
    ----------
    bigg_ids : dict
        key: ec_number. value: corresponding bigg ids.
    del_ec : list
        empty list which is appended to here containing depicrated ec numbers

    Returns
    -------
    dict
        key: ec number. value: all collected data in program by this process.
    """

    try:
        mol_weights = {}
        for ec_number in bigg_ids:  # WARNING: joblib may require list
            mol_weights[ec_number] = {}
            print('Currently processing BiGG id: ' + ec_number)
            mol_weights[ec_number]['bigg ids'] = bigg_ids[ec_number]
            try:
                genes = CS.ecnumber_to_genes(ec_number)
            except HTTPError as err:
                if err.code == 404:
                    print('Excepted: No entry for ec number: '+ec_number)
                    continue
                else:
                    raise
            if genes:
                loop = 'best'
                searching = True
                while searching:
                    try:
                        loopHandler(mol_weights, ec_number, genes, loop)
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
    finally:
        return mol_weights


if __name__ == '__main__':
    sub_dict_1 = {}
    sub_dict_2 = {}
    sub_dict_3 = {}
    sub_dict_4 = {}

    mol_weights = {}
    if len(sys.argv) == 1:
        brenda_parameters = openJson('JSONs/brenda_parameters.json')
    else:
        brenda_parameters = openJson(sys.argv[1])
    simplified_brenda = {}
    for bigg_id in brenda_parameters:
        simplified_brenda[bigg_id] = brenda_parameters[bigg_id][0]
    optimized_bigg = {}
    for k, v in simplified_brenda.items():
        optimized_bigg[v] = optimized_bigg.get(v, [])
        optimized_bigg[v].append(k)
    counter = 0
    for ec_number in optimized_bigg:
        if counter % 4 == 0:
            sub_dict_1[ec_number] = optimized_bigg[ec_number]
        if counter % 4 == 1:
            sub_dict_2[ec_number] = optimized_bigg[ec_number]
        if counter % 4 == 2:
            sub_dict_3[ec_number] = optimized_bigg[ec_number]
        if counter % 4 == 3:
            sub_dict_4[ec_number] = optimized_bigg[ec_number]
        counter = counter + 1
    try:
        with Pool(processes=4) as pool:
            del_ec1 = []
            del_ec2 = []
            del_ec3 = []
            del_ec4 = []
            mw_1 = pool.apply_async(mainSubprocess, (sub_dict_1, del_ec1,))
            mw_2 = pool.apply_async(mainSubprocess, (sub_dict_2, del_ec2,))
            mw_3 = pool.apply_async(mainSubprocess, (sub_dict_3, del_ec3,))
            mw_4 = pool.apply_async(mainSubprocess, (sub_dict_4, del_ec4,))
            pool.close()
            pool.join()
            for ec in del_ec1:
                mw_1.pop(ec, None)
            for ec in del_ec2:
                mw_2.pop(ec, None)
            for ec in del_ec3:
                mw_3.pop(ec, None)
            for ec in del_ec4:
                mw_4.pop(ec, None)
    finally:
        mol_weights.update(mw_1.get())
        mol_weights.update(mw_2.get())
        mol_weights.update(mw_3.get())
        mol_weights.update(mw_4.get())
        if len(sys.argv) > 1:
            write('Unit Tests/multiprocessing_sub_output1.json', mw_1.get())
            write('Unit Tests/multiprocessing_sub_output3.json', mw_3.get())
        mol_weights_to_write = {}
        for ec_number in mol_weights:
            for bigg_id in mol_weights[ec_number]['bigg ids']:
                mol_weights_to_write[bigg_id] = {}
                mol_weights_to_write[bigg_id]['ec_number'] = ec_number
                mol_weights_to_write[bigg_id].update(mol_weights[ec_number])
        write('JSONs/molecular_weights.json', mol_weights_to_write)
