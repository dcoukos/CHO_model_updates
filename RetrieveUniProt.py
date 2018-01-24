#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 16:03:28 2018

@author: dimitricoukos
"""

from urllib.error import HTTPError
from Bio.KEGG import REST, Enzyme
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
    elif 'DME' in gene_dict:
        return 'dme: '+ gene_dict['DME']
    elif 'SCE' in gene_dict:
        return 'sce: '+ gene_dict['SCE']
    elif 'ECO' in gene_dict:
        return 'eco: '+ gene_dict['ECO']

def returnAddress(gene):
    org, address = gene.split(':')
    return lower(org) + ': ' + address.split('(')[0].strip()

if __name__ == '__main__':
    uniprot_ids = {}
    brenda_parameters = openJson('JSONs/brenda_parameters.json')
    for bigg_id in brenda_parameters:
        #IN PROGRESS:
        print('Currently processing BiGG id: ' + bigg_id)
        try:
            ec_number = brenda_parameters[bigg_id][0]
            text = REST.kegg_get('ec:' + ec_number).read()
            passed_genes = False
            passed_end = False
            genes = []
            for line in text.rstrip().split('\n'):
                section = line[:12].strip()
                if section == 'GENES':
                    gene = line.split(' ')[-2:]
                    line = gene[0]+' '+gene[1]
                    passed_genes = True
                if section == 'DBLINKS':
                    passed_end == True
                if passed_genes == True and passed_end == False:
                    genes.append(line.strip())
            close_genes = list(filter(lambda x: x.split(':') in [
                'CGE',
                'MMU',
                'RNO',
                'HSA',
                'DME',
                'SCE',
                'ECO'], genes ))
        except HTTPError as err:
            if err.code == 404:
                uniprot_ids[bigg_id] = ''
            else:
                raise
        try:
            text = REST.kegg_get(returnBestAddress(close_genes))
            index = text.index('UniProt')
            uniprot_id = text[index+9:index+15]
            uniprot_ids[bigg_id] = uniprot_id
            print('id ' + bigg_id + ' added to uniprot_ids')
        except HTTPError or ValueError as err:
            if err.code == 404:
                for gene in close_genes:
                    try:
                        text = REST.kegg_get(returnBestAddress(close_genes))
                        index = text.index('UniProt')
                        uniprot_id = text[index+9:index+15]
                        uniprot_ids[bigg_id] = uniprot_id
                        print('id ' + bigg_id + ' added to uniprot_ids')
                        break
                    except HTTPError as err:
                        continue
                    except ValueError:
                        continue

    write('JSONs/UniProt_Ids.json', uniprot_ids)
