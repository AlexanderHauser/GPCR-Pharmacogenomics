"""
Created in Jan 2017

@author: Alexander Hauser <alexshauser@gmail.com>

extracts residue contacts between structures based on:
Arpeggio: http://bleoberis.bioc.cam.ac.uk/arpeggioweb/
interaction types are described here:
https://bitbucket.org/harryjubb/arpeggio/src/master/README.md?fileviewer=file-view-default
"""

import os, sys, glob
import pandas as pd
import subprocess
import numpy as np
import random
import requests

def structureInfo(PDB):

    url = "http://test.gpcrdb.org/services/structure/" + PDB + "/"
    structure = requests.get(url).json()

    return structure

def get_protein_residues(EntryName):

    url = "http://test.gpcrdb.org/services/residues/" + EntryName
    res_data = requests.get(url).json()

    return res_data

def GPCRdb_numbering_translate(df):
    ## translate classB to classA numbering scheme
    # class B 1x50 = class A 1x46
    table = {'TM1':-4, 'TM2':-7, 'TM3': -4, 'TM4': 0, 'TM5': 4, 'TM6': -5, 'TM7': -4, 'H8': 0, 'ICL1': 0, 'ICL2': 1, 'ICL3':0}
    if df['Segment'] in table and str(df['GPCRdb']) != 'nan':
        print df['GPCRdb_short']
        return str(df['GPCRdb_short'].split('x')[0]) + 'x' + str(int(df['GPCRdb_short'].split('x')[1]) + table[df['Segment']])

interactions = pd.DataFrame(columns=['EntryName','Residue','GPCRdb','GPCRdb_short','Segment','PDB'])
pdbs = ['5UZ7','5VAI'] # CLASS B
pdbs = ['4ZWJ','5DGY'] #'4PXF' # ARRESTIN
pdbs = ['3SN6','5G53'] # CLASS A GPCR-Gprotein COMPLEX

## Do arepeggio calculations first by using local install of arpeggio
# changed default 5.0 A to 4.5 A distance
## e.g.:
## arpeggio 5G53.pdb -s /A//

for pdb in pdbs:
    contactfile = glob.glob("../data/interactions_Gprotein/" + pdb + "/*.contacts")[0]
    # contactfile = glob.glob("../data/interactions_arrestin/" + pdb + "/*.contacts")[0]

    contacts = pd.read_csv(contactfile, sep='\s+', header=None)

    info = structureInfo(pdb)

    # Extract Res and chain info
    contacts['chain1'] = contacts[0].str.split('/', expand=True)[0]
    contacts['res1'] = contacts[0].str.split('/', expand=True)[1]
    contacts['chain2'] = contacts[1].str.split('/', expand=True)[0]
    contacts['res2'] = contacts[1].str.split('/', expand=True)[1]

    # Look for cross chain interactions:
    interchain = contacts[(contacts.chain1 != contacts.chain2) & (contacts.chain1 != 'P') & (contacts.chain2 != 'P')]

    ## take receptor residues only
    ## chain R in Arpeggio for G protein interactions
    ## chain A in Arpeggio for Arrestin interactions
    ## specific interactions defined here: https://bitbucket.org/harryjubb/arpeggio/src/master/README.md?fileviewer=file-view-default

    if 'arrestin' in contactfile:
        receptor_residues = list(set(list(interchain[(interchain.chain1=='A')].res1.unique()) + list(interchain[(interchain.chain2=='A')].res2.unique())))
        entryName = info['protein']
    else:
        entryName = info['protein']
        receptor_residues = list(set(list(interchain[(interchain.chain1=='R')].res1.unique()) + list(interchain[(interchain.chain2=='R')].res2.unique())))
        if not receptor_residues:
            receptor_residues = list(set(list(interchain[(interchain.chain1=='A')].res1.unique()) + list(interchain[(interchain.chain2=='A')].res2.unique())))

    # map residues on receptor --> gives arrestin specific residue positions
    protein_residue_info = get_protein_residues(entryName)

    for res in receptor_residues:
        pos = interactions.shape[0]+1

        gpcrdb = str(protein_residue_info[int(res)-1]['display_generic_number'])
        segment = str(protein_residue_info[int(res)-1]['protein_segment'])
        interactions.loc[pos,'Segment'] = segment

        if 'x' in gpcrdb:
            interactions.loc[pos,'GPCRdb'] = gpcrdb
            interactions.loc[pos,'GPCRdb_short'] = gpcrdb.split('.')[0] + 'x' + gpcrdb.split('x')[1]

        interactions.loc[pos,'EntryName'] = entryName
        interactions.loc[pos,'Residue'] = int(res)
        interactions.loc[pos,'PDB'] = pdb

## translate GPCRdb class A translation
# interactions['GPCRdb_classA'] = interactions.apply(GPCRdb_numbering_translate, axis=1)
interactions = interactions.sort_values(by='Residue').reset_index(drop=True)
interactions.to_csv("../data/classA_gprotein_contacts.csv")
