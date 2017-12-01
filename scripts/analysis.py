"""
Created in Dec 2016

@author: Alexander Hauser <alexshauser@gmail.com>

- GPCR Pharmacogenomics analysis
"""
import pandas as pd

import os, glob
from collections import OrderedDict

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec

import numpy as np
import re, math
from scipy import stats

import time
import datetime
from tqdm import *

now = datetime.datetime.now()
timestamp = str(now.year)[2:]+['0' if now.month < 10 else ''][0]+str(now.month)+['0' if now.day < 10 else ''][0]+str(now.day)

plt.style.use('ggplot')

### ============== FUNCTIONS ==============
def r2(x, y):
    return stats.pearsonr(x, y)[0] ** 2

def insert_string(string, insert, index):
    return string[:index] + insert + string[index:]

def get_significances(DF, group, value_column):

    from itertools import combinations
    from scipy.stats import ttest_ind
    from scipy.stats import f_oneway

    segment_groups = {}
    # for group in relativeRegionTable.groupby('Region')['RelativeSNP']:
    #     segment_groups[group[0]] = list(group[1])

    for group in DF.groupby(group)[value_column]:
        segment_groups[group[0]] = list(group[1])
    # independent t-tests:
    # Pairwise: Is a similar to b? Is a similar to c? Is b similar to c?
    # list1 = list(relativeRegionTable[relativeRegionTable['Region']=='TM1']['RelativeSNP'])
    # list2 = list(relativeRegionTable[relativeRegionTable['Region']=='ICL1']['RelativeSNP'])
    significance = pd.DataFrame()
    for segment1, segment2 in combinations(segment_groups.keys(), 2):
        t, p = ttest_ind(segment_groups[segment1], segment_groups[segment2])
        star = 'NaN'
        if p>=0.05:
            star = 'ns'
        if p<=0.05:
            star = '*'
        if p<=0.01:
            star = '**'
        if p<=0.001:
            star = '***'
        if p<=0.0001:
            star = '****'
        significance.loc[segment1,segment2]=star
        significance.loc[segment2,segment1]=star
    return significance

def cDNA_to_codon(cDNA, pos, Tcon):
    table={
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
    'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
    'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
    'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
    'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
    'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
    'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
    'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G'}

    positions = [304, 148, 311, 243, 150, 150, 139]
    Tcons = ['c.910A>C','c.442C>G','c.932T>A','c.727G>A','c.449G>A','c.448C>T','c.416G>T']
    for pos, Tcon in zip(positions,Tcons):
        Tpos = int(re.findall(r'\d+', Tcon)[0])
        Taa = Tcon[-1]
        WT = cDNA[(pos-1)*3:(pos-1)*3+3]
        mutantcDNA = cDNA[:Tpos-1] + Taa + cDNA[Tpos:]
        MT = mutantcDNA[(pos-1)*3:(pos-1)*3+3]
        print 'WT:', table[WT], WT
        print 'MT:', table[MT], MT
        print '--\n'
### ================= END =================

def mastertable():
    dfm = pd.read_csv("../data/ExAC_GPCRs_raw.csv")
    dfm = dfm.drop_duplicates(subset=["EntryName","Family","Ligandtype","Class","Uniprot","gene","GPCRdb","SequenceNumber","Segment","GPCRdbWT","NMaa","MutationType","Allele Count","Allele Frequency","Allele Number","Number of Homozygotes","sift_word","sift_score"])

    drugs_raw = pd.read_csv('../data/drug_data.csv')

    ## Exclude gnrr2_human as only a putative GPCR (5TM, not in GPCRdb)
    drugs_raw = drugs_raw[drugs_raw.EntryName!='gnrr2_human']

    ## Subset variant df on clinical targets
    drugTargets = drugs_raw[(drugs_raw.PMID.notnull()) & (drugs_raw['Status']=='approved')]['EntryName'].unique()
    trialTargets = np.array(list(set(drugs_raw['EntryName'].unique()) - set(drugs_raw[drugs_raw['Status']=='approved']['EntryName'].unique())))
    # drugs_raw[(drugs_raw.Status=='in trial') & (drugs_raw.ClinicalStatus.isin(ongoing)) & (drugs_raw.EntryName.isin(trialTargets))]['EntryName'].nunique()
    dfm_drugs = dfm[dfm.EntryName.isin(drugTargets)]

    ###  Create MasterTable
    ## LOAD DATASETS

    ## Mutations
    mutations = pd.read_csv("../data/Mutants_orthologs.csv")
    mutations['foldchange_abs'] = abs(mutations['foldchange'])
    mutations['NMaa'] = mutations['Mutantaa']
    mutations_pos = pd.DataFrame({'foldchangeMaxAbs' : mutations.groupby( ['SequenceNumber','EntryName','NMaa'] )['foldchange_abs'].max(),
    'foldchangeMaxPos': mutations.groupby( ['SequenceNumber','EntryName','NMaa'] )['foldchange'].max(),
    'foldchangeMaxNeg': mutations.groupby( ['SequenceNumber','EntryName','NMaa'] )['foldchange'].min()}).reset_index()

    ## Disease --> get only those with highest score for a given position!
    diseasemap = pd.read_csv("../data/Disgenet_disease_ass.csv")
    diseasemap = diseasemap[['SequenceNumber','EntryName','diseaseId','score','diseaseName']]
    diseasemap = diseasemap.sort_values(by=['SequenceNumber','EntryName','score'], ascending=False).groupby( ['SequenceNumber','EntryName'] ).first().reset_index()

    ## PTM
    ptms = pd.read_csv('../data/PTM_gpcrs.csv')
    ptms['Type'] = ptms['Type'].str.split("(", expand=True)[0].str.strip()
    ptms = ptms[['SequenceNumber','EntryName','Type']].drop_duplicates()
    ptms['PTMsite'] = 'yes'

    ## LB coorect for orthologs/different sequence numbers!!!
    # merge all human structures based on sequenceNumbers and all orthologs based on GPCRdb
    lb_structure = pd.read_csv('../data/GPCRdb_ligand_interactions.csv')
    lb_structure.rename(columns={'Residue': 'SequenceNumber'}, inplace=True)
    # To include also non-human crystals! --> get the right SequenceNumber to merge on humans!
    lb_structure_human = lb_structure[lb_structure.EntryName.str.contains('_human')][['EntryName', 'SequenceNumber','Interaction']].drop_duplicates()
    lb_structure_human['LB_structure_human'] = 'interacting'

    lb_structure_ortho = lb_structure[~lb_structure.EntryName.str.contains('_human')][['EntryName', 'GPCRdb','Interaction']].drop_duplicates()
    lb_structure_ortho['LB_structure_ortho'] = 'interacting'
    lb_structure_ortho['EntryName'] = lb_structure_ortho['EntryName'].str.split("_", expand=True)[0]+'_human'

    lb_fam = pd.read_csv('../data/generic_LB_site.csv')
    lb_fam['LB_fam'] = 'interacting'
    lb_fam['EntryName'] = lb_fam['EntryName'].str.split("_", expand=True)[0]+'_human'
    lb_fam['SequenceNumber'] = lb_fam['SequenceNumber'].astype(int)
    lb_fam = lb_fam[['LB_fam', 'EntryName', 'SequenceNumber']].drop_duplicates()

    ## Arrestin from arpeggio
    arrestin = pd.read_csv("../data/arrestin_contacts.csv")
    arrestin_positions = arrestin.dropna().GPCRdb.unique() # based on class!

    ## Gprotein from arpeggio
    # Gprot_positions = ['8x47', '8x48', '8x49', '8x51', '3x50', '3x53', '3x54', '3x55', '3x56', '5x61', '5x64', '5x65', '5x66', '5x67', '5x68', '5x69', '5x71', '5x72', '5x74', '5x75', '6x25', '6x26', '6x28', '6x29', '6x32', '6x33', '6x36', '6x37', '6x40','7x55', '7x56','34x50', '34x51', '34x52', '34x53', '34x54', '34x55', '34x56', '34x57'] # Class A + accessible

    Gprot_positions_df = pd.read_csv('../data/classA_gprotein_contacts.csv')
    Gprot_positions= Gprot_positions_df[Gprot_positions_df.GPCRdb_short.notnull()].GPCRdb_short.unique()

    Gprot_positions_B_df = pd.read_csv('../data/classB_gprotein_contacts.csv')
    Gprot_positions_B = Gprot_positions_B_df[Gprot_positions_B_df.GPCRdb_short.notnull()].GPCRdb_short.unique()

    ## =================== Big merge ================================
    dfm_drugs['GPCRdb_short'] = dfm_drugs['GPCRdb'].str.split('.', expand=True)[0]+'x'+dfm_drugs['GPCRdb'].str.split('x', expand=True)[1]

    ## Signalling interface class specific!
    def cross_class_position(GPCRdb_short, class_name):
        # translates cross_class positions to class A generic positions:
        GPCRdb_split = GPCRdb_short.split('x')

        if 'Class C' in class_name:
            translate_C = {'1': -4, '2':4, '3': -4, '4': 10, '5': 0, '6': -2, '7': 5, '8': 0, '12': 0, '34': 0, '45': 0}
            updated_position = int(GPCRdb_split[1]) + translate_C[GPCRdb_split[0]]
            return GPCRdb_split[0] + 'x' + str(updated_position)

        elif 'Class F' in class_name:
            translate_F = {'1': 3, '2': 1, '3': 0, '4': 0, '5': -4, '6': 1, '7': 0, '8': 0, '12': 0, '34': 0, '45': 0}
            updated_position = int(GPCRdb_split[1]) + translate_F[GPCRdb_split[0]]
            return GPCRdb_split[0] + 'x' + str(updated_position)

        elif 'Class B' in class_name:
            # class B 1x50 = class A 1x46
            translate_B = {'1': -4, '2': -7, '3': -4, '4': 0, '5': 4, '6': -5, '7': -4, '8': 0, '12': 0, '34': 0, '45': 0}
            updated_position = int(GPCRdb_split[1]) + translate_B[GPCRdb_split[0]]
            return GPCRdb_split[0] + 'x' + str(updated_position)

        else:
            # other classes for now and 'Class A' in class_name
            return GPCRdb_short

    def cross_class_position_CtoB(GPCRdb_short, class_name):
        # translates class C positions to class B generic positions:
        GPCRdb_split = GPCRdb_short.split('x')

        if 'Class C' in class_name:
            translate_C = {'1': 0, '2': 11, '3': 0, '4': 10, '5': -4, '6': 3, '7': 9, '8': 0, '12': 0, '34': 0, '45': 0}
            updated_position = int(GPCRdb_split[1]) + translate_C[GPCRdb_split[0]]
            return GPCRdb_split[0] + 'x' + str(updated_position)

    def cross_class_position_FtoB(GPCRdb_short, class_name):
        # translates class F positions to class B generic positions:
        GPCRdb_split = GPCRdb_short.split('x')

        if 'Class F' in class_name:
            translate_F = {'1': 7, '2': 8, '3': 4, '4': 0, '5': -8, '6': 6, '7': 4, '8': 0, '12': 0, '34': 0, '45': 0}
            updated_position = int(GPCRdb_split[1]) + translate_F[GPCRdb_split[0]]
            return GPCRdb_split[0] + 'x' + str(updated_position)

    def arrestin_interface(sub_df):
        class_name = sub_df.Class
        if 'x' in str(sub_df.GPCRdb_short):

            translated_generic_position = cross_class_position(sub_df.GPCRdb_short, class_name)

            if translated_generic_position in arrestin[arrestin.GPCRdb.notnull()].GPCRdb_short.unique():
                interaction = 'putative'
            else:
                interaction = None
        else:
            interaction = None

        return interaction

    def gprotein_interface(sub_df):
        ## different G protein interaction set for class A and B
        ## uniset of residue positions for Class C and F

        class_name = sub_df.Class
        if 'x' in str(sub_df.GPCRdb_short):

            translated_generic_position = cross_class_position(sub_df.GPCRdb_short, class_name)
            translated_generic_position_CtoB = cross_class_position_CtoB(sub_df.GPCRdb_short, class_name)
            translated_generic_position_FtoB = cross_class_position_FtoB(sub_df.GPCRdb_short, class_name)

            if 'Class A' in class_name and translated_generic_position in Gprot_positions:
                interaction = 'putative'
            elif 'Class B' in class_name and translated_generic_position in Gprot_positions_B:
                interaction = 'putative'
            elif 'Class C' in class_name and translated_generic_position in Gprot_positions:
                interaction = 'putative'
            elif 'Class C' in class_name and translated_generic_position_CtoB in Gprot_positions_B:
                interaction = 'putative'
            elif 'Class F' in class_name and translated_generic_position in Gprot_positions:
                interaction = 'putative'
            elif 'Class F' in class_name and translated_generic_position_FtoB in Gprot_positions_B:
                interaction = 'putative'
            else:
                interaction = None
        else:
            interaction = None

        return interaction

    dfm_drugs['GProteinInteraction'] = dfm_drugs.apply(gprotein_interface, axis=1)
    dfm_drugs['ArrestinInteraction'] = dfm_drugs.apply(arrestin_interface, axis=1)

    dfm_drugs_master = dfm_drugs.merge(mutations_pos, on=['SequenceNumber','EntryName','NMaa'], how='left').merge(diseasemap, on=['SequenceNumber','EntryName'], how='left').merge(ptms, on=['SequenceNumber','EntryName'], how='left').merge(lb_fam, on=['SequenceNumber','EntryName'], how='left').merge(lb_structure_human, on=['SequenceNumber','EntryName'], how='left').merge(lb_structure_ortho, on=['GPCRdb','EntryName'], how='left')

    ## CLASS A ONLY!
    ## make new calculation on additional functional sites by including:

    ## AJs paper: Diverse activation pathways in class A GPCRs converge near the G-protein-coupling region
    activation_pathway = ['3x46','6x37','7x53']

    ## PIF motif from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3644390/
    pif_motif = {'5x50':'P', '3x40':'I', '6x44':'F'}

    ## microswitches (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3343417/table/T1/)
    ## additional motifs/ conserved positions from the same table Asn1.50, Asp2.50, Pro2.59, Cys3.25, Asp3.32, Glu/Gln3.37, Leu3.40, Tyr3.51, Tyr3.60, Trp4.50, Phe6.44, Cys6.47, Asn7.45, Ser7.46
    micro_switches = {'3x49':'D,E', '7x43':'K,Y', '3x50':'R', '5x47':'F', '5x50':'P', '5x58':'Y', '6x30':'E', '6x34':'T', '6x48':'W', '6x50':'P', '7x49':'N', '7x50':'P', '7x53':'Y', '3x40':'I', '6x44':'F'} #includes the PIF motif!

    def micro_switch_func(df):
        if str(df.GPCRdb_short) in micro_switches:
            if str(df.WTaa) in micro_switches[df.GPCRdb_short]:
                return 'yes'
            else:
                return None
        else:
            return None

    ## SODIUM allosteric sites (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4106411/#SD1)
    sodium_pocket = {'1x50':'N','1x53':'V','2x46':'L','2x47':'A','2x49':'A','2x50':'D','3x39':'S','3x43':'L','6x44':'F','6x48':'W','7x45':'N','7x46':'S','7x49':'N','7x50':'P','7x53':'Y'}

    def sodium_pocket_func(df):
        if str(df.GPCRdb_short) in sodium_pocket:
            if str(df.WTaa) in sodium_pocket[df.GPCRdb_short]:
                return 'yes'
            else:
                return None
        else:
            return None

    dfm_drugs_master['ActivationPathway'] = dfm_drugs_master[dfm_drugs_master.Class=='Class A (Rhodopsin)'].GPCRdb_short.apply(lambda x: 'yes' if x in activation_pathway else None)
    dfm_drugs_master['MicroSwitch'] = dfm_drugs_master[['GPCRdb_short','WTaa']].apply(micro_switch_func, axis=1)
    dfm_drugs_master['SodiumPocket'] = dfm_drugs_master[['GPCRdb_short','WTaa']].apply(sodium_pocket_func, axis=1)

    func_put = dfm_drugs_master[(dfm_drugs_master.polyphen_score>=0.1) | (dfm_drugs_master.sift_score<=0.05)]
    func_known = dfm_drugs_master[(dfm_drugs_master.LB_fam=='interacting') | (dfm_drugs_master.PTMsite=='yes') | (dfm_drugs_master.GProteinInteraction=='putative') | (dfm_drugs_master.ArrestinInteraction=='putative') | (dfm_drugs_master['MicroSwitch']=='yes') | (dfm_drugs_master['SodiumPocket']=='yes') | (dfm_drugs_master['ActivationPathway']=='yes')]

    dfm_drugs_master['Functional'] = dfm_drugs_master.index.isin(func_put.index)
    dfm_drugs_master['func_known'] = dfm_drugs_master.index.isin(func_known.index)
    dfm_drugs_master['func_putative'] = dfm_drugs_master.index.isin(func_put.index)
    dfm_drugs_master['SegmentAgg'] = dfm_drugs_master['Segment'].apply(lambda x: 'TM' if 'TM' in x else 'EC-Loop' if 'ECL' in x else 'IC-Loop' if 'ICL' in x else 'H8' if 'H8' in x else 'terminus' if 'term' in x else 'None')
    [dfm_drugs_master.func_putative==True]
    dfm_drugs_master['VarantType'] = dfm_drugs_master['Allele Frequency'].apply(lambda x: 'CV' if x >= 0.001 else 'RV')
    dfm_drugs_master['codon_WT'] = dfm_drugs_master['k3_run'].str.split(':',expand=True)[0]

    dfm_drugs_master = dfm_drugs_master[dfm_drugs_master.columns[~dfm_drugs_master.columns.str.contains('Unnamed:')]]

    ## Remove Duplicates?
    dup = ['EntryName','SequenceNumber','NMaa','Protein Consequence','hgvsc_refine','pos_bin']
    dfm_drugs_master = dfm_drugs_master[~dfm_drugs_master.duplicated(dup, keep='first')]

    dfm_drugs_master.to_csv("../data/dfm_drugs_masterTable.csv")

    ## ==================== MERGING ====================
    ## Total number of known functional sites per receptor/ Functional Table!
    ## ptm + arrestin + gprotein + LB-family
    table = pd.read_csv("../data/all_GPCR_residues.csv")
    table.rename(columns={'amino_acid': 'WTaa'}, inplace=True)
    table = table.merge(ptms, on=['EntryName', 'SequenceNumber'], how='left').merge(lb_fam, on=['EntryName', 'SequenceNumber'], how='left')
    table['GPCRdb_short'] = table['GPCRdb'].str.split('.', expand=True)[0]+'x'+table['GPCRdb'].str.split('x', expand=True)[1]

    table['GProteinInteraction'] = table.apply(gprotein_interface, axis=1)
    table['ArrestinInteraction'] = table.apply(arrestin_interface, axis=1)

    table['ActivationPathway'] = table[table.Class=='Class A (Rhodopsin)'].GPCRdb_short.apply(lambda x: 'yes' if x in activation_pathway else None)
    table['MicroSwitch'] = table[['GPCRdb_short','WTaa']].apply(micro_switch_func, axis=1)
    table['SodiumPocket'] = table[['GPCRdb_short','WTaa']].apply(sodium_pocket_func, axis=1)

    table['Functional'] = table.index.isin(table[(table.LB_fam=='interacting') | (table.PTMsite=='yes') | (table.GProteinInteraction=='putative') | (table.ArrestinInteraction=='putative') | (table.ActivationPathway=='yes') | (table.MicroSwitch=='yes') | (table.SodiumPocket=='yes')].index)
    functionalTable = pd.DataFrame({'count_known_func' : table.groupby( ['EntryName','Functional'] ).size()}).reset_index()
    functionalTable = functionalTable[functionalTable.Functional==True]
    functionalTable.to_csv('../data/functionalTable.csv')

    table['PTMsite'] = table['PTMsite'].replace(np.nan, False)
    table['LB_fam'] = table['LB_fam'].replace('yes', True)
    table['GProteinInteraction'] = table['GProteinInteraction'].replace('yes', True)
    table['ArrestinInteraction'] = table['ArrestinInteraction'].replace('yes', True)
    table['ActivationPathway'] = table['ActivationPathway'].replace('yes', True)
    table['MicroSwitch'] = table['MicroSwitch'].replace('yes', True)
    table['SodiumPocket'] = table['SodiumPocket'].replace('yes', True)

    table = table[table.columns[~table.columns.str.contains('Unnamed:')]]
    drugs_raw = pd.read_csv("../data/drug_data.csv")
    drugTargets = drugs_raw[(drugs_raw.PMID.notnull()) & (drugs_raw['Status']=='approved')]['EntryName'].unique()
    table.to_csv("../data/all_GPCR_residues_functional_anno.csv")
    table[table.EntryName.isin(drugTargets)].to_csv('../data/drugtargets_functional_anno.csv')

    ## =================== Big merge END ================================
    dfm_drugs_master = pd.read_csv("../data/dfm_masterTable.csv")
    dfm_drugs_master = pd.read_csv("../data/dfm_masterTable_trials.csv")
    ## Venn maer: http://www.interactivenn.net/

    ## Sift score above threshold
    dfm_drugs_master[((dfm_drugs_master.polyphen_score>=0.1) | (dfm_drugs_master.sift_score<=0.05))].drop_duplicates().index.values

    ## Ligand interacting
    dfm_drugs_master[(dfm_drugs_master.LB_fam=='interacting')].index.values

    ## Changed property
    scored = dfm_drugs_master[((dfm_drugs_master.polyphen_score>=0.1) | (dfm_drugs_master.sift_score<=0.05))]
    dfm_drugs_master[dfm_drugs_master.MutationType=='changed'].drop_duplicates().index.values

    ## Either of any
    dfm_drugs_master[((dfm_drugs_master.GProteinInteraction=='putative') | (dfm_drugs_master.ArrestinInteraction=='putative')) | (dfm_drugs_master.score>=0.05) | (dfm_drugs_master.LB_fam=='interacting') | (dfm_drugs_master.maxFoldChange>=0.05)]

    ## Signalling interacting
    dfm_drugs_master[((dfm_drugs_master.GProteinInteraction=='putative') | (dfm_drugs_master.ArrestinInteraction=='putative'))].index.values

    ## MicroSwitch
    dfm_drugs_master[(dfm_drugs_master['MicroSwitch']=='yes') | (dfm_drugs_master['SodiumPocket']=='yes') | (dfm_drugs_master['ActivationPathway']=='yes')].index.values

    ## PTMsite
    dfm_drugs_master[(dfm_drugs_master.PTMsite=='yes')].index.values

    ## disease known
    dfm_drugs_master[(dfm_drugs_master.score>=0.005)].index.values

    ## InVitro Mutant
    dfm_drugs_master[(dfm_drugs_master.foldchangeMaxAbs>=5)].index.values

    ## clinical associations function
    dfm_drugs_master[(dfm_drugs_master.SequenceNumber==7) & (dfm_drugs_master.EntryName=='glpr1_human')]

    ## Has structure with a drug bound:
    drug_structures = pd.read_csv('../data/drug_structures.csv')
    scored[(scored.LB_fam=='interacting') & (scored.EntryName.isin(drug_structures.protein.unique()))]

    scored[(scored.maxFoldChange>=1.5) & (scored.EntryName.isin(drug_structures.protein.unique()))]
    dfm_drugs_master[(dfm_drugs_master.MutationType=='changed') & (scored.EntryName.isin(drug_structures.protein.unique())) & ((dfm_drugs_master.score>=0.1) | (dfm_drugs_master.maxFoldChange>=1.5))]

    ## How much on known funcitonal sites Introduction
    round(dfm_drugs_master[(dfm_drugs_master.LB_fam=='interacting') | (dfm_drugs_master.PTMsite=='yes') | (dfm_drugs_master.GProteinInteraction=='putative') | (dfm_drugs_master.ArrestinInteraction=='putative') | (dfm_drugs_master['MicroSwitch']=='yes') | (dfm_drugs_master['SodiumPocket']=='yes')].shape[0]/float(dfm_drugs_master.shape[0])*100,2)

    ## How much on putative funcitonal sites Introduction
    # AND
    round(dfm_drugs_master[(dfm_drugs_master.LB_fam=='interacting') | (dfm_drugs_master.PTMsite=='yes') | (dfm_drugs_master.GProteinInteraction=='putative') | (dfm_drugs_master.ArrestinInteraction=='putative') | (dfm_drugs_master['MicroSwitch']=='yes') | (dfm_drugs_master['SodiumPocket']=='yes') | ((dfm_drugs_master.polyphen_score>=0.1) | (dfm_drugs_master.sift_score<=0.05))].shape[0]/float(dfm_drugs_master.shape[0])*100,2)

    putative_functional.to_csv("../data/masterTable_putative_functional.csv")

    dfm_drugs_master['Region'] = dfm_drugs_master['Segment'].apply(lambda x: 'EC-Loop' if 'ECL' in str(x) else 'TM' if 'TM' in str(x) else 'IC-Loop' if 'ICL' in str(x) else 'N-term' if 'N-term' in str(x) else 'H8' if 'H8' in str(x) else 'C-term' if 'C-term' in str(x) else 'Other')
    functional_count = pd.DataFrame({'count' : dfm_drugs_master.groupby( ['EntryName','Region','func_putative'] ).size()}).reset_index()
    functional_count['func_putative'] = np.where(functional_count['func_putative'], 'putative functional', 'unknown')
    ## Export
    functional_count.to_csv("../data/functional_count.csv")
    dfm_drugs_master.to_csv("../data/dfm_drugs_masterTable.csv")

    ## Average Rare and Common variants
    dfm_drugs_master[dfm_drugs_master['VarantType']=='CV'].groupby(["EntryName"]).size().mean()

    ## MV and invitro mutations
    mutations.head()
    mutations.rename(columns={'Mutantaa': 'NMaa'}, inplace=True)
    overlap = dfm_drugs_master.merge(mutations, on=['EntryName','SequenceNumber','WTaa', 'NMaa'])
    overlap[['EntryName','SequenceNumber','WTaa', 'NMaa']].drop_duplicates().shape

    ## polymorphism-induced alteration in drug response in at least
    ## XX % of targeted GPCRs
    dfm_drugs_master[((dfm_drugs_master.LB_fam=='interacting') | (dfm_drugs_master.PTMsite=='yes') | (dfm_drugs_master.GProteinInteraction=='putative') | (dfm_drugs_master.ArrestinInteraction=='putative')) & ((dfm_drugs_master.polyphen_score>=0.1) & (dfm_drugs_master.sift_score<=0.05))]['ac'].mean()

    ## For instance, X of the YYY positions that constitute the binding pocket for (Drug) in the X receptor shows polymorphism
    dfm_drugs_master[(dfm_drugs_master.GProteinInteraction=='putative')].GPCRdb.value_counts()

    # Gprotein For instance, XX of the 108 GPCR drug targets harbour a mutation in position 3.X
    dfm_drugs_master[(dfm_drugs_master.GProteinInteraction=='putative')].GPCRdb.value_counts()
    # Arrestin
    dfm_drugs_master[(dfm_drugs_master.ArrestinInteraction=='putative')].GPCRdb.value_counts()

    ## Functional sites:
    dfm_drugs_master[(dfm_drugs_master['ArrestinInteraction'] == 'putative') | (dfm_drugs_master['GProteinInteraction'] == 'putative') | (dfm_drugs_master.LB_fam=='interacting') | (dfm_drugs_master.PTMsite=='yes') | (dfm_drugs_master['MicroSwitch']=='yes') | (dfm_drugs_master['SodiumPocket']=='yes')]

    ## unique receptor positions
    dfm_drugs_master[['EntryName','SequenceNumber']].drop_duplicates().shape

def drugs_pairing_analysis():
    ### TARGET DIFFERENTIATION

    drugs_raw = pd.read_csv('../data/drug_data.csv')
    ## approved primary
    drugs_raw[(drugs_raw['Status']=='approved') & (drugs_raw.TargetCategory=='primary')]['EntryName'].nunique()
    ## in trials
    len(trialTargets)
    ## in trials and primary
    len(np.array(list(set(drugs_raw[(drugs_raw.TargetCategory=='primary')]['EntryName'].unique()) - set(drugs_raw[drugs_raw['Status']=='approved']['EntryName'].unique()))))
    ## in trials and primary and ongoing trials
    len(np.array(list(set(drugs_raw[(drugs_raw.TargetCategory=='primary') & (drugs_raw.ClinicalStatus.isin(['Completed', 'recruiting', 'Ongoing', 'Not open yet', 'suspended']))]['EntryName'].unique()) - set(drugs_raw[drugs_raw['Status']=='approved']['EntryName'].unique()))))
    # )

def multi_MV_positions():
    """
    multiple single protein position mutations possible and also supported
    in the data.
    Which positions have multiple nucleotide changes leading to the same or
    different amino acid mutation?
    """
    proteins = snp['EntryName'].unique()
    for protein in proteins:
        SMs = list(snp[snp['EntryName']==protein]['site'])
        doubles = list(set([x for x in SMs if SMs.count(x) == 2]))
        triples = list(set([x for x in SMs if SMs.count(x) == 3]))
        quadruples = list(set([x for x in SMs if SMs.count(x) == 4]))
        quintet = list(set([x for x in SMs if SMs.count(x) == 5]))
        # no six plus mutation sites!
        if doubles:
            print protein, doubles

def allele_frequencies():
    ## Overview of distributions of allele_frequencies in dataset

    # Common Variant (CV) >= 0.001
    # Rare Variant (RV) <= 0.001
    snp = pd.read_csv("../data/dfm_masterTable.csv", low_memory=False)

    ## Update MutationType
    AA_GROUPS = OrderedDict([
        ('hp',  ('A', 'C', 'F', 'I', 'L', 'M', 'V', 'W', 'Y')),
        ('ar',  ('F', 'H', 'W', 'Y')),
        ('pu', ('S','T','N','Q')),
        ('pro', ('P','G')),
        ('neg', ('D', 'E')),
        ('pos', ('H', 'K', 'R')),
    ])
    # ('sma'), ('G','A','V','C')),

    three_letter = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
      'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
      'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
      'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
      'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

    for index, rows in tqdm(enumerate(snp.iterrows())):
        WT = snp[index:index+1]['WTaa'].values[0]
        NM = snp[index:index+1]['NMaa'].values[0]
        if len(set([i for i in AA_GROUPS if WT in AA_GROUPS[i]]).intersection([i for i in AA_GROUPS if NM in AA_GROUPS[i]]))>0:
            Mtype = 'similar'
        else:
            Mtype = 'changed'

        snp.loc[index,'MutationType'] = Mtype

    ## Aggregate allele variants for the same amino acid position:
    snp['PositionAF'] = snp.groupby(['EntryName','SequenceNumber'])['Allele Frequency'].cumsum()

    ## calculate Variant Type for the aggregated protein position allele frequency
    snp['VarantType'] = snp['PositionAF'].apply(lambda x: 'CV' if x >= 0.001 else 'RV')
    f = lambda x: 1 if x>0 else 0 if x ==0 else -1
    snp['SegmentType'] = snp['Segment'].apply(lambda x: 'Loop' if 'CL' in str(x) else 'TM' if 'TM' in str(x) else 'Other')

    Singletons = snp[snp['Allele Count']==1]['MutationType'].value_counts().to_dict()
    AF_001 = snp[(snp['Allele Frequency']<=0.0001) & (snp['Allele Count']>1)]['MutationType'].value_counts().to_dict()
    AF_01 = snp[(snp['Allele Frequency']>0.001) & (snp['Allele Frequency']<=0.01)]['MutationType'].value_counts().to_dict()
    AF_1 = snp[(snp['Allele Frequency']>0.01) & (snp['Allele Frequency']<=0.1)]['MutationType'].value_counts().to_dict()
    AF_10 = snp[(snp['Allele Frequency']>0.1)]['MutationType'].value_counts().to_dict()

    stacked_plot = pd.DataFrame(snp['MutationType'].value_counts().to_dict(), index=['All'])
    stacked_plot.loc['Singletons'] = pd.Series(Singletons)
    stacked_plot.loc['<0.01 %'] = pd.Series(AF_001)
    stacked_plot.loc['0.1 %'] = pd.Series(AF_01)
    stacked_plot.loc['1 %'] = pd.Series(AF_1)
    stacked_plot.loc['>10 %'] = pd.Series(AF_10)

    sns.set(style="white", color_codes=True, font_scale=1.4)
    f,(ax,ax2) = plt.subplots(1,2,sharey=True, facecolor='w')
    reds = np.array([[0.914,0.11,0.141,1.0],[0.612,0.141,0.1378,1.0]])
    stacked_plot.plot(kind='barh', ax=ax, stacked=True, legend=False, color=reds)
    stacked_plot.plot(kind='barh', ax=ax2, stacked=True, color=reds).legend(bbox_to_anchor=(0.95, 0.2))

    ## hide the spines between ax and ax2
    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax.yaxis.tick_left()
    ax.tick_params(labelright='off')
    ax2.yaxis.tick_right()
    d = .015 # how big to make the diagonal lines in axes coordinates
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((1-d,1+d), (-d,+d), **kwargs)
    ax.plot((1-d,1+d),(1-d,1+d), **kwargs)

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d,+d), (1-d,1+d), **kwargs)
    ax2.plot((-d,+d), (-d,+d), **kwargs)

    ax.set_xlim(0,300)
    ax2.set_xlim(1000,15000)
    # ax2.yaxis.set_major_formatter(plt.NullFormatter())

    ax.set_xticks([0, 100, 200])
    ax2.set_xticks([5e3, 10e3, 15e3])

    plt.gca().invert_yaxis()
    plt.xlabel('Number missense variations')
    plt.ylabel('Allele frequencies')
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=0)
    ax.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    plt.savefig("../figures/overview_frequencies.svg")
    plt.show()

    ## Across all proteins for generic positions
    snp['GenericPositionAF'] = snp.groupby(['GPCRdb'])['Allele Frequency'].cumsum()
    snp.sort_values(by='Allele Frequency')[['EntryName','Allele Frequency', 'VarantType','Protein Consequence']]

def receptor_structural_regions():
    ## absolute and relative abundance of MVs in secondary structure elements

    def get_sublist(EntryName, Region):
        value = 0
        try:
            value = segment_lengths[segment_lengths['EntryName']==EntryName][Region].values[0]
        except Exception as msg:
            print msg
        return value

    snp = pd.read_csv("../data/ExAC_GPCRs_raw.csv")
    snp = pd.read_csv("../data/dfm_masterTable.csv")
    snp = snp[snp['Class']!='G-Protein']
    # Aggregate lengths for selected structural regions
    segment_lengths = pd.read_csv("../data/GPCR_segment_lengths.csv")
    segment_lengths['Loop'] = segment_lengths[segment_lengths.columns[segment_lengths.columns.str.contains('CL')]].sum(axis=1)
    segment_lengths['TM'] = segment_lengths[segment_lengths.columns[segment_lengths.columns.str.contains('TM')]].sum(axis=1)
    segment_lengths['EC-Loop'] = segment_lengths[segment_lengths.columns[segment_lengths.columns.str.contains('ECL')]].sum(axis=1)
    segment_lengths['IC-Loop'] = segment_lengths[segment_lengths.columns[segment_lengths.columns.str.contains('ICL')]].sum(axis=1)
    segment_lengths['H8'] = segment_lengths[segment_lengths.columns[segment_lengths.columns.str.contains('H8')]].sum(axis=1)

    # Count number of SNPs in assigned regions
    snp['Region'] = snp['Segment'].apply(lambda x: 'EC-Loop' if 'ECL' in str(x) else 'TM' if 'TM' in str(x) else 'IC-Loop' if 'ICL' in str(x) else 'N-term' if 'N-term' in str(x) else 'H8' if 'H8' in str(x) else 'C-term' if 'C-term' in str(x) else 'Other')

    snp_count = pd.DataFrame({'snp_count_CV' : snp[snp['Allele Frequency']>=0.001].groupby( ['EntryName','Region','MutationType'] )['SequenceNumber'].nunique(),
    'snp_count_RV' : snp[(snp['Allele Frequency']<0.001)].groupby( ['EntryName','Region','MutationType'] )['SequenceNumber'].nunique()}).reset_index() # RV
    snp_count = snp_count[snp_count['Region']!='Other']
    snp_count['RegionLength'] = map(get_sublist,snp_count['EntryName'],snp_count['Region'])
    snp_count['RelativeSNP_RV'] = snp_count.apply(lambda x: round(float(x['snp_count_RV'])/x['RegionLength'],2), axis=1)
    snp_count['RelativeSNP_CV'] = snp_count.apply(lambda x: round(float(x['snp_count_CV'])/x['RegionLength'],2), axis=1)

    snp_count.Region = snp_count.Region.astype("category")
    labels = ['C-term', 'EC-Loop', 'TM', 'IC-Loop', 'H8', 'N-term']
    snp_count.Region.cat.set_categories(labels, inplace=True)

    f, (ax,ax2) = plt.subplots(figsize=(26, 12), ncols=2)
    sns.violinplot(x="Region", y="RelativeSNP_RV", hue="MutationType", palette="gray", split=True, data=snp_count, inner='box', cut=0, bw=.3, ax=ax)
    sns.violinplot(x="Region", y="RelativeSNP_CV", hue="MutationType", palette="gray", split=True, data=snp_count, inner='box', cut=0, bw=.3, ax=ax2)
    # sns.boxplot(x="Region", y="RelativeSNP_RV", data=snp_count, ax=ax, palette="Blues_r")
    # sns.boxplot(x="Region", y="RelativeSNP_CV", data=snp_count, ax=ax2, palette="Reds_r")
    # ax.set_title('rare variants')
    # ax2.set_title('common variants')

    ax.set_ylim(ax.get_ylim()[::-1])
    ax2.set_ylim(ax2.get_ylim()[::-1])
    ax.xaxis.tick_top()
    ax2.xaxis.tick_top()
    ax.annotate('rare variants', xy=(0.5, 0.5), xytext=(0.5, 0.5))
    ax2.annotate('common variants', xy=(0.5, 0.5), xytext=(0.5, 0.5))

    # ax.set_ylim([0,0.4])
    # ax2.set_ylim([0,0.4])
    plt.show()

    ## ======

    snp_count = pd.DataFrame({'snp_count' : snp[snp['Allele Frequency']<=0.001].groupby( ['EntryName','Region','Ligandtype'] )['SequenceNumber'].nunique(),'ac_sum' : snp[snp['Allele Frequency']<=0.001].groupby( ['EntryName','Region','Ligandtype'] )['ac'].sum()}).reset_index()
    snp_count['RegionLength'] = map(get_sublist,snp_count['EntryName'],snp_count['Region'])
    snp_count.head()
    snp_count['snp_density'] = snp_count.apply(lambda x: round(float(x['snp_count'])/x['RegionLength'],2), axis=1)
    snp_count['ac_density'] = snp_count.apply(lambda x: round(float(x['ac_sum'])/x['snp_count'],2), axis=1)
    df = snp_count[['Region','Ligandtype','snp_density']]


    ## DATA PREP
    snp_count.ix[(snp_count.Ligandtype=='Protein'),'Ligandtype']='Peptide'
    snp_count[['EntryName','Ligandtype']].drop_duplicates().Ligandtype.value_counts()
    sns.set(style="ticks", color_codes=True, font_scale=1.4)
    ## Functional sites absolute counts


    functional_count = pd.DataFrame({'count' : dfm_drugs_master.groupby( ['EntryName','Region','func_putative'] ).size()}).reset_index()
    functional_count['RegionLength'] = map(get_sublist,functional_count['EntryName'],functional_count['Region'])
    functional_count['snp_density'] = functional_count.apply(lambda x: round(float(x['count'])/x['RegionLength'],2), axis=1)
    functional_count.ix[(functional_count.func_putative==True),'func_putative']='Polymorphism in putative functional site'
    functional_count.ix[(functional_count.func_putative==False),'func_putative']='Polymorphism in uncharacterised sites'
    functional_count.to_csv("../data/functional_count.csv")

    ## VIZ
    f, (ax,ax2) = plt.subplots(figsize=(13, 12), nrows=2, sharex=True)
    ax2 = sns.boxplot(x="Region", y='snp_density',hue="func_putative", data=functional_count, palette=['#3aa74d','#dddddd'], ax=ax2, linewidth=0.7, notch=False)
    ax2.xaxis.set_tick_params(labelsize=19)
    ax2.set_xlabel('')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.xaxis.tick_top()
    ax2.set_ylabel('fraction of missense variants')
    l = ax2.legend()
    l.set_title('')
    ## Rebuttal Ligandtype:
    # sns.boxplot(x="Region", y='snp_density',hue="Ligandtype", data=snp_count[(snp_count.Ligandtype.isin(['Peptide', 'Protein', 'Lipid', 'Aminergic',  'Nucleotide'])) & (snp_count.Region.isin(['N-term','C-term','EC-Loop', 'H8', 'IC-Loop', 'TM']))], palette="BrBG", ax=ax, linewidth=0.7, notch=False)
    # # ['Peptide', 'Protein', 'Lipid', 'Aminergic', 'Amino acid', 'Alicarboxylic acid', 'Nucleotide', 'Melatonin']
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # ax.spines['bottom'].set_visible(False)
    # ax.set_xlabel('')
    # ax.set_ylabel('missense variation density in segment')
    # ax.set_ylim([0,0.65])
    # l = ax.legend()
    # l.set_title('')
    # bottom
    sns.set(style="ticks", color_codes=True, font_scale=1.4)
    ax = sns.boxplot(x="Region", y='count',hue="func_putative", data=functional_count, palette=['#3aa74d','#dddddd'], ax=ax, linewidth=0.7, notch=False)
    # ['Peptide', 'Protein', 'Lipid', 'Aminergic', 'Amino acid', 'Alicarboxylic acid', 'Nucleotide', 'Melatonin']
    ax.xaxis.set_tick_params(labelsize=19)
    ax.set_xlabel('')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.xaxis.tick_top()
    ax.set_ylabel('no. of missense variants')

    l = ax.legend()
    l.set_title('')
    # ax2.xaxis.set_ticks_position('both')

    hspace = 0.25   # the amount of height reserved for white space between subplots
    plt.subplots_adjust(hspace=hspace)
    plt.savefig("../figures/mv_distributions.svg", bbox_inches='tight')
    plt.show()

def Gprotein_interface_matrix():
    ## interface matrix to prioritize interface positions

    def Gprotein_contacts(x):

        try:
            Gprotein_pos = interface_matrix[interface_matrix.GPCRdb==x].notnull().any(axis=0)[interface_matrix[interface_matrix.GPCRdb==x].notnull().any(axis=0)].keys().drop('GPCRdb')
            structures = interface_matrix[interface_matrix.GPCRdb==x][Gprotein_pos].values[0]
            return len(Gprotein_pos), ','.join(list(set(structures))), ','.join(list(set(Gprotein_pos)))
        except:
            return 0, 'Gt?', 'Gt?'

    interface_matrix = pd.read_excel("../data/GPCR-GProtein_interface_matrix.xlsx")
    dfm_drugs_master = pd.read_csv("../data/dfm_drugs_masterTable.csv")
    dfm_drugs_master['codon_WT'] = dfm_drugs_master['k3_run'].str.split(':',expand=True)[0]
    export=['number_of_contacts','structures','positions_on_the_G_protein','codon_WT','Transcript Consequence']

    SB = dfm_drugs_master[dfm_drugs_master.EntryName.isin(['cckar_human']) & ((dfm_drugs_master.ArrestinInteraction=='putative') | (dfm_drugs_master.GProteinInteraction=='putative'))]

    SB['number_of_contacts'], SB['structures'], SB['positions_on_the_G_protein'] = zip(*SB["GPCRdb"].map(Gprotein_contacts))
    SB.sort_values(by=['Allele Frequency'], ascending=False)[export]

def drugs_MVs():
    ## For each drug get the number of SNPS in functional sites
    ## Drugs vs. MVs for input in R vis

    drugs_raw = pd.read_csv("../data/drug_data.csv")
    dfm_drugs_master = pd.read_csv("../data/dfm_masterTable.csv")

    drugTargets = drugs_raw[(drugs_raw.PMID.notnull()) & (drugs_raw['Status']=='approved')]['EntryName'].unique()
    drugs = pd.DataFrame({'Approved.Drugs' : drugs_raw[(drugs_raw.PMID.notnull()) & (drugs_raw['Status']=='approved')].groupby( ['EntryName'] )['Drug Name'].nunique()}).reset_index()

    relative_table = pd.read_csv("../data/relative_table.csv")
    relative_table = relative_table[relative_table.EntryName.isin(drugTargets)]

    exac_GPCR = pd.read_csv("../data/exac_GPCRs.csv")

    lof = pd.read_csv("../data/SNP_LoF_ExAC_ENSG.csv")
    lof['ac_het'] = lof['Allele Count']-lof['Number of Homozygotes']
    lof_target = pd.DataFrame({'LoF_positions' : lof.groupby( ['EntryName'] ).size(),
            'ac_hom' : lof.groupby( ['EntryName'] )['Number of Homozygotes'].sum(),
            'ac' : lof.groupby( ['EntryName'] )['Allele Count'].sum(),
            'min_ind' : lof.groupby( ['EntryName'] )['ac_het'].max(),
            'max_ac' : lof.groupby( ['EntryName'] )['Allele Count'].max()}).reset_index()
    lof_target['max_ind'] = lof_target['ac']-lof_target['ac_hom']
    lof_target['LoF_count_min_freq'] = lof_target['min_ind']/60760*100
    lof_target['LoF_count_max_freq'] = lof_target['max_ind']/60760*100
    lof_target = lof_target[lof_target.EntryName.isin(drugTargets)]
    # sns.lmplot(x="n_lof", y="ac", data=lof_target, fit_reg=False)

    merge = pd.merge(left=drugs,right=relative_table, how='right', left_on='EntryName', right_on='EntryName').merge(lof_target, on='EntryName')
    merge = pd.merge(left=merge,right=exac_GPCR, how='left', on = ['EntryName'])
    # merge[pd.isnull(merge).any(axis=1)]
    # merge = merge.dropna()
    merge['ProteinID'] = merge['ProteinID'].str.upper()

    functionalTable = pd.read_csv('../data/functionalTable.csv')
    pfdf = pd.DataFrame({'n_mis_functional' : dfm_drugs_master[dfm_drugs_master.func_known==True].groupby( ['EntryName'] )['SequenceNumber'].nunique()}).reset_index() # functional

    Func = pfdf.merge(functionalTable, on = 'EntryName')
    Func['Relative_functional'] = Func.apply(lambda x: round(float(x['n_mis_functional'])/x['count_known_func'],2), axis=1)

    merge = merge.merge(Func, on = 'EntryName')
    # merge['RelativeSNP_func'] = merge.apply(lambda x: round(float(x['n_mis_functional'])/x['seqLen'],2), axis=1)
    merge.sort_values(by='Relative_functional')
    merge = merge.dropna()
    merge['ProteinID'] = merge.EntryName.str.split('_', expand=True)[0].str.upper()
    merge.to_csv('../data/ratio_funcational_sites.csv')

    # Excel export
    from pandas import ExcelWriter
    writer = ExcelWriter('~/Downloads/Supp7_receptors_known.xlsx')
    merge[['EntryName','Relative_functional']].to_excel(writer, 'known', index=False)
    writer.save()

    # missing:
    dfm_drugs_master[dfm_drugs_master.EntryName=='ada2a_human']

def drugs_CNVs():
    ## Drugs vs. CNVs for input in R vis

    drugs_raw = pd.read_csv("../data/drug_data.csv")

    drugTargets = drugs_raw[(drugs_raw.PMID.notnull()) & (drugs_raw['Status']=='approved')]['EntryName'].unique()
    drugs = pd.DataFrame({'Approved.Drugs' : drugs_raw[(drugs_raw.PMID.notnull()) & (drugs_raw['Status']=='approved')].groupby( ['EntryName'] )['Drug Name'].nunique()}).reset_index()

    relative_table = pd.read_csv("../data/relative_table.csv")
    cnf = pd.read_csv('../data/CNV/cnv_GPCRs.csv')

    merge = pd.merge(left=drugs,right=cnf, how='left', left_on='EntryName', right_on='EntryName')
    merge['CNVs'] = merge['del'] + merge['dup']
    merge = merge.dropna()

    merge.to_csv('../data/Drugs_vs_CNVs.csv')

def estimated_economic_burden():
    ## Calculations of estimated economic burden_p

    # (NHS sales cost per month) x (% people affected by a drug from variants in putative functional sites from 1k genotype data) x (correction factor for side effects as not all costs can be associated with ADRs)
    # in numbers it could be something like: 10mio GBP x 0.32 x 0.1 = 0.32 mio GBP per month (for one drug)
    # correction factor for sales/person + other effects including ADR and dose/efficacy effects.

    nhs_raw_target = pd.read_csv("../data/NHS_raw.csv")
    nhs_raw = nhs_raw_target[['actual_cost', 'date', 'drugClass', 'drugCode', 'drugName', 'drugNameQuery', 'items', 'quantity', 'section']].drop_duplicates()
    nhs_raw['year'] = pd.DatetimeIndex(nhs_raw['date']).year
    nhs_raw['section'] = nhs_raw['section'].str.split(': ', expand=True)[1]

    # nhs_raw.replace(to_replace="Drugs Used In Park'ism/Related Disorders", value="Drugs Used In Parkinson/Related Disorders", inplace=True)
    # nhs_raw.replace(to_replace="Drugs Used In Psychoses & Rel.Disorders", value="Drugs Used In Psychoses & Relelated Disorders", inplace=True)
    # nhs_raw.replace(to_replace="Sex Hormones & Antag In Malig Disease", value="Hormones & Antagonists In Malignant Disease", inplace=True)
    # nhs_raw.replace(to_replace="Hypothalamic&Pituitary Hormones&Antioest", value="Hypothalamic & Pituitary Hormones", inplace=True)
    # nhs_raw['section'] = nhs_raw['section'].str.lower()
    # nhs_raw.replace(to_replace="antispasmod.&other drgs alt.gut motility", value="drugs altering gut motility", inplace=True)
    # nhs_raw.replace(to_replace='antihist, hyposensit & allergic emergen', value="antihistamines and antiallergics", inplace=True)

    nhs = pd.DataFrame({
        'yearly_items' : nhs_raw.groupby( ['drugNameQuery', 'section', 'year'] )['items'].sum(),
        'yearly_cost' : nhs_raw.groupby( ['drugNameQuery', 'section', 'year'] )['actual_cost'].sum()}).reset_index()
    drugsdf = pd.read_csv('../data/drugsdf_putative.csv')
    drugsdf_homo = pd.read_csv('../data/drugsdf_putative_homo.csv')
    drugsdf_known = pd.read_csv('../data/drugsdf_known.csv')
    drugsdf_known_homo = pd.read_csv('../data/drugsdf_known_homo.csv')
    # drugsdf = drugsdf_known

    burden_p = drugsdf.merge(nhs, left_on = ['drug'], right_on = ['drugNameQuery'])
    burden_p = burden_p[(burden_p.year<2017) & (burden_p.year>2012)] # only 4 complete years!
    burden_p = burden_p[['drug', u'individuals', u'percentage', u'yearly_cost', u'yearly_items', 'section']]
    burden_p['estimated_burden_p'] = burden_p['yearly_cost'] * burden_p['percentage']/100 * 1

    burden_p_homo = drugsdf_homo.merge(nhs, left_on = ['drug'], right_on = ['drugNameQuery'])
    burden_p_homo = burden_p_homo[(burden_p_homo.year<2017) & (burden_p_homo.year>2012)]
    burden_p_homo = burden_p_homo[['drug', u'individuals', u'percentage', u'yearly_cost', u'yearly_items', 'section']]
    burden_p_homo['estimated_burden_p_homo'] = burden_p_homo['yearly_cost'] * burden_p_homo['percentage']/100 * 1

    burden_k = drugsdf_known.merge(nhs, left_on = ['drug'], right_on = ['drugNameQuery'])
    burden_k = burden_k[(burden_k.year<2017) & (burden_k.year>2012)]
    burden_k = burden_k[['drug', u'individuals', u'percentage', u'yearly_cost', u'yearly_items', 'section']]
    burden_k['estimated_burden_k'] = burden_k['yearly_cost'] * burden_k['percentage']/100 * 1

    burden_k_homo = drugsdf_known_homo.merge(nhs, left_on = ['drug'], right_on = ['drugNameQuery'])
    burden_k_homo = burden_k_homo[(burden_k_homo.year<2017) & (burden_k_homo.year>2012)]
    burden_k_homo = burden_k_homo[['drug', u'individuals', u'percentage', u'yearly_cost', u'yearly_items', 'section']]
    burden_k_homo['estimated_burden_k_homo'] = burden_k_homo['yearly_cost'] * burden_k_homo['percentage']/100 * 1

    burden_agg = pd.DataFrame({
    'av_year_p' : burden_p.groupby( ['drug','section','percentage'] )['yearly_cost'].mean(),
    'av_year_p_homo' : burden_p_homo.groupby( ['drug','section','percentage'] )['yearly_cost'].mean(),
    'av_year_k' : burden_k.groupby( ['drug','section','percentage'] )['yearly_cost'].mean(),
    'av_year_k_homo' : burden_k_homo.groupby( ['drug','section','percentage'] )['yearly_cost'].mean()}).reset_index()

    # burden_agg.sort_values(by='estimated_burden_1.0')

    # for i in [0.1,0.5,1.0]:
    #     burden_agg['estimated_burden_'+str(i)] = burden_agg['av_year'] * burden_agg['percentage'] * i

    for i in ['p','p_homo','k','k_homo']:
        burden_agg['estimated_burden_'+str(i)] = burden_agg['av_year_'+i] * burden_agg['percentage']/100 * 1

    #
    # burden_section = pd.DataFrame({'0.1' : burden_agg.groupby( ['section'] )['estimated_burden_0.1'].sum(),
    # '0.5' : burden_agg.groupby( ['section'] )['estimated_burden_0.5'].sum(),
    # '1.0' : burden_agg.groupby( ['section'] )['estimated_burden_1.0'].sum()}).reset_index()
    # burden_section.sort_values(by='0.1').tail()

    burden_section = pd.DataFrame({'p' : burden_agg.groupby( ['section'] )['estimated_burden_p'].sum(),
    'p_homo' : burden_agg.groupby( ['section'] )['estimated_burden_p_homo'].sum(),
    'k' : burden_agg.groupby( ['section'] )['estimated_burden_k'].sum(),
    'k_homo' : burden_agg.groupby( ['section'] )['estimated_burden_k_homo'].sum()}).reset_index()

    # pd.melt(combined_regions, id_vars=['EntryName','Region'], value_vars=['1-Sift','Polyphen'])
    # top10 = list(burden_section.sort_values(by='0.1').tail(n=10)['section'])
    # burden_section['section'] = burden_section['section'].apply(lambda x: x if x in top10 else 'other')

    # Rreshape = pd.melt(burden_section, id_vars=['section'], value_vars=['0.1','0.5','1.0'])
    common_sections = burden_section[burden_section.k > 10000000].dropna().section.unique()
    Rreshape = pd.melt(burden_section, id_vars=['section'], value_vars=['p','p_homo','k','k_homo'])

    Rreshape.ix[(Rreshape.value<5000000), 'section']='other'
    Rreshape = pd.DataFrame({'value' : Rreshape.groupby( ['section', 'variable'] )['value'].sum()}).reset_index()
    print(Rreshape.section.nunique())

    final = []
    var_translate = {'k_homo':'known-homozygous','k':'known-all variants','p_homo':'putative-homozygous','p':'putative-all variants'}
    for section in Rreshape.section.unique():
        data = []
        for var in ['k_homo','k','p_homo','p']:
            try:
                y = Rreshape[(Rreshape.section==section) & (Rreshape.variable==var)]['value'].values[0]
            except:
                y = 0
            # x = Rreshape[index:index+1]['variable'].values[0]
            data.append({'y': int(y), 'x': var_translate[var]})
        key = section
        final.append({'values': data, 'key': key})
    print final

    # Rreshape = Rreshape[(Rreshape.value>0) & (Rreshape.section.isin(common_sections))].reset_index() #
    Rreshape.to_csv("../data/economic_burden.csv")

    print Rreshape.groupby('variable').sum()/1000000000

    from pandas import ExcelWriter
    burden_agg = pd.DataFrame({
    'putative_all_variants' : burden_p.groupby( ['drug','section'] )['estimated_burden_p'].mean(),
    'putative_homozygous' : burden_p_homo.groupby( ['drug','section'] )['estimated_burden_p_homo'].mean(),
    'known_all_variants' : burden_k.groupby( ['drug','section'] )['estimated_burden_k'].mean(),
    'known_homozygous' : burden_k_homo.groupby( ['drug','section'] )['estimated_burden_k_homo'].mean(),
    'average_yearly_cost' : burden_k_homo.groupby( ['drug','section'] )['yearly_cost'].mean(),
    'average_yearly_items' : burden_k_homo.groupby( ['drug','section'] )['yearly_items'].mean()}).reset_index()
    burden_agg.sort_values(by='drug')
    writer = ExcelWriter('~/Downloads/Supp11_estimated_economic_burden.xlsx')
    burden_agg.to_excel(writer, index=False)
    writer.save()

    ### oprm
    oprm_drugs = nhs_raw_target[nhs_raw_target.drugTarget=='oprm_human'].drugNameQuery.unique()
    burden_agg[burden_agg.drug.isin(oprm_drugs)].av_year_p_homo.sum()
    burden_agg[burden_agg.drug.isin(oprm_drugs)].estimated_burden_p_homo.sum()


    print burden_agg.estimated_burden_k_homo.sum()/1000000
    print burden_agg.estimated_burden_k.sum()/1000000
    print burden_agg.estimated_burden_p_homo.sum()/1000000
    print burden_agg.estimated_burden_p.sum()/1000000

    ### ED7 plot
    burden_agg.index = burden_agg['drug']
    burden_agg[['estimated_burden_1.0']].sort_values(by='estimated_burden_1.0', ascending=True).tail(n=15).plot.barh(color='#db7f2d')
    plt.show()
    plt.savefig("../figures/top_economicburden.svg")

def export():
    ## export to Excel tables
    from pandas import ExcelWriter

    ## PERMUATION TEST (SREENI)
    dfm_drugs_master[['Ligandtype','EntryName','GPCRdb','SequenceNumber','Segment','SegmentAgg','func_known','func_putative','ac']].to_csv('../permutation_input.csv')

    ## MASTER All
    writer = ExcelWriter('~/Downloads/Supp4_Master_all_variants.xlsx')
    dfm_drugs_master[["EntryName","Family","Ligandtype","Class","Uniprot","gene","GPCRdb","SequenceNumber","Segment","GPCRdbWT","NMaa","MutationType","Allele Count","Allele Frequency","Allele Number","Number of Homozygotes","sift_word","sift_score","polyphen_word","polyphen_score","GProteinInteraction","ArrestinInteraction","ActivationPathway","MicroSwitch","SodiumPocket","foldchangeMaxAbs","diseaseId","score","diseaseName","Type","PTMsite","LB_structure","LB_fam"]].to_excel(writer, 'All', index=False)
    writer.save()

    ### Functional sites
    func_sites = ExcelWriter('~/Downloads/Supp5_functional_sites.xlsx')
    ligand_binding = dfm_drugs_master[(dfm_drugs_master.LB_fam=='interacting')]
    ligand_binding[["EntryName","Family","Ligandtype","Class","Uniprot","gene","GPCRdb","SequenceNumber","Segment","GPCRdbWT","NMaa","MutationType","Allele Count","Allele Frequency","Allele Number","Number of Homozygotes","sift_word","sift_score","polyphen_word","polyphen_score","foldchangeMaxAbs","diseaseId","score","diseaseName","Type","PTMsite","LB_structure","LB_fam","GProteinInteraction","ArrestinInteraction","ActivationPathway","MicroSwitch","SodiumPocket"]].to_excel(func_sites, 'Ligand Binding', index=False)

    signalling_binding = dfm_drugs_master[(dfm_drugs_master.GProteinInteraction=='putative') | (dfm_drugs_master.ArrestinInteraction=='putative')]
    signalling_binding[["EntryName","Family","Ligandtype","Class","Uniprot","gene","GPCRdb","SequenceNumber","Segment","GPCRdbWT","NMaa","MutationType","Allele Count","Allele Frequency","Allele Number","Number of Homozygotes","sift_word","sift_score","polyphen_word","polyphen_score","foldchangeMaxAbs","diseaseId","score","diseaseName","Type","PTMsite","LB_structure","LB_fam","GProteinInteraction","ArrestinInteraction","ActivationPathway","MicroSwitch","SodiumPocket"]].to_excel(func_sites, 'Gprotein-Arrestin', index=False)

    PTMsites = dfm_drugs_master[(dfm_drugs_master.PTMsite=='yes')]
    PTMsites[["EntryName","Family","Ligandtype","Class","Uniprot","gene","GPCRdb","SequenceNumber","Segment","GPCRdbWT","NMaa","MutationType","Allele Count","Allele Frequency","Allele Number","Number of Homozygotes","sift_word","sift_score","polyphen_word","polyphen_score","foldchangeMaxAbs","diseaseId","score","diseaseName","Type","PTMsite","LB_structure","LB_fam","GProteinInteraction","ArrestinInteraction","ActivationPathway","MicroSwitch","SodiumPocket"]].to_excel(func_sites, 'PTM sites', index=False)

    switches = dfm_drugs_master[(dfm_drugs_master['MicroSwitch']=='yes') | (dfm_drugs_master['SodiumPocket']=='yes') | (dfm_drugs_master['ActivationPathway']=='yes')]
    switches[["EntryName","Family","Ligandtype","Class","Uniprot","gene","GPCRdb","SequenceNumber","Segment","GPCRdbWT","NMaa","MutationType","Allele Count","Allele Frequency","Allele Number","Number of Homozygotes","sift_word","sift_score","polyphen_word","polyphen_score","foldchangeMaxAbs","diseaseId","score","diseaseName","Type","PTMsite","LB_structure","LB_fam","GProteinInteraction","ArrestinInteraction","ActivationPathway","MicroSwitch","SodiumPocket"]].to_excel(func_sites, 'Microswitches and sodium pocket', index=False)
    func_sites.save()

    ### Mutations
    MAF_merge = dfm_drugs_master[(dfm_drugs_master.foldchangeMaxAbs>=0)].sort_values(by='foldchangeMaxAbs',ascending=False)[['EntryName','SequenceNumber','GPCRdb','WTaa','NMaa','Number of Homozygotes','Allele Count', 'sift_score', 'polyphen_score', 'foldchangeMaxAbs', 'foldchangeMaxNeg', 'foldchangeMaxPos', 'func_known', 'func_putative']]
    writer = ExcelWriter('~/Downloads/Supp5_in_vitro_mutations.xlsx')
    MAF_merge.to_excel(writer, index=False)
    writer.save()
