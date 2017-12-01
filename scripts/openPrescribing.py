"""
Created in Jan 2017

@author: Alexander Hauser <alexshauser@gmail.com>

Extract data from OpenPrescribing from the NHS England per Drug
API: https://openprescribing.net/api/
https://openprescribing.net/
All BNF sections: https://openprescribing.net/bnf/
"""

import pandas as pd
import requests
import json
import urllib2

import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import *

def query(query):

    url = "https://openprescribing.net/api/1.0/bnf_code/?format=json&q=" + query
    results = requests.get(url).json()

    return results

def spending(code):

    url = "https://openprescribing.net/api/1.0/spending/?format=json&date=2015-04-01&code=" + code
    spendingdata = requests.get(url).json()

    return spendingdata

def get_NHS_data():
    drugs = pd.read_csv("../data/drug_data.csv")
    drugnames = drugs[(drugs.Status=="approved")]['Drug Name'].unique()
    drugtargets = drugs[(drugs.Status=="approved")]['EntryName'].unique()

    # =====
    # DRUGS
    # =====
    NHS = pd.DataFrame(columns=['drugNameQuery','drugTarget','Subfamily','drugClass','drugName','Approval','section','drugCode','actual_cost', 'date', 'items', 'quantity'])
    for drug in tqdm(drugnames):
        results = query(drug)

        if not results:
            drug_alias = drugs[(drugs['Drug Name']==drug)]['DrugAliases'].values[0]
            if str(drug_alias) != 'nan':
                results = query(drug_alias)

        for result in results:
            if result['type']=='chemical':
                code = result['id']
                spendingdata = spending(code)
                temp = pd.DataFrame(spendingdata)

                for target in list(drugs[drugs['Drug Name']==drug]['EntryName'].unique()):
                    temp['drugCode'] = code
                    temp['section'] = result['section']
                    temp['drugName'] = result['name']
                    temp['drugNameQuery'] = drug
                    temp['drugTarget'] = target
                    temp['Subfamily'] = drugs[drugs['EntryName']==target]['Subfamily'].values[0]
                    temp['drugClass'] = drugs[drugs['Drug Name']==drug]['Drug Class'].values[0]
                    temp['Approval'] = drugs[drugs['Drug Name']==drug]['Approval'].values[0]
                    NHS = NHS.append(temp)

    missing_drugs = list(set(drugnames)-set(NHS.drugNameQuery.unique()))
    NHS.to_csv("../data/NHS_raw.csv")

def analysis():
    drugs = pd.read_csv("../data/drug_data.csv")
    drugnames = drugs[(drugs.Status=="approved")]['Drug Name'].unique()
    drugtargets = drugs[(drugs.Status=="approved")]['EntryName'].unique()

    ### ==================================
    NHS = pd.read_csv('../data/170829_NHS_raw.csv')
    NHS_section = pd.DataFrame({'actual_cost': NHS.groupby( ['section','date'] )['actual_cost'].sum()}).reset_index()
    NHS_section.to_csv("../data/drug_sales_section.csv")

    NHS_drugClass = pd.DataFrame({'actual_cost': NHS.groupby( ['drugClass','date'] )['actual_cost'].sum()}).reset_index()
    NHS_drugClass.to_csv("../data/drug_sales_drugClass.csv")

    relative_table = pd.read_csv("../data/170420_relative_table.csv")
    relative_table.rename(columns={'EntryName': 'drugTarget'}, inplace=True)
    relative_table = relative_table[['drugTarget','RelativeSNP','snp_count','seqLen']]

    # =======
    # TARGETS
    # =======
    # Loop over all receptors:
    NHS_targets = pd.DataFrame(columns=['drugTarget','Subfamily','actual_cost', 'date', 'items'])

    for target in drugtargets:
        subNHS = NHS[NHS['drugTarget'].str.contains(target)]
        by_month = pd.DataFrame({'actual_cost' : subNHS.groupby( ['date'] )['actual_cost'].sum(),
            'items' : subNHS.groupby( ['date'] )['items'].sum()}).reset_index()
        by_month['drugTarget'] = target
        by_month['Subfamily'] = drugs[drugs['EntryName']==target]['Subfamily'].values[0]

        NHS_targets = NHS_targets.append(by_month)
    NHS_targets.to_csv("../data/drug_sales_target_months.csv")

    # Over targets
    NHS_targets_agg = pd.DataFrame({'actual_cost_5y': NHS_targets.groupby( ['drugTarget','Subfamily'] )['actual_cost'].sum(),
        'items' : NHS_targets.groupby( ['drugTarget','Subfamily'] )['items'].mean()}).reset_index()
    NHS_targets_agg.sort_values(by='actual_cost_5y').tail()

    NHS_targets_agg = NHS_targets_agg.merge(relative_table, on=['drugTarget'])
    NHS_targets_agg.to_csv("../data/drug_sales_target.csv")

    # Over families
    NHS_targets_agg = pd.DataFrame({'actual_cost_5y': NHS_targets.groupby( ['Subfamily'] )['actual_cost'].sum()}).reset_index()

    # ========== LOF ===========
    lof = pd.read_csv("../data/170131_SNP_LoF_ExAC.csv")
    lof_target = pd.DataFrame({'positions' : lof.groupby( ['EntryName'] ).size(),
            'ac_hom' : lof.groupby( ['EntryName'] )['Number of Homozygotes'].sum(),
            'ac' : lof.groupby( ['EntryName'] )['Allele Count'].sum(),
            'min_ind' : lof.groupby( ['EntryName'] )['Allele Count'].max(),
            'max_ac' : lof.groupby( ['EntryName'] )['Allele Count'].max()}).reset_index()
    lof_target['max_ind'] = lof_target['ac']-lof_target['ac_hom']

    lof_target.rename(columns={'EntryName': 'drugTarget'}, inplace=True)
    NHS_agg = NHS.merge(relative_table, on=['drugTarget']).merge(lof_target, on=['drugTarget'])

    putative_functional = pd.read_csv("../data/dfm_drugs_masterTable.csv") #170824_dfm_masterTable.csv
    pfdf = pd.DataFrame({'n_mis_func' : putative_functional.groupby( ['EntryName'] )['SequenceNumber'].nunique()}).reset_index()
    NHS_agg = NHS_agg.merge(pfdf, left_on = 'drugTarget', right_on = 'EntryName')
    NHS_agg['RelativeSNP_func'] = NHS_agg.apply(lambda x: round(float(x['n_mis_func'])/x['seqLen'],2), axis=1)
    # NHS_agg.sort_values(by='RelativeSNP_func')
    NHS_agg = NHS_agg.dropna()

    DrugMeanDensity = pd.DataFrame({'meanDensity' : NHS_agg.groupby( ['drugNameQuery'] )['RelativeSNP'].mean(),
        'maxDensity' : NHS_agg.groupby( ['drugNameQuery'] )['RelativeSNP'].max(),
        'maxDensity_func' : NHS_agg.groupby( ['drugNameQuery'] )['RelativeSNP_func'].max(),
        'MV_count_max' : NHS_agg.groupby( ['drugNameQuery'] )['snp_count'].max(),
        'MV_count_max_func' : NHS_agg.groupby( ['drugNameQuery'] )['n_mis_func'].max(),
        'LoF_count_max' : NHS_agg.groupby( ['drugNameQuery'] )['max_ind'].max(),
        'LoF_count_min' : NHS_agg.groupby( ['drugNameQuery'] )['min_ind'].max(),
        'monthly_items' : NHS_agg.groupby( ['drugNameQuery'] )['items'].mean(),
        'monthly_cost' : NHS_agg.groupby( ['drugNameQuery'] )['actual_cost'].mean()}).reset_index()
    DrugMeanDensity['LoF_count_min_freq'] = DrugMeanDensity['LoF_count_min']/60760*100
    DrugMeanDensity.to_csv("../data/170829_drugs_sales.csv")

    # ========== CNV ===========
    cnv_target = pd.read_csv("../data/CNV/170216_cnv_GPCRs.csv")
    cnv_target.rename(columns={'EntryName': 'drugTarget'}, inplace=True)
    cnv_target['cnv'] = cnv_target['del'] + cnv_target['dup']

    NHS = pd.read_csv("../data/170829_NHS_raw.csv")
    relative_table.rename(columns={'EntryName': 'drugTarget'}, inplace=True)
    NHS_agg = NHS.merge(relative_table, on=['drugTarget']).merge(cnv_target, on=['drugTarget'])

    DrugMeanDensity = pd.DataFrame({'meanDensity' : NHS_agg.groupby( ['drugNameQuery'] )['RelativeSNP'].mean(),
        'maxcnv' : NHS_agg.groupby( ['drugNameQuery'] )['cnv'].max(),
        'monthly_items' : NHS_agg.groupby( ['drugNameQuery'] )['items'].mean(),
        'monthly_cost' : NHS_agg.groupby( ['drugNameQuery'] )['actual_cost'].mean()}).reset_index()
    DrugMeanDensity.to_csv("../data/drugs_sales_cnv.csv")

    #### GPCR SALES/ ITEMS compared to all in 2016:
    url = "https://openprescribing.net/api/1.0/spending/?format=json"
    all_spending = requests.get(url).json()
    all_spending_df = pd.read_json(json.dumps(all_spending))

    all_items_2016 = all_spending_df[(all_spending_df['date'] >= '2016-01-01') & (all_spending_df['date'] < '2017-01-01')]['items'].sum()
    all_costs_2016 = all_spending_df[(all_spending_df['date'] >= '2016-01-01') & (all_spending_df['date'] < '2017-01-01')]['actual_cost'].sum()

    GPCR_items_2016 = NHS[NHS['date'].str.contains('2016-')][['drugName','items','actual_cost']].drop_duplicates()['items'].sum()
    GPCR_costs_2016 = NHS[NHS['date'].str.contains('2016-')][['drugName','items','actual_cost','date']].drop_duplicates()['actual_cost'].sum()

    print round(float(GPCR_items_2016)/all_items_2016*100) # GPCR ITEMS % of all
    print round(float(GPCR_costs_2016)/all_costs_2016*100) # GPCR SALES % of all
