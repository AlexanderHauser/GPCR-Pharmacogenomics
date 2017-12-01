"""
Created in Dec 2016

@author: Alexander Hauser <alexshauser@gmail.com>

- List of functions for datasets that have been
employed thorughout this project
- Relevant APIs for mainly GPCRdb data retrieval
"""

import pandas as pd
import requests
import re, math, glob

from matplotlib import gridspec
import seaborn as sns
import matplotlib.pyplot as plt

from pandas import ExcelWriter

import time
import datetime
from tqdm import *
now = datetime.datetime.now()
timestamp = str(now.year)[2:]+['0' if now.month < 10 else ''][0]+str(now.month)+['0' if now.day < 10 else ''][0]+str(now.day)

### ================ APIs ================
def get_protein_residues(EntryName):
    ## get protein residue numbering

    url = "http://test.gpcrdb.org/services/residues/" + EntryName
    res_data = requests.get(url).json()

    return res_data

def get_protein_SequenceNumber(EntryName, SequenceNumber):

    url = "http://test.gpcrdb.org/services/protein/" + EntryName
    res_data = requests.get(url).json()
    seq = res_data['sequence'][int(SequenceNumber)-1]

    return seq

def get_protein_GPCRdb(EntryName,GPCRdb):
    url = "http://gpcrdb.org/services/alignment/protein/" + EntryName + '/' + GPCRdb + '/'
    res_data = requests.get(url).json()

    return res_data[EntryName]

def get_all_GPCRs():
    ## retrieve all GPCRs from GPCRdb

    url = "http://gpcrdb.org/services/proteinfamily/?page=000"
    res_data = requests.get(url).json()
    df = pd.read_json(json.dumps(res_data))
    df[df.slug.str.len()==15]
    df.drop(['parent'], inplace=True, axis=1)

def get_all_structures():
    ## retrieve all GPCR structures from GPCRdb

    url = "http://gpcrdb.org/services/structure/"
    ## LOCAL SERVER:
    # url = "http://0.0.0.0:8000/services/structure/"
    structures = requests.get(url).json()
    unique_structures = list(set([i['protein'] for i in structures]))

    return unique_structures

def get_alignments(proteins):
    ## get a sequence alignment with residue properties

    ## Just an input of receptors
    snp = pd.read_csv("../data/ExAC_GPCRs_raw.csv", low_memory=False)
    proteins = ','.join(snp[snp.Class=='Class A (Rhodopsin)'].EntryName.unique())

    ## Define segment order here:
    # segments_order = ['TM1', 'ICL1', 'TM2', 'ECL1', 'TM3', 'ICL2', 'TM4', 'ECL2', 'TM5', 'TM6', 'TM7', 'H8'] #
    segments_order = ['TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7']

    ## Define segment range (start,end of segment) here:
    segment_ranges = {'TM1':['1x30','1x60'], 'TM2':['2x37','2x66'], 'TM3':['3x21','3x56'], 'TM4':['4x39','4x63'], 'TM5':['5x36','5x68'], 'TM6':['6x30','6x61'], 'TM7':['7x30','7x56']}

    positions = {}
    for segment in segment_ranges.keys():
        start = int(segment_ranges[segment][0].split('x')[1])
        end = int(segment_ranges[segment][1].split('x')[1])
        positions[segment] = [segment.split('TM')[1]+'x'+str(i) for i in range(start, end + 1)]

    ## remove bulge from class A TM7
    positions['TM7'].remove('7x44')

    GP = pd.DataFrame(columns=[])

    ## needs to be careful with proteinsets and their segment lengths!!
    for SEG in segments_order:
        GPCRdb_num = positions[SEG]
        if len(GPCRdb_num) > 0:
            url = "http://test.gpcrdb.org/services/alignment/protein/" + proteins + "/" + ','.join(GPCRdb_num) + "/statistics/"

            alignment_data = requests.get(url).json()

            for protein in alignment_data:

                if protein != "statistics":
                    for i, num in enumerate(GPCRdb_num):
                        try: # in case there is no aa at this position
                            GP.loc[protein, num] = alignment_data[protein][i]
                        except:
                            GP.loc[protein, num] = '-'

                elif protein == "statistics":
                    for prop in alignment_data[protein]:
                        for i, num in enumerate(GPCRdb_num):
                            try:
                                GP.loc[prop, num] = alignment_data[protein][prop][i]
                            except:
                                GP.loc[prop, num] = 0

    identity = GP[GP.index.isin(['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'])]
    identity.loc['MaxIdentity'] = identity.astype(int).max()
    # max_identity_values = list(identity.astype(int).max())
    # identity['MaxIdentity':].to_dict()

    fasta = "../msa.fasta"
    with open(fasta, "a") as fname:
        for protein in proteins.split(','):
            seq = ''.join(list(GP.loc[protein]))
            fname.write(">" + protein + "\n" + seq + "\n")

def eva_api():
    ## API for EVA 1000 Genomes - not used

    url = 'https://www.ebi.ac.uk/eva/webservices/rest/v1/segments/19:35358460-35360485/variants?species=hsapiens_grch37&studies=PRJEB6930'
    data = requests.get(url).json()

    for i in data['response']:
        print i

def ensemble_API_variants(hgvsc):
    ## API for ENSEMBLE hbgsc - not used

    url = "http://rest.ensembl.org/vep/human/hgvs/" + hgvsc + "?content-type=application/json"
    variant_data = requests.get(url).json()

    return variant_data

    return
### ================ END ================

### ============== DATASETS =============
def GPCR_res_table():
    ## Makes a all GPCR residue Table

    ## Cany input of protein entry names possible here:
    relativeTable = pd.read_csv("../data/relative_table.csv")
    table = pd.DataFrame(columns = ['EntryName', 'Family', u'amino_acid', u'display_generic_number', u'protein_segment',
       u'sequence_number'])
    for EntryName in tqdm(relativeTable.EntryName.unique()):
        url = "http://www.gpcrdb.org/services/residues/" + EntryName
        res_data = requests.get(url).json()
        df = pd.read_json(json.dumps(res_data))
        df['EntryName'] = EntryName
        df['Class'] = relativeTable[relativeTable.EntryName==EntryName].Class.values[0]
        df['Family'] = relativeTable[relativeTable.EntryName==EntryName].Family.values[0]
        table = pd.concat([table,df], ignore_index=True)
    table.rename(columns={'display_generic_number': 'GPCRdb', 'amino_acid': 'WTaa', 'protein_segment': 'Segment', 'sequence_number':"SequenceNumber"}, inplace=True)
    table.to_csv("../data/all_GPCR_residues.csv")

def receptor_sequence_length_table():
    ## Makes a dataframe of receptor sequence lengths

    df = pd.read_csv("../data/EST_ENSG_GENE.csv", low_memory=False)

    GPCRs = pd.DataFrame({'count' : df.groupby( ['EntryName','Uniprot'] ).size()}).reset_index()
    GPCRs['seqLen'] = 0

    for index, res in enumerate(GPCRs.iterrows()):

        seqLen = 0
        accession = GPCRs[index:index+1]['Uniprot'].values[0]

        url = "http://test.gpcrdb.org:80/services/protein/accession/" + accession
        data = requests.get(url).json()

        SeqLength = len(data['sequence'])
        GPCRs.loc[index,'seqLen'] = SeqLength

    GPCRs.drop(['count'],inplace=True, axis=1)
    GPCRs.to_csv("../data/"+timestamp+"_GPCR_lengths.csv")

def receptor_segment_lengths(receptor=True):
    ## Makes a dataframe of receptor segment lengths

    df = pd.read_csv("../data/ExAC_GPCRs_raw.csv", low_memory=False)

    if receptor:
        df = df[df['Class']!='G-Protein']
        protein_type = "GPCR"
    else:
        df = df[df['Class']=='G-Protein']
        protein_type = "Gprotein"

    Proteins = pd.DataFrame({'count' : df.groupby( ['EntryName','Uniprot','Family','Ligandtype','Class'] ).size()}).reset_index()

    for index, protein in tqdm(enumerate(Proteins.iterrows())):
        entry_name = Proteins[index:index+1]['EntryName'].values[0]

        url = "http://test.gpcrdb.org/services/residues/" + entry_name
        res_data = requests.get(url).json()
        # segments = {'TM1':0,'TM2':0,'TM3':0,'TM4':0,'TM5':0,'TM6':0,'TM7':0,'N-term':0,'ECL1':0,'ECL2':0,'ECL3':0,'C-term':0,'H8':0,'ICL1':0,'ICL2':0,'ICL3':0,'ICL4':0}
        segments = {}
        if index == 0:
            agg = pd.DataFrame()

        for res in res_data:
            if not res['protein_segment'] in segments.keys():
                segments[res['protein_segment']] = 1
            else:
                segments[res['protein_segment']] += 1

        new_df = pd.DataFrame(segments, index=[index])
        agg = agg.append(new_df)
        agg.is_copy = False

    Proteins = Proteins.merge(agg, left_index=True, right_index=True, how='outer')

    Proteins.to_csv("../data/"+timestamp+"_"+protein_type+"_segment_lengths.csv")

def receptor_ligand_interactions():
    ## get all structuree receptor-liagnd interactions

    # Incooporate different levels of interaction depending on the availability of data (receptor, familylevel, ligandtype ...)
    snp = pd.read_csv("../data/ExAC_GPCRs_raw.csv", low_memory=False)
    # GPCRs = pd.DataFrame({'count' : snp.groupby( ['EntryName','Uniprot','Family','Ligandtype','Class'] ).size()}).reset_index()
    GPCRs = pd.read_csv("../data/GPCR_lengths.csv")
    snp['ProteinID'] = snp.EntryName.str.split('_', expand=True)[0].str.upper()

    interactions = pd.DataFrame(columns=['EntryName','ProteinID','Uniprot','Family','Ligandtype','Class'])
    unique_structures = get_all_structures()
    unique_structures = [i for i in unique_structures if '_hcmva' not in i and 'todpa' not in i]

    # for index, protein in tqdm(enumerate(GPCRs.iterrows())):
    for entry_name in tqdm(unique_structures):
        # entry_name = GPCRs[index:index+1]['EntryName'].values[0]
        # Uniprot = GPCRs[index:index+1]['Uniprot'].values[0]
        # Family = GPCRs[index:index+1]['Family'].values[0]
        # Ligandtype = GPCRs[index:index+1]['Ligandtype'].values[0]
        # Class = GPCRs[index:index+1]['Class'].values[0]
        ProteinID = entry_name.split('_')[0].upper()

        # Take from Family, as there are non-human crystallized receptors:
        # http://www.gpcrdb.org/family/001_002_021_001/
        url = "http://www.gpcrdb.org:80/interaction/ajax/" + entry_name
        ## LOCAL SERVER:
        # url = "http://0.0.0.0:8000/interaction/ajax/" + entry_name
        try:
            data = requests.get(url).json()
            protein_residue_info = get_protein_residues(entry_name)

            Uniprot = snp[snp.ProteinID==ProteinID].Uniprot.values[0]
            Family = snp[snp.ProteinID==ProteinID].Family.values[0]
            # Ligandtype = snp[snp.ProteinID==ProteinID].Ligandtype.values[0]
            Class = snp[snp.ProteinID==ProteinID].Class.values[0]
        except:
            print url, entry_name

        for residue in data:
            pos = interactions.shape[0]+1

            interactions.loc[pos, 'GPCRdb'] = protein_residue_info[int(residue)-1]['display_generic_number']
            interactions.loc[pos, 'Residue'] = residue
            interactions.loc[pos, 'Interaction'] = list(set([i[1] for i in data[residue]]))
            interactions.loc[pos, 'EntryName'] = entry_name
            interactions.loc[pos, 'ProteinID'] = ProteinID

            interactions.loc[pos, 'Uniprot'] = Uniprot
            interactions.loc[pos, 'Family'] = Family
            # interactions.loc[pos, 'Ligandtype'] = Ligandtype
            interactions.loc[pos, 'Class'] = Class
    """
    Receptors in GPCRdb missing ligand interactions: 'ox1r_human', 'opsd_human', 'acm1_human', 'calcr_human', 'aa1r_human'
    """
    interactions.to_csv("../data/"+timestamp+"_ligand_interactions.csv")

def receptor_invitro_mutations():
    ## in vitro pharmacological mutations

    # url = "http://test.gpcrdb.org:80/interaction/ajax/" + accession
    # data = requests.get(url).json()
    # extract data directly from the SQL table: mutation_experiment, which includes precalculated foldchanges etc.
    # mutations = MutationExperiment.objects.filter(protein__entry_name=slug).order_by('residue__sequence_number').prefetch_related('residue')
    # get the right links and export as csv to facilitate better data anlysis
    # extract all 'positions' (per receptor) to have a pharmacological effect (where, how many?)
    # mutations_raw = pd.read_csv("../data/GPCRdb_mutation_raw.csv")
    mutations = pd.read_csv("../data/170125_GPCRdb_mutation.csv")
    ## This way its easier to interpret: high positive change means increased binding/potency, negative means decreased binding
    mutations['foldchange'] = mutations['foldchange']*-1
    mutations['proteinID'] = mutations['EntryName'].str.split('_', expand=True)[0]
    ## ONLY HUMAN RECEPTORS:
    # mutations = mutations[mutations.EntryName.str.contains('_human')]
    ## ALL RECEPTORS IF THEY ARE EQUIVALENT AT THE HUMAN RECEPTOR:
    for index, protein in tqdm(enumerate(mutations.iterrows())):
        EntryName = mutations[index:index+1]['EntryName'].values[0]
        WTaa =  mutations[index:index+1]['WTaa'].values[0]
        SequenceNumber = mutations[index:index+1]['SequenceNumber'].values[0]
        GPCRdb = mutations[index:index+1]['GPCRdb'].values[0]
        if not 'human' in EntryName:
            if GPCRdb and str(GPCRdb) != 'nan':
                humanWT = get_protein_GPCRdb(EntryName, GPCRdb.split('.')[0]+'x'+GPCRdb.split('x')[1])

            else:
                # Also look for same sequence positions for mutations outsite GPCRdb numbering!
                humanWT = get_protein_SequenceNumber(EntryName, SequenceNumber)

            if humanWT == WTaa:
                mutations.loc[index, 'SequenceNumber'] = SequenceNumber
                mutations.loc[index, 'EntryName'] = EntryName.split('_')[0]+'_human'
    mutations[mutations.EntryName.str.contains('human')]

    # ### Differences between orthologs:
    # # a = pd.DataFrame({'WTpos' : mutations.groupby( ['protein','GPCRdb'] )['WTaa'].nunique()}).reset_index()
    # indices = []
    # def func(sub_df):
    #     ## if there is no human in the invitro data to compare with, take the snp data and search for the ortholog
    #
    #     subdfHuman = sub_df[sub_df.EntryName.str.contains('human')]['WTaa'].unique()
    #     EntryName = sub_df.EntryName.values[0]
    #     SequenceNumber = sub_df.SequenceNumber.values[0]
    #     GPCRdb = sub_df.GPCRdb.values[0]
    #
    #     # if len(subdfHuman) > 0:
    #         # humanWT = sub_df[sub_df.EntryName.str.contains('human')]['WTaa'].unique()[0]
    #     if GPCRdb:
    #         humanWT = get_protein_GPCRdb(EntryName, GPCRdb.split('.')[0]+'x'+GPCRdb.split('x')[1])
    #     else:
    #         # Also look for same sequence positions for mutations outsite GPCRdb numbering!
    #         humanWT = get_protein_SequenceNumber(EntryName, SequenceNumber)
    #
    #     returndf = sub_df[sub_df['WTaa']==humanWT]
    #     indices.append(returndf.index)
    #
    # ## where the WTaa is the same in the ortholog species for the same GPCRdb number
    # mutations[~mutations.EntryName.str.contains('human')].groupby(['EntryName','SequenceNumber']).apply(lambda x: func(x))
    # indices2 = []
    # for i in np.array(indices).flatten():
    #     for n in i:
    #         indices2.append(n)
    # mutations = mutations[mutations.index.isin(indices2)]
    # mutations['EntryName'] = mutations['EntryName'].str.split('_', expand=True)[0]+'_human'
    # mutations['SequenceNumber'] = mutations['SequenceNumber'].astype(int)

    mutations.to_csv("../data/mutants_orthologs.csv")
    # ===========================================================================
    mutations = pd.read_csv("../data/mutants_orthologs.csv")

    df_foldchange = pd.DataFrame({'foldchange_min' : df.groupby( ['EntryName','Mutantaa'] )['foldchange'].min(),
        'foldchange_max' : df.groupby( ['GPCRdb','Mutantaa'] )['foldchange'].max(),'Protein_count' : df.groupby( ['GPCRdb','Mutantaa'] )['EntryName'].nunique(),
        'ligand_count' : df.groupby( ['GPCRdb','Mutantaa'] )['Ligand'].nunique()},
        'foldchangeAbs_max' : df.groupby( ['EntryName','Mutantaa'] )['foldchange_abs'].max()).reset_index()

    ### Difference in genetic variation of two sets of increase/decrease reported sets?
    # foldthreshold = 5
    # decreased = mutations[mutations.foldchange<-foldthreshold].GPCRdb.unique()
    # increased = mutations[mutations.foldchange>foldthreshold].GPCRdb.unique()
    # decreased_specific = list(set(decreased)-set(increased))
    # increased_specific = list(set(increased)-set(decreased))
    # print snp[(snp.GPCRdb.isin(increased_specific))].shape[0]/len(increased)
    # print snp[(snp.GPCRdb.isin(decreased_specific))].shape[0]/len(decreased)

    # Vis of subset of mutants with high range of effect on receptor ligand coupling
    # positions = mutations[mutations.foldchange<-foldthreshold].GPCRdb.value_counts().head(n=20).keys()
    # sub = mutations[mutations.GPCRdb.isin(positions)]
    # sns.boxplot(x="GPCRdb", y="foldchange", data=sub, palette="Greens_r")
    # plt.ylim(-100,100)
    # plt.xticks(rotation=90)
    # plt.show()

    ## SNP Pos Set
    def get_sublist_snp(EntryName, Residue):
        values = snp[(snp['EntryName']==EntryName) & (snp['SequenceNumber']==Residue)].index.tolist()
        return values

    relativeTable_invitro = pd.DataFrame(columns=['EntryName','Region','RelativeSNP'])

    # folds = ['decrease','increase']
    folds = ['all']
    foldthreshold = 1.5
    for fold in folds:
        # if fold == 'decrease':
        #     submutations = mutations[mutations.foldchange<-foldthreshold]
        # else:
        #     submutations = mutations[mutations.foldchange>foldthreshold]
        submutations = mutations[abs(mutations.foldchange)>foldthreshold]

        indices_raw_snp = map(get_sublist_snp, submutations['EntryName'],submutations['SequenceNumber'])
        indices_snp = [item for sublist in indices_raw_snp for item in sublist]
        snp_invitro = snp.ix[indices_snp]

        for index, receptor in enumerate(submutations[submutations['EntryName'].str.contains('_')]['EntryName'].unique()):
            number_of_invitromutations = len(submutations[submutations['EntryName']==receptor].GPCRdb.unique())
            number_of_snp_interactions = len(snp_invitro[snp_invitro['EntryName']==receptor].GPCRdb.unique())
            # number_of_cancer_interactions = len(cancer_interactions[cancer_interactions['EntryName']==receptor])
            if len(snp[snp['EntryName']==receptor])>=1:
                relativeTable_invitro.loc[index, 'EntryName'] = receptor
                relativeTable_invitro.loc[index, 'Region'] = 'InvitroMutation' # _'+fold
                relativeTable_invitro.loc[index, 'NumberOfSNPs'] = number_of_snp_interactions
                relativeTable_invitro.loc[index, 'NumberOfinvitroMutants'] = number_of_invitromutations
                relativeTable_invitro.loc[index, 'RelativeSNP'] = float(number_of_snp_interactions)/number_of_invitromutations
                # polyphen_mean = snp_interactions[snp_interactions['EntryName']==receptor].polyphen_score.mean()
                # sift_score2_mean = snp_interactions[snp_interactions['EntryName']==receptor].sift_score2.mean()
                # relativeTable_invitro.loc[index, 'sift_score2_mean'] = sift_score2_mean
                # relativeTable_invitro.loc[index, 'polyphen_mean'] = polyphen_mean

    relativeTable_invitro = relativeTable_invitro[relativeTable_invitro.RelativeSNP!='NaN']
    relativeTable_invitro = relativeTable_invitro[relativeTable_invitro.NumberOfinvitroMutants>=5]
    relativeTable_invitro['RelativeSNP'] = relativeTable_invitro['RelativeSNP'].astype(float)
    relativeTable_invitro.to_csv("../data/relativeTable_InVitroMutants.csv")
    # ==========================================================================
    ### How many SNPs have already been experimentally 'validated'?
    drugs_raw = pd.read_excel("../data/drug_data.csv")
    drugTargets = drugs_raw[drugs_raw['Status']=='approved']['EntryName'].unique()

    # mutations[(mutations.EntryName == 'sctr_human') & (mutations.SequenceNumber==66) & (mutations.WTaa=='C') & (mutations.Mutantaa == 'S')]
    dfm_drugs_master = pd.read_csv("../data/dfm_drugs_masterTable.csv")
    dfm_drugs_master.rename(columns={'NMaa': 'Mutantaa'}, inplace=True)
    invitroSNPs = mutations.merge(dfm_drugs_master, on=['EntryName','SequenceNumber','Mutantaa','WTaa'])
    invitroSNPs = invitroSNPs[invitroSNPs.columns[~invitroSNPs.columns.str.contains('Unnamed:')]]
    invitroSNPs.drop_duplicates(inplace=True)
    ## Make an aggegration over positions and ommit different experiments/ligands; mean over the folfchange at the position
    invitroSNPs['foldchange_abs'] = abs(invitroSNPs['foldchange'])
    invitroSNPs = invitroSNPs[invitroSNPs.EntryName.isin(drugTargets)]
    invitroSNPsdf = pd.DataFrame({'foldchange_max' : invitroSNPs.groupby(['EntryName','SequenceNumber','Mutantaa'])['foldchange_abs'].max(), 'foldchangePos_max' : invitroSNPs.groupby(['EntryName','SequenceNumber','Mutantaa'])['foldchange'].max()}).reset_index()
    # invitroSNPsdf[abs(invitroSNPsdf.foldchange_max)>1.5].shape

    ## WHat are the MAF at those positions?
    # MAF_merge = invitroSNPsdf.merge(snp[['EntryName','SequenceNumber','Mutantaa','WTaa']], on=['EntryName','SequenceNumber'])
    MAF_merge = invitroSNPsdf.merge(dfm_drugs_master[['EntryName','SequenceNumber','GPCRdb','WTaa','Mutantaa','Number of Homozygotes','Allele Count', 'sift_score', 'polyphen_score', 'func_known', 'func_putative']], on=['EntryName','SequenceNumber','Mutantaa'])
    MAF_merge[MAF_merge.func_known==True].shape
    MAF_merge[(MAF_merge.func_putative==True)].shape
    MAF_merge[(MAF_merge.func_putative==True) & (MAF_merge.foldchange_max>=5)]
    MAF_merge[MAF_merge.foldchange_max>=5].shape # positive change in efficacy --> more ADR?

    ## Export
    from pandas import ExcelWriter
    writer = ExcelWriter('../SuppTable_in_vitro_mutations.xlsx')
    MAF_merge.to_excel(writer, index=False)
    writer.save()

def table_view():
    ## Table overview of MV, LOF and CNV mutations
    ## HIERACHICAL CLUSTERING output to R

    drugs_raw = pd.read_excel("../data/drug_data.csv")
    dfm_drugs_master = pd.read_csv("../data/dfm_drugs_masterTable.csv")


    drugTargets = drugs_raw[(drugs_raw['Status']=='approved') & (drugs_raw.PMID.notnull())]['EntryName'].unique()
    trialTargets = np.array(list(set(drugs_raw['EntryName'].unique()) - set(drugs_raw[drugs_raw['Status']=='approved']['EntryName'].unique())))

    relative_table_all = pd.read_csv("../data/relative_table.csv")
    relative_table = relative_table_all[relative_table_all.EntryName.isin(drugTargets)]
    relative_tableT = relative_table_all[relative_table_all.EntryName.isin(trialTargets)]

    # Update
    exac_GPCR_raw = pd.read_csv("../data/170216_exac_GPCRs.csv") # exac_GPCRs
    exac_GPCR = exac_GPCR_raw[exac_GPCR_raw.EntryName.isin(drugTargets)]
    exac_GPCRT = exac_GPCR_raw[exac_GPCR_raw.EntryName.isin(trialTargets)]

    cnf_gpcr_raw = pd.read_csv("../data/CNV/170216_cnv_GPCRs.csv")
    cnf_gpcr = cnf_gpcr_raw[cnf_gpcr_raw.EntryName.isin(drugTargets)]
    cnf_gpcrT = cnf_gpcr_raw[cnf_gpcr_raw.EntryName.isin(trialTargets)]

    lof = pd.read_csv("../data/170131_SNP_LoF_ExAC_ENSG.csv")
    lof['ac_het'] = lof['Allele Count']-lof['Number of Homozygotes']
    lof_target_raw = pd.DataFrame({'LoF_positions' : lof.groupby( ['EntryName'] ).size(),
            'ac_hom' : lof.groupby( ['EntryName'] )['Number of Homozygotes'].sum(),
            'ac' : lof.groupby( ['EntryName'] )['Allele Count'].sum(),
            'min_ind' : lof.groupby( ['EntryName'] )['ac_het'].max(),
            'max_ac' : lof.groupby( ['EntryName'] )['Allele Count'].max()}).reset_index()
    lof_target_raw['max_ind'] = lof_target_raw['ac']-lof_target_raw['ac_hom']
    lof_target_raw['LoF_count_min_freq'] = lof_target_raw['min_ind']/60760*100
    lof_target_raw['LoF_count_max_freq'] = lof_target_raw['max_ind']/60760*100
    lof_target = lof_target_raw[lof_target_raw.EntryName.isin(drugTargets)]
    lof_targetT = lof_target_raw[lof_target_raw.EntryName.isin(trialTargets)]

    ### EXPORT AS TABLE / Supplementary!
    cols1 = ['EntryName','Family', 'Class', 'seqLen',
       'snp_count', 'RelativeSNP', 'LoF_positions', 'LoF_count_min_freq', 'del', 'dup']
    Master = relative_table.merge(lof_target, on='EntryName', how='left').merge(cnf_gpcr, on='EntryName', how='left')
    # Master = Master[Master.EntryName.isin(drugs_raw['EntryName'].unique())][cols1]
    Master = Master[Master.EntryName.isin(drugTargets)]
    MasterT = relative_tableT.merge(lof_targetT, on='EntryName', how='left').merge(cnf_gpcrT, on='EntryName', how='left')
    MasterT = MasterT[MasterT.EntryName.isin(drugs_raw['EntryName'].unique())][cols1]
    MasterT = MasterT[MasterT.EntryName.isin(trialTargets)]

    ### HIERACHICAL CLUSTERING ON Z-Scores
    colsZS = ['snp_count','RelativeSNP','LoF_positions','LoF_count_min_freq','del','dup']
    for col in colsZS:
        col_zscore = col + '_zscore'
        Master[col_zscore] = (Master[col] - Master[col].mean())/Master[col].std(ddof=0)
        MasterT[col_zscore] = (MasterT[col] - MasterT[col].mean())/MasterT[col].std(ddof=0)
    zz = ["snp_count_zscore", "RelativeSNP_zscore", "LoF_positions_zscore", "LoF_count_min_freq", "del_zscore", "dup_zscore"]

    ## EXPORT TO R
    Master.index = Master.EntryName
    MasterT.index = MasterT.EntryName
    exp = Master[zz]
    expT = MasterT[zz]
    Master = Master.fillna(0)
    MasterT = MasterT.fillna(0)
    Master['mean zscores'] = Master[zz].apply(lambda x: x.values.mean(), axis=1)
    MasterT['mean zscores'] = MasterT[zz].apply(lambda x: x.values.mean(), axis=1)
    Master[zz].to_csv("../data/zscores.csv")

    Master.sort_values(by='mean zscores').tail(n=10)
    MasterT.sort_values(by='mean zscores').tail(n=10)
    Master['CNVs'] = Master['del'] + Master['dup']

    ## EXCEL Export:
    for z in zz:
        Master[z].replace(0, 'NA',inplace=True)
        MasterT[z].replace(0, 'NA',inplace=True)
    Master.rename(columns={"snp_count_zscore": 'MV positions count (z-score)', "RelativeSNP_zscore": 'MV density (z-score)', "LoF_positions_zscore": 'LoF positions count (z-score)', "LoF_count_min_freq": 'LoF min pop freq (z-score)', "del_zscore": 'number of observed deletions (z-score)', "LoF_count_min_freq_zscore": "minimum pop frequency (z-score)", "dup_zscore": 'number of observed duplications (z-score)', 'snp_count': 'MV_positions', 'del': 'deletions', 'dup': 'duplications', 'LoF_count_min_freq': 'Lof min pop freq', 'RelativeSNP': 'MV_density'}, inplace=True)
    MasterT.rename(columns={"snp_count_zscore": 'MV positions count (z-score)', "RelativeSNP_zscore": 'MV density (z-score)', "LoF_positions_zscore": 'LoF positions count (z-score)', "LoF_count_min_freq": 'LoF min pop freq (z-score)', "del_zscore": 'number of observed deletions (z-score)', "LoF_count_min_freq_zscore": "minimum pop frequency (z-score)", "dup_zscore": 'number of observed duplications (z-score)', 'snp_count': 'MV_positions', 'del': 'deletions', 'dup': 'duplications', 'LoF_count_min_freq': 'Lof min pop freq', 'RelativeSNP': 'MV_density'}, inplace=True)
    writer = ExcelWriter('../SuppTable_drug_targets_zscores.xlsx')
    Master.sort_values(by='mean zscores', ascending=False).to_excel(writer,sheet_name='GPCR drug targets (n=108)', index=False)
    MasterT.sort_values(by='mean zscores', ascending=False).to_excel(writer,sheet_name='GPCR trial targets (n=66)', index=False)
    writer.save()
    # Master = Master.merge(exp, on = 'EntryName')
    # MasterT = MasterT.merge(exp, on = 'EntryName')

    ## VIS
    plt.style.use('seaborn-white')
    fig = plt.figure(figsize=(12, 6))
    gridtranslate = {"n_syn": '#', "syn_z": "z-score", "n_mis": '#', "mis_z": "z-score", "RelativeSNP": "density", "snp_rv": 'rare', "snp_cv": 'common', "lof_z": "z-score", 'del': "#deletions", "dup": "#duplications", "LoF_positions": 'LoF_n_pos', 'LoF_count_min_freq': 'pop. min. freq.'}
    colors = ["#ed1c24","#ed1c24","#5992b0","#5992b0","#8d54a2","#8d54a2"]
    # gridspace = ["n_syn", "syn_z", "n_mis", "mis_z", "RelativeSNP", "snp_rv", "snp_cv","n_lof", "lof_z", 'del', 'dup']
    gridspace = ["n_mis", "RelativeSNP", "LoF_positions", 'LoF_count_min_freq', 'del', 'dup']
    gs = gridspec.GridSpec(1, len(gridspace), width_ratios=[1]*len(gridspace))
    for i, grid in enumerate(gridspace):
        ax = plt.subplot(gs[i])
        ax.title.set_text(gridtranslate[grid])
        # ax0.set_suptitle("Loss of Function mutations", fontsize=16, fontweight="bold")
        if (grid == 'del' or grid == 'dup'):
            cnf_gpcr = cnf_gpcr.sort_values(by=[grid])
            topcnf = cnf_gpcr[cnf_gpcr[grid]>0].head(n=10).append(cnf_gpcr.tail(n=10))
            axbar = ax.barh(range(0,20), topcnf[grid].values, alpha=0.8, color=colors[i])
            labels = topcnf.ProteinID.values
            max_val = round(topcnf[grid].values.max())
            # min_val = round(topcnf[grid].values.min())
            min_val = 0
        elif grid == "RelativeSNP" or grid == "snp_rv" or grid == "snp_cv":
            relative_table = relative_table.sort_values(by=[grid])
            toprel = relative_table.head(n=10).append(relative_table.tail(n=10))
            axbar = ax.barh(range(0,20), toprel[grid].values, alpha=0.8, color=colors[i])
            # labels = toprel["Gene names  (primary )"].values
            labels = [i.split('_human')[0].upper() for i in toprel["EntryName"].values]
            max_val = toprel[grid].values.max()
            min_val = 0
        elif grid == 'LoF_positions' or grid == 'LoF_count_min_freq':
            lof_target = lof_target.sort_values(by=[grid])
            toprel = lof_target[lof_target[grid]>0].head(n=10).append(lof_target.tail(n=10))
            axbar = ax.barh(range(0,20), toprel[grid].values, alpha=0.8, color=colors[i])
            # labels = toprel["Gene names  (primary )"].values
            labels = toprel["EntryName"].str.split('_',expand=True)[0].str.upper().values
            max_val = round(toprel[grid].values.max(),0)
            # min_val = toprel[grid].values.min()
            min_val = round(toprel[grid].values.min(),0)
        else:
            exac_GPCR = exac_GPCR.sort_values(by=[grid])
            topset = exac_GPCR.head(n=10).append(exac_GPCR.tail(n=10))
            axbar = ax.barh(range(0,20), topset[grid].values, alpha=0.8, color=colors[i])
            labels = topset.ProteinID.str.upper().values
            max_val = round(topset[grid].values.max(),0)
            # min_val = round(topset[grid].values.min(),0)
            min_val = 0

        ax.set_xticks([min_val, max_val])
        ax.set_xlim([min_val,max_val])
        ax.set_ylim([-0.5,19.5])
        ax.spines['bottom'].set_color('#EDEDED')
        ax.spines['top'].set_color('#EDEDED')
        ax.spines['right'].set_color('#EDEDED')
        ax.spines['left'].set_color('#EDEDED')
        # ax.set_xticklabels(labels, rotation=40)
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=60)
        rects = axbar.patches
        plt.setp(ax.get_yticklabels(), visible=False)
        for rect, label in zip(rects, labels):
            height = rect.get_height()
            ax.text(rect.get_x(), rect.get_y()+0.05, label, fontsize=10, family='sans-serif', style='italic', ha='left',va='bottom')
    fig.tight_layout()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.075, hspace=0)
    plt.savefig("../figures/table1.svg", bbox_inches='tight')
    plt.show()

def ExAC_GPCR_data():
    ## EXAC data refinement
    ## update with SequenceNumbers
    ## follows relative_table()

    ## Raw file from EXAC ressource webserver
    snp = pd.read_csv("../data/170216_exac_GPCRs.csv")
    snp = snp[snp['Class']!='G-Protein']

    ## GPCR IDs
    enst = pd.read_csv("../data/ENST_Protein.csv")
    ENSTs = enst.ENST.tolist()
    uni_ensemble = pd.read_csv("../data/EST_ENSG_GENE.csv")
    searchfor = uni_ensemble.ENSG.unique()

    ## EXAC RAW file (~8gb) from publication
    exac = pd.read_csv('../data/exac.csv', dtype=types)

    ## MV only
    exac_raw = exac[(exac.category=='missense_variant')] # discripencay between exac online (with use == False) and the R script (github), which excludes use = False
    exac_raw_GPCRs = exac_raw[exac_raw.gene.str.contains('|'.join(searchfor))]
    exac_raw_GPCRs = exac_raw_GPCRs.reset_index()

    ## Split amino_acids column
    exac_raw_GPCRs['WTaa'] = exac_raw_GPCRs.amino_acids.str.split('/', expand=True)[0]
    exac_raw_GPCRs['NMaa'] = exac_raw_GPCRs.amino_acids.str.split('/', expand=True)[1].str.split(',', expand=True)[0]
    ## Clean multiple ENSG cells e.g.:
    ## exac_raw_GPCRs[exac_raw_GPCRs.gene=='ENSG00000268170,ENSG00000221937'] index = 32365
    exac_raw_GPCRs['gene'] = exac_raw_GPCRs.apply(lambda x: ','.join([gene for gene in x['gene'].split(',') if gene in searchfor]), axis=1)

    ## Clean multiple ENSTs for known ENSTs
    exac_raw_GPCRs['hgvsc_refine'] = exac_raw_GPCRs.apply(lambda x: ','.join([hgv for hgv in x['hgvsc'].split(',') if hgv.split('.')[0] in ENSTs]), axis=1)
    # exac_raw_GPCRs[exac_raw_GPCRs.hgvsc_refine.str.contains(',')]
    # exac_raw_GPCRs[exac_raw_GPCRs.hgvsc_refine != '')]

    ## compare GPCRdb with each position and only keep if match with correct AA
    ## loop over table --> match with snp and pull data
    ## in case of multiple sequence Numbers --> take the best match only
    ## Manually try to find the best ENST matches for GPCRdb WT proteins
    ## TESTING:
    # snp.rename(columns={'Sequence Number': 'SequenceNumber'}, inplace=True)
    # a = exac_raw_GPCRs[exac_raw_GPCRs.gene == 'ENSG00000182885'][['feature','WTaa', 'SequenceNumber']]
    # # a['SequenceNumber'] = a['SequenceNumber'].astype(int)
    # b = snp[snp.Uniprot=='Q86Y34'][['WTaa', 'SequenceNumber']].sort_values(by='SequenceNumber')
    # b.merge(a, how='left', on=['WTaa', 'SequenceNumber'])
    # a[a.feature.str.contains('ENST00000337539')]

    ## Map Protein Sequence Numbers onto exac_RAW:
    ## FOR EACH ROW CHECK FOR THE CORRECT FEATURE
    exac_raw_GPCRs = exac_raw_GPCRs[exac_raw_GPCRs.hgvsc_refine != '']
    exac_raw_GPCRs = exac_raw_GPCRs.reset_index()

    for index, rows in tqdm(enumerate(exac_raw_GPCRs.iterrows())):
        hgvsc = exac_raw_GPCRs[index:index+1]['hgvsc_refine'].values[0]
        sequenceNumbers = []
        for feat in hgvsc.split(','):
            cDNAposition = int(re.findall(r'\d+',feat.split(':c.')[1])[0])
            # For the same run devide the hgvsc by 3 to get from cDNA to protein sequence site
            sequenceNumber = str(int((math.ceil(cDNAposition/3.0))))
            if sequenceNumber not in sequenceNumbers:
                sequenceNumbers.append(str(sequenceNumber))
        exac_raw_GPCRs.loc[index,'SequenceNumber'] = ','.join(sequenceNumbers)
    exac_raw_GPCRs['SequenceNumber'] = exac_raw_GPCRs['SequenceNumber'].astype(int)

    ## Merge SNP and exacRAW!
    exac_raw_GPCRs['feature_use'] = exac_raw_GPCRs['hgvsc_refine'].str.split('.', expand=True)[0]
    snp = snp.reset_index()
    for index, rows in tqdm(enumerate(snp.iterrows())):
        uniprot = snp[index:index+1]['Uniprot'].values[0]
        try:
            ID = enst[enst.uniprot==uniprot].ENST.values[0]
        except:
            ID = 'nan'
        snp.loc[index,'feature_use'] = ID

    dfm = snp.merge(exac_raw_GPCRs, how='left', on=['feature_use','WTaa','SequenceNumber', 'NMaa'])
    dfm = dfm.dropna(subset=['hgvsc_refine'])
    dfm.to_csv("../data/ExAC_GPCRs_raw.csv")

    ### ========================================================================
    ## +++
    ## API
    ## +++
    ## Export hgvsc:
    exac_raw_GPCRs['SequenceNumber'] = 'nan'
    ## translate hgvsc to protein sequence number
    for index, snp in tqdm(enumerate(exac_raw_GPCRs.iterrows())):
        hgvsc = exac_raw_GPCRs[index:index+1]['hgvsc'].values[0]
        try:
            variant = ensemble_API_variants(hgvsc)
            sequenceNumber = variant[0]['transcript_consequences'][0]['protein_start']
            exac_raw_GPCRs.loc[index,'SequenceNumber'] = sequenceNumber
        except:
            print hgvsc

    ## +++
    ## GET DATA FROM:
    ## VEP (Variant Effect Predictor or try Biomart (http://www.ensembl.org/biomart/)
    ## +++
    ## Export to http://www.ensembl.org/Homo_sapiens/Tools/IDMapper/
    exac_raw_GPCRs['hgvsc'].str.split('.',expand=True)[0].unique()
    mapping = pd.read_csv("../ENSTs_mapping.csv")
    to_remove = mapping[mapping[' New stable ID'].str.contains("<retired>")]['Old stable ID'].str.split('.',expand=True)[0].unique()
    # From Variant Effet Predictor (online)
    exac_raw_GPCRs[~exac_raw_GPCRs['hgvsc'].str.contains('|'.join(to_remove))][['hgvsc']].to_csv("~/Downloads/test.txt", header=False, index=False, quoting=csv.QUOTE_NONE, quotechar='',escapechar=',')

def relative_table():
    ## Make relative Table and Quality Control

    # based on the exac_gpcr file
    GPCRs = pd.read_csv("../data/GPCR_lengths.csv")
    exac_GPCR = pd.read_csv("../data/exac_GPCRs.csv")

    for index, protein in tqdm(enumerate(exac_GPCR.iterrows())):

        try:
            Uniprot = exac_GPCR[index:index+1]['Uniprot'].values[0]
            seqLen = GPCRs[GPCRs['Uniprot']==Uniprot]['seqLen'].values[0]
            # Class = GPCRs[GPCRs['Uniprot']==Uniprot]['Class'].values[0]
            # Family = GPCRs[GPCRs['Uniprot']==Uniprot]['Family'].values[0]

            exac_GPCR.loc[index,'seqLen'] = seqLen
            # exac_GPCR.loc[index,'Class'] = Class
            # exac_GPCR.loc[index,'Family'] = Family
        except:
            print '\nNo information for: ', exac_GPCR[index:index+1]['Uniprot'].values[0]

            exac_GPCR['RelativeSNP'] = exac_GPCR.apply(lambda x: round(float(x['n_mis'])/x['seqLen'],2), axis=1)
            exac_GPCR.reset_index(inplace=True)
            exac_GPCR.to_csv("../data/170131_relative_table_exac.csv")

            ## QUALITY CONTROL AND VALIDATION!
            relative_table = pd.read_csv("../data/170210_relative_table_gag.csv")
            merge = pd.merge(left=relative_table, right=exac_GPCR, how='right', left_on='Uniprot', right_on='Uniprot')
            merge['deltaRelative'] = abs(merge['RelativeSNP_x']-merge['RelativeSNP_y'])
            merge[merge.deltaRelative>0.05][['gene','RelativeSNP_y','RelativeSNP_x','n_mis','snp_count']]

            ## Less than expected
            outliers = merge[(merge.snp_count<merge.n_mis) & (merge.deltaRelative>0.03)][['gene','RelativeSNP_y','RelativeSNP_x','n_mis','snp_count','Uniprot']]
            ## update those!
            for Uniprot in outliers.Uniprot.unique():
                relative_table.loc[relative_table[relative_table.Uniprot==Uniprot].index,'RelativeSNP']=round(merge[merge.Uniprot==Uniprot].RelativeSNP_y.values[0],2)
                relative_table.to_csv("../data/170210_relative_table_gag.csv")

                ## More than expected
                merge[(merge.snp_count>merge.n_mis) & (merge.deltaRelative>0.05)][['gene','RelativeSNP_y','RelativeSNP_x','n_mis','snp_count']]
                ## Manual update of outliers!

def PharmGKB_network():
    ## Prepare for Disease-Drug Network

    dds = pd.read_excel('../data/TableS2.xlsx', sheetname=1)
    dds['Level of Evidence'] = dds['Level of Evidence'].astype(str)
    dds = dds[(dds.Category=='FDA appoved drug') & ((dds['Level of Evidence'].str.contains('3')) | (dds['Level of Evidence'].str.contains('2')))]
    disease_ontologies = pd.read_excel('../data/TableS2.xlsx', sheetname=3)
    dds_targ = pd.concat([pd.Series(row['Target'], row['Diseases'].split(';')) for _, row in dds.iterrows()]).reset_index()
    dds_out = pd.concat([pd.Series(row['Chemical Name'], row['Diseases'].split(';')) for _, row in dds.iterrows()]).reset_index()
    dds_out.rename(columns={'index': 'Disease',0: 'Drug'}, inplace=True)
    dds_out['DrugAttribute'] = 'Drug'
    dds_out['DiseaseAttribute'] = 'Disease'
    dds_out = dds_out.merge(disease_ontologies, on = 'Disease')

    ### ONLINE VIS
    ## http://app.rawgraphs.io/
    dds_variant = pd.concat([pd.Series(row['Chemical Name'], row['Variants'].split(';')) for _,row in dds.iterrows()]).reset_index()
    dds_variant.rename(columns={'index': 'Variant',0: 'Drug'}, inplace=True)
    dds_out = pd.concat([pd.Series(row['Chemical Name'], row['Diseases'].split(';')) for _, row in dds.iterrows()]).reset_index()
    dds_out.rename(columns={'index': 'Disease',0: 'Drug'}, inplace=True)
    dds_out = dds_out.merge(disease_ontologies, on = 'Disease')
    dds_out = dds_out.merge(dds_variant, on = 'Drug')
    dds_out['VariantAA'] = dds_out['Variant'].str.split(' ', expand=True)[1]
    dds_out['VariantRS'] = dds_out['Variant'].str.split(' ', expand=True)[0]

    def af(x):
        raw = dds[dds.Variants.str.contains(x)].Variants.values[0]
        if not ';' in raw:
            return dds[dds.Variants.str.contains(x)].AlleleFrequency.values[0]
        else:
            for i,rs in enumerate(raw.split(';')):
                if rs.split(' (')[0] == x:
                    raw_AF = dds[dds.Variants==raw].AlleleFrequency.values[0]
                    return raw_AF.split(';')[i]

    dds_out['AF'] = dds_out['VariantRS'].apply(lambda x: af(x))

    # LOAD variantRS IDs to proteinIDs from Biomart
    VarTranslate =pd.read_csv('../data/170326_variantTable_translate.csv')
    dds_out = dds_out.merge(VarTranslate, on = 'VariantRS')
    dds_out['ProteinID'] = dds_out['EntryName'].str.upper().str.split('_', expand=True)[0]
    dds_out['Variant-output'] = dds_out['ProteinID']+' '+dds_out['VariantAA']
    dds_out['SequenceNumber'] = dds_out['VariantAA'].apply(lambda x: int(re.findall(r'\d+',x)[0]))
    dds_out.sort_values(by='AF',ascending=False).tail()
    dds_out.to_csv("../data/disease_drug_links.csv",index=False)

    ## =========================================================================

    SegmentList = list(dfm_drugs_master[dfm_drugs_master.EntryName=='adrb2_human'].Segment.unique())
    non_ExAC = {'rs1801253': 'H8', 'rs3787429': 'ECL2', 'rs3787430': 'ECL2', 'rs6166':'N-term', 'rs6280': 'N-term'}
    for index, rows in enumerate(dds_out.iterrows()):
        SN = dds_out[index:index+1]['SequenceNumber'].values[0]
        EntryName = dds_out[index:index+1]['EntryName'].values[0]
        VariantRS = dds_out[index:index+1]['VariantRS'].values[0]
        category = dds_out[index:index+1]['Category'].values[0]
        Segment = dfm_drugs_master[(dfm_drugs_master.SequenceNumber==SN) & (dfm_drugs_master.EntryName==EntryName)].Segment.unique()
        if len(Segment)>0:
            dds_out.loc[index,'Segment'] = Segment[0]
        else:
            Segment = non_ExAC[VariantRS]
            dds_out.loc[index,'Segment'] = Segment
        # ------
        dds_out.loc[index,'Category2'] = int(list(dds_out['Category'].unique()).index(category))
        dds_out.loc[index,'Variant-output2'] = 'abcdefghijklmnopqrstuvq'[(int(SegmentList.index(Segment)))]+'-' + dds_out[index:index+1]['Variant-output'].values[0]

    ## Add numbers for drugnames and variants in front to get the order --> export as svg and remove the text!?
    dds_out['Category2'] = dds_out['Category2'].astype(str)
    dds_out['Drug2'] = dds_out['Category2']+'-'+dds_out['Drug']
    dds_out['Disease2'] = dds_out['Category2']+'-'+dds_out['Disease']

    dds_out.drop_duplicates().to_csv("../data/suppT1_clinicalMVs.csv",index=False)

def exac_drugscore():
    ## Calculate the drugscore for all drugs targeting variant GPCRs

    drugs_raw = pd.read_excel("../data/drug_data.csv")
    dfm_drugs_master = pd.read_csv("../data/dfm_drugs_masterTable.csv")

    drugs = drugs_raw[(drugs_raw.PMID.notnull()) & (drugs_raw['Status']=='approved')]['Drug Name'].unique()
    functional_anno = pd.read_csv('../data/Drugs_vs_SNPs.csv')

    drugscore = {}
    for drug in tqdm(drugs):
        targets = drugs_raw[drugs_raw['Drug Name']==drug].EntryName.unique()

        drugscore_exac = functional_anno[functional_anno.EntryName.isin(targets)].Relative_functional.sum()
        drugscore[drug] = drugscore_exac

    drugsdf = pd.DataFrame(drugscore.items(), columns=['drug','drugscore_known_exac'])
    drugsdf.sort_values(by='drugscore_known_exac', ascending=False)

    writer = ExcelWriter('../SuppTable_drugscore_exac.xlsx')
    drugsdf.to_excel(writer, 'drugscore_known_exac', index=False)
    writer.save()

def G1000K_male_female(rsIds):
    ## Look at male/female rations for specific rsIds

    sample_info = pd.read_csv("../data/Genome1K_samplesinfo.tsv", sep='\t')
    rsIds = 'rs6166' # ass with Infertility;Ovarian hyperstimulation syndrome

    df = genotypes[genotypes.ID==rsIds]
    rs_exac = dfm_drugs_master[dfm_drugs_master.func_known==True].id.unique()
    df = genotypes[genotypes.ID.isin(rs_exac)]

    ## Heterocygous individuals:
    individuals = df.columns[10:][([df[col].str.contains("1").any() for col in df.columns[10:]])]
    females_with_variant = sample_info[(sample_info['Sample name'].isin(individuals)) & (sample_info.Sex=='female')]
    total_females = sample_info[(sample_info['Sample name'].isin(df.columns[10:])) &  (sample_info.Sex=='female')]

    ## percentage of females with that variant allele
    print float(females_with_variant.shape[0])/total_females.shape[0]*100

def G1000K_genotypes():

    ## Combine all chromosome files
    chromosomes = glob.glob("../data/GPCR_eSNPs/gpcr*_eSNPs*")
    df_from_each_file = (pd.read_csv(chromosome+'/'+chromosome.split('/')[-1]+'_GT.txt', sep='\t') for chromosome in chromosomes)
    concatenated_df = pd.concat(df_from_each_file, ignore_index=True)
    concatenated_df.to_csv("~/Downloads/GPCR_eSNPs/170227_1k_GT.csv")

    ## Load data (reanalysis from here)
    genotypes = pd.read_csv("../data/GPCR_eSNPs/170227_1k_GT.csv")
    gene = pd.read_csv('../data/GPCR_eSNPs/GPCR37_info.txt', sep='\t')
    variant_type = pd.read_csv('../data/GPCR_eSNPs/variant_types.csv')
    VEP_variant_type = pd.read_csv('../data/GPCR_eSNPs/VEP_rsIDs.txt', sep='\t')

    ## Drugtargets:
    drugs_raw = pd.read_excel("../data/drug_data.csv")
    drugTargets = drugs_raw[(drugs_raw.PMID.notnull()) & (drugs_raw['Status']=='approved')]['EntryName'].unique()
    translate = pd.read_csv("../data/EST_ENSG_GENE.csv")
    allGPCRs = translate.ENSG.unique()
    drugTargets_ENSGs = translate[translate.EntryName.isin(drugTargets)].ENSG.unique()

    ## Retrieve only MVs
    ## from BioMart: http://www.ensembl.org/biomart/martview/69651c2699acee365afa83b88e1dcb62
    # variant_type = variant_type[variant_type['Variant consequence']=='missense_variant']
    # genotypes = genotypes[genotypes.ID.isin(variant_type['Variant Name'].unique())]
    # df_RsID_drugtargets = variant_type[variant_type['Gene stable ID'].isin(drugTargets_ENSGs)]
    # RsID_drugtargets = df_RsID_drugtargets['Variant Name'].unique()
    # number_of_genes = pd.merge(left=ind, right=variant_type, how='left', left_on='ID', right_on='Variant Name')['Gene stable ID'].nunique() # of genes

    ## from variant effect predictor: http://www.ensembl.org/Homo_sapiens/Tools/VEP/Results?db=core;tl=ybOfJbau5bRjLL9a-2854089
    VEP_variant_type = VEP_variant_type[VEP_variant_type['Consequence']=='missense_variant']
    genotypes = genotypes[genotypes.ID.isin(VEP_variant_type['#Uploaded_variation'].unique())]
    df_RsID_drugtargets = VEP_variant_type[VEP_variant_type['Gene'].isin(drugTargets_ENSGs)]
    RsID_drugtargets = df_RsID_drugtargets['#Uploaded_variation'].unique()

    ## GPCRs with MVS per individual
    list_number_of_mutations = []
    list_number_of_genes = []

    drugtargets_li_numb_mut = []
    drugtargets_li_numb_genes = []
    targets = {}
    li_clinical_variants = []
    variants = {}

    ## FROM TABLE S2
    rsID_clinical_variants = [u'rs1042713', u'rs1799971', u'rs1042714', u'rs1801253',
       u'rs1801252', u'rs6280', u'rs3787429', u'rs3787430', u'rs6314',
       u'rs6318', u'rs4994', u'rs6166', u'rs6165', u'rs10305420',
       u'rs6923761', u'rs1042636']

    for individual in tqdm(genotypes.columns[9:].unique()):
        ind = genotypes[genotypes[individual]!='0|0'][['ID',individual]] # --> Missense mutations in GPCRs

        merge = pd.merge(left=ind, right=VEP_variant_type, how='left', left_on='ID', right_on='#Uploaded_variation')

        ## all GPCRs
        number_of_mutations = ind.ID.nunique() # of MVs
        number_of_genes = merge[merge['Gene'].isin(allGPCRs)].Gene.nunique() # of genes
        list_number_of_mutations.append(number_of_mutations)
        list_number_of_genes.append(number_of_genes)

        ## Clinical variants
        ## FROM TABLE S2: Reported variant-drug response associations
        num_clinical_variants = ind[ind.ID.isin(rsID_clinical_variants)].ID.nunique()
        li_clinical_variants.append(num_clinical_variants)

        for i in ind[ind.ID.isin(rsID_clinical_variants)].ID.unique():
            if not i in variants:
                variants[i] = 0
            variants[i] += 1

        ## drugtargets
        number_of_mutations_drugtargets = ind[ind.ID.isin(RsID_drugtargets)].ID.nunique() # of MVs

        merge_drugs = pd.merge(left=ind[ind.ID.isin(RsID_drugtargets)], right=VEP_variant_type, how='left', left_on='ID', right_on='#Uploaded_variation')
        number_of_drugtargets = merge_drugs[merge_drugs['Gene'].isin(allGPCRs)].Gene.nunique() # of genes

        for i in merge_drugs[merge_drugs['Gene'].isin(allGPCRs)].Gene.unique():
            proteinName = translate[translate.ENSG==i].EntryName.values[0].split('_human')[0].upper()
            if not proteinName in targets:
                targets[proteinName] = 0
            targets[proteinName] += 1

        drugtargets_li_numb_mut.append(number_of_mutations_drugtargets)
        drugtargets_li_numb_genes.append(number_of_drugtargets)

    ## Create a figure instance (Figure 2)
    fig = plt.figure()
    gridspace = [drugtargets_li_numb_mut, drugtargets_li_numb_genes, li_clinical_variants, li_clinical_variants_targets]
    gridlabels = ['Count MVs', 'Count drug targets', 'Count clinical MVs', 'Count clinical MVs in drug target']
    gs = gridspec.GridSpec(2, len(gridspace), width_ratios=[1]*len(gridspace))
    for i, grid in enumerate(gridspace):
        ax = plt.subplot(gs[i])
        # ax = fig.add_subplot(111)
        bp = ax.boxplot([grid], patch_artist=True, showfliers=False) #, li_clinical_variants, li_clinical_variants_targets
        ax.set_ylabel(gridlabels[i], fontsize=16)
        # ax.set_xticklabels(['# of\nmissense variants', '# of polymorphic\ndrugtargets', '# of clinically\nknown MVs'])
        # ax.set_title('missense mutations per individual\nin GPCR drugtargets (n=108)')
        for box in bp['boxes']:
            box.set( color='#000000', linewidth=0.1)
            box.set( facecolor = '#ed1c24' )
        for median in bp['medians']:
            median.set(color='#000000', linewidth=1.5)
        # ax.set_yticks([0,20,40,60,80,100])
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if i <= 1:
            plt.ylim(0,110)
    plt.show()
    plt.savefig("../figures/1k.svg", bbox_inches='tight')

def G1000K_NHS_drugs():
    ## Number of effected individuals per drug

    nhs = pd.read_csv("../data/170829_drugs_sales.csv")
    drugs_raw = pd.read_excel("../data/drug_data.csv")

    drugs = drugs_raw[(drugs_raw.PMID.notnull()) & (drugs_raw['Status']=='approved')]['Drug Name'].unique()

    # nhs_drugs = nhs.merge(drugs, left_on='drugNameQuery', right_on='Drug Name', how='left')
    dfm_drugs_master = pd.read_csv("..../data/dfm_drugs_masterTable.csv")
    genotypes = pd.read_csv("..../data/GPCR_eSNPs/170227_1k_GT.csv")
    gene = pd.read_csv('..../data/GPCR_eSNPs/GPCR37_info.txt', sep='\t')
    variant_type = pd.read_csv('..../data/GPCR_eSNPs/variant_types.csv')
    VEP_variant_type = pd.read_csv('..../data/GPCR_eSNPs/VEP_rsIDs.txt', sep='\t')

    rs_1k_exac = set(dfm_drugs_master[dfm_drugs_master.Functional==True].id.unique()) & set(genotypes.ID.unique())

    ######
    ## FOR EACH DRUG
    ######
    li_target_func_IDs = []
    li_target_func_IDs_homo = []
    variants = {}
    variants_homo = {}
    # for drug in tqdm(nhs_drugs.drugNameQuery.unique()):
    for drug in tqdm(drugs):
        # targets = nhs_drugs[nhs_drugs.drugNameQuery==drug].EntryName.unique()
        targets = drugs_raw[drugs_raw['Drug Name']==drug].EntryName.unique()
        target_func_IDs = dfm_drugs_master[(dfm_drugs_master.EntryName.isin(targets)) & (dfm_drugs_master.id !='.') & ((dfm_drugs_master.func_known==True) | (dfm_drugs_master.func_putative==True))].id.unique()
        #  | (dfm_drugs_master.func_putative==True)
        df = genotypes[genotypes.ID.isin(target_func_IDs)]
        individuals = df.columns[9:][([df[col].str.contains("1").any() for col in df.columns[9:]])]
        homozygous = df.columns[9:][([(df[col]=="1|1").any() for col in df.columns[9:]])]
        li_target_func_IDs.append(len(individuals))
        li_target_func_IDs_homo.append(len(homozygous))

        variants[drug] = len(individuals)
        variants_homo[drug] = len(homozygous)

    drugsdf = pd.DataFrame(variants.items(), columns=['drug','individuals'])
    drugsdf['percentage'] = drugsdf['individuals']/2504*100
    print drugsdf['percentage'].mean()
    drugsdf = drugsdf.sort_values(by='percentage', ascending=False)
    drugsdf.to_csv("../data/drugsdf_putative.csv")

    drugsdf_homo = pd.DataFrame(variants_homo.items(), columns=['drug','individuals'])
    drugsdf_homo['percentage'] = drugsdf_homo['individuals']/2504*100
    print drugsdf_homo['percentage'].mean()
    drugsdf_homo = drugsdf_homo.sort_values(by='percentage', ascending=False)
    drugsdf_homo.to_csv("../data/drugsdf_putative_homo.csv")

    from pandas import ExcelWriter
    writer = ExcelWriter('../SuppTable_drugs_target_func_sites.xlsx')
    drugsdf.to_excel(writer, index=False)
    writer.save()

    plt.figure()
    plt.boxplot([li_target_func_IDs], 1, 'o')
    plt.show()
    ### Create a figure instance
    sns.set(style="white", color_codes=True, font_scale=1)
    fig = plt.figure()
    gridspace = [li_target_func_IDs]
    gridlabels = ['drug on number of individuals\nwith a mutated functional site']
    gs = gridspec.GridSpec(1, len(gridspace), width_ratios=[1]*len(gridspace))
    for i, grid in enumerate(gridspace):
        ax = plt.subplot(gs[i])
        # ax = fig.add_subplot(111)
        bp = ax.boxplot([grid], patch_artist=True) #, li_clinical_variants, li_clinical_variants_targets
        ax.set_ylabel(gridlabels[i], fontsize=12)
        # ax.set_xticklabels(['# of\nmissense variants', '# of polymorphic\ndrugtargets', '# of clinically\nknown MVs'])
        # ax.set_title('missense mutations per individual\nin GPCR drugtargets (n=108)')
        for box in bp['boxes']:
            box.set( color='#000000', linewidth=1)
            box.set( facecolor = '#ef3636' )
        for median in bp['medians']:
            median.set(color='#000000', linewidth=2)
        for flier in bp['fliers']:
            flier.set(marker='o', color='#000000', linewidth=1, alpha=1)
        ax.set_ylim(0,40)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
    plt.show()
    plt.savefig("../figures/1k_drug_individuals.svg", bbox_inches='tight')

def G1000K_functionaL_variants():
    ## functional variant analysis of the 1000 Genomes Project

    # We find that on average, X% of the 2,504 individuals in the 1,000 Genomes Project carry a variant GPCR drug target with at least one missense mutation in a known functional site (Figure 6B; Z% in known or putative functional site; Table SX).
    drugs_raw = pd.read_excel("../data/drug_data.csv")
    drugsTargets = drugs_raw[(drugs_raw['Status']=='approved')].EntryName.unique()

    dfm_drugs_master = pd.read_csv("../data/dfm_drugs_masterTable.csv")

    genotypes = pd.read_csv("../data/GPCR_eSNPs/170227_1k_GT.csv")
    gene = pd.read_csv('../data/GPCR_eSNPs/GPCR37_info.txt', sep='\t')
    variant_type = pd.read_csv('../data/GPCR_eSNPs/variant_types.csv')
    VEP_variant_type = pd.read_csv('../data/GPCR_eSNPs/VEP_rsIDs.txt', sep='\t')

    ######
    ## FOR EACH Individual
    ######
    rs_1k_exac = set(dfm_drugs_master[dfm_drugs_master.func_known==True].id.unique()) & set(genotypes.ID.unique())

    ## Repeat for both putative and known functional sites
    li_receptor_func_IDs = []
    li_receptor_func_IDs_homo = []
    variants = {}
    variants_homo = {}
    for receptor in tqdm(dfm_drugs_master.EntryName.unique()):
        ## known sites:
        # target_func_IDs = dfm_drugs_master[(dfm_drugs_master.EntryName==receptor) & (dfm_drugs_master.id !='.') & (dfm_drugs_master.func_known==True)].id.unique()

        ## known and putative functional sites:
        target_func_IDs = dfm_drugs_master[(dfm_drugs_master.EntryName==receptor) & (dfm_drugs_master.id !='.') & ((dfm_drugs_master.func_known==True) | (dfm_drugs_master.func_putative==True))].id.unique()

        df = genotypes[genotypes.ID.isin(target_func_IDs)]
        individuals = df.columns[9:][([df[col].str.contains("1").any() for col in df.columns[9:]])]
        homozygous = df.columns[9:][([(df[col]=="1|1").any() for col in df.columns[9:]])]
        li_receptor_func_IDs.append(len(individuals))
        li_receptor_func_IDs_homo.append(len(homozygous))

        variants[receptor] = len(individuals)
        variants_homo[receptor] = len(homozygous)

    drugsdf = pd.DataFrame(variants_homo.items(), columns=['EntryName','individuals'])
    drugsdf['percentage'] = drugsdf['individuals']/2504*100
    print drugsdf['percentage'].mean()
    drugsdf = drugsdf.sort_values(by='percentage', ascending=False)
    drugsdf.to_csv("../data/170829_receptor_putative_homo.csv") # _homo

    ## Export
    drugsdf_k = pd.read_csv("../data/170829_receptor_known.csv") # _homo
    drugsdf_k_homo = pd.read_csv("../data/170829_receptor_known_homo.csv") #
    drugsdf_p = pd.read_csv("../data/170829_receptor_putative.csv") # _homo
    drugsdf_p_homo = pd.read_csv("../data/170829_receptor_putative_homo.csv") # _homo

    writer = ExcelWriter('../SuppTable_target_func_sites.xlsx')
    drugsdf_p_homo.to_excel(writer, 'drugsdf_p_homo', index=False)
    drugsdf_k_homo.to_excel(writer, 'drugsdf_k_homo', index=False)
    drugsdf_p.to_excel(writer, 'drugsdf_p', index=False)
    drugsdf_p_homo.to_excel(writer, 'drugsdf_p_homo', index=False)
    writer.save()

def PTMs():
    ## Preparing PTM sites

    ESTs = pd.read_csv("../data/EST_ENSG_GENE.csv")

    ### dbPTM
    dbptm = pd.read_csv("../data/170201_dbPTM3.txt", sep='\t', names=['EntryName','Uniprot','SequenceNumber','Description','Sites','Source','AA','Type'])
    dbptm['EntryName'] = dbptm.EntryName.str.lower()
    dbptm_gpcrs = dbptm[(dbptm.Uniprot.isin(ESTs.Uniprot.unique())) | (dbptm.EntryName.isin(ESTs.EntryName.unique()))]
    dbptm_gpcrs = dbptm_gpcrs.reset_index()
    dbptm_gpcrs = dbptm_gpcrs.merge(ESTs, on='Uniprot', how='left') # ~1400 datapoints
    dbptm_gpcrs.rename(columns={'EntryName_y': 'EntryName'}, inplace=True)

    ### PhosphoSitePlus
    ## README in PTMVar
    ## LTP_LIT: low-throughput experimental techniques
    ## MS2_LIT: proteomic MS experiments / MS_LIT > 0
    ## CST_CS: shotgun proteomic experiments
    datasets = ["Acetylation_site_dataset", "Methylation_site_dataset", "O-GalNAc_site_dataset", "O-GlcNAc_site_dataset", "Phosphorylation_site_dataset", "Sumoylation_site_dataset", "Ubiquitylation_site_dataset"]
    psp_GPCRs = pd.DataFrame(columns=[u'GENE', u'PROTEIN', u'ACC_ID', u'HU_CHR_LOC', u'MOD_RSD', u'SITE_GRP_ID', u'ORGANISM', u'MW_kD', u'DOMAIN', u'SITE_+/-7_AA',
       u'LT_LIT', u'MS_LIT', u'MS_CST', u'CST_CAT'])
    for dataset in datasets:
        ptm = pd.read_csv("../data/PSP_MAR_21_2017/" + dataset , sep='\t', comment='#')
        ptm_gpcrs = ptm[ptm.ACC_ID.isin(ESTs.Uniprot.unique())].reset_index()
        ptm_gpcrs['Type'] = dataset.split('_')[0].replace('O-GlcNAc','O-linked Glycosylation')
        psp_GPCRs = pd.concat([psp_GPCRs, ptm_gpcrs], ignore_index=True)
    psp_GPCRs['Source'] = 'PSP'
    psp_GPCRs['AA'] = psp_GPCRs['MOD_RSD'].str[0]
    psp_GPCRs['SequenceNumber'] = psp_GPCRs['MOD_RSD'].apply(lambda x: int(re.findall(r'\d+', x)[0]))
    psp_GPCRs = psp_GPCRs.merge(ESTs, left_on='ACC_ID', right_on='Uniprot', how='left')
    ptm_gpcrs = psp_GPCRs.merge(dbptm_gpcrs, how='outer') # , on=['EntryName', 'SequenceNumber', 'AA', 'Type'])
    ptm_gpcrs = ptm_gpcrs.drop_duplicates(subset=['EntryName', 'SequenceNumber', 'AA', 'Type'], keep='last').reset_index()

    ### Ben/Screenis Dataset (old and depreciated, no control on quality)
    # ptm = pd.read_csv("unimod_human.txt", sep='\t')
    # [i+'_HUMAN' for i in uni_ensemble.ProteinID.unique()]
    # ptm_gpcrs = ptm[ptm.name.isin([i+'_HUMAN' for i in uni_ensemble.ProteinID.unique()])]
    # ptm_gpcrs.name.value_counts()

    ## Choose either starting set and do the PTM vs. SNP analysis
    ## Map GPCRdb onto SequenceNumber
    ptm_gpcrs['GPCRdb'] = 'none'

    for entryName in tqdm(ptm_gpcrs.EntryName.unique()):
        ## map to GPCRdbEntryNames:
        uniprot = ptm_gpcrs[ptm_gpcrs.EntryName==entryName]['Uniprot'].values[0]

        protein_residue_info = get_protein_residues(entryName)
        for x in ptm_gpcrs[ptm_gpcrs.EntryName==entryName].iterrows():
            index = x[0]
            aa = ptm_gpcrs[index:index+1]['AA'].values[0]
            residue = ptm_gpcrs[index:index+1]['SequenceNumber'].values[0]
            try:
                GPCRdbWT = protein_residue_info[int(residue)-1]['amino_acid']
            except:
                GPCRdbWT = 'none'
            if aa == GPCRdbWT:
                ptm_gpcrs.loc[index,'GPCRdb'] = protein_residue_info[int(residue)-1]['display_generic_number']
                ptm_gpcrs.loc[index,'EntryName'] = entryName
                ptm_gpcrs.loc[index,'use'] = 'yes'
            else:
                ptm_gpcrs.loc[index,'use'] = 'no'
                print entryName, residue, aa, GPCRdbWT
    ptm_gpcrs.drop(['index_x', 'index_y', 'EntryName_x', 'MW_kD', 'ORGANISM'], inplace=True, axis=1)
    ptm_gpcrs = ptm_gpcrs[ptm_gpcrs.use=='yes']
    ptm_gpcrs.to_csv('../data/PTM_gpcrs.csv')

def export(list_of_proteins):
    ## EXPORT a list of proteins from the MasterTable and sort them by prioties for funcitonal testing
    from pandas import ExcelWriter

    dfm_drugs_master = pd.read_csv("../data/dfm_drugs_masterTable.csv")
    subcols = ['EntryName','Uniprot','Class','Family','Segment','GPCRdb','SequenceNumber','id','WTaa','NMaa','MutationType','Allele Frequency','Allele Count','Number of Homozygotes','func_known','func_putative','sift','polyphen','score','diseaseName','Type',  'MicroSwitch','ActivationPathway', 'PTMsite', 'foldchangeMaxAbs', 'ArrestinInteraction','GProteinInteraction', 'LB_structure', 'LB_fam']

    def que(x):
        if x['LB_structure_human'] == 'interacting' or x['LB_structure_ortho'] == 'interacting':
            return 'interacting'
        else:
            'NaN'

    dfm_drugs_master['LB_structure'] = dfm_drugs_master.apply(que, axis=1)

    writer = ExcelWriter('~/Downloads/export.xlsx', engine='xlsxwriter')
    workbook = writer.book

    for protein in list_of_proteins:
        subset = dfm_drugs_master[dfm_drugs_master.EntryName==protein][subcols]
        sheetname = protein.split('_')[0].upper()

        number_rows = len(subset.index)
        # Define our range for the color formatting
        color_range = "F2:F{}".format(number_rows+1)

        subset.sort_values(by=['Allele Count','func_known', 'func_putative'], ascending=False)[subcols].to_excel(writer, sheet_name=sheetname, index=False)
        ## Iterate through each column and set the width == the max length in that column.
        worksheet = writer.sheets[sheetname]

        ## Make changes to formatting:
        ## https://xlsxwriter.readthedocs.io/format.html
        for idx, col in enumerate(subset):  # loop through all columns
            series = subset[col]
            max_len = max((
            series.astype(str).map(len).max(),  # len of largest item
            len(str(series.name))  # len of column name/header
            )) + 1  # adding a little extra space
            worksheet.set_column(idx, idx, max_len)
            for row in range(0,number_rows):
                ## SIFT
                cell = subset[row:row+1]['sift'].values[0]
                format = workbook.add_format()
                format.set_pattern(1)  # for solid fill.
                if 'deleterious' in cell:
                    format.set_bg_color('#e66767')
                else:
                    format.set_bg_color('#68d479')
                worksheet.write(1+row,16,cell,format) # sift

                ## polyphen
                cell = subset[row:row+1]['polyphen'].values[0]
                format = workbook.add_format()
                format.set_pattern(1)  # for solid fill.
                if 'probably_damaging' in cell:
                    format.set_bg_color('#e66767')
                elif 'possibly_damaging' in cell:
                    format.set_bg_color('#ffdd75')
                else:
                    format.set_bg_color('#68d479')
                worksheet.write(1+row,17,cell,format) # polyphen

                ## func_known
                cell = subset[row:row+1]['func_known'].values[0]
                format = workbook.add_format()
                format.set_pattern(1)  # for solid fill.
                if cell == False:
                    format.set_font_color('#e66767')
                else:
                    format.set_font_color('#68d479')
                worksheet.write(1+row,14,str(cell),format)

                ## func_putative
                cell = subset[row:row+1]['func_putative'].values[0]
                format = workbook.add_format()
                format.set_pattern(1)  # for solid fill.
                if cell == False:
                    format.set_font_color('#e66767')
                else:
                    format.set_font_color('#68d479')
                worksheet.write(1+row,15,str(cell),format)
    writer.save()


    def functional_output():
        all_subset = dfm_drugs_master[dfm_drugs_master.EntryName.isin(list_of_proteins)]
        all_subset_func = all_subset[(all_subset.func_known==True)]

        ## ALL
        ALL = all_subset[(all_subset.func_putative==True) & ((all_subset.PTMsite=='yes') | (LB_subset.LB_fam=='interacting') | (SB_subset.ArrestinInteraction=='putative') | (SB_subset.GProteinInteraction=='putative')) ].sort_values(by=['Allele Frequency'], ascending=False)[subcols].reset_index().drop_duplicates(subset=['EntryName','NMaa', 'SequenceNumber','Allele Frequency'], keep='last')
        ALL.to_excel(writer, sheet_name='ALL', index=False)

        ## LB mutations
        LB_subset = dfm_drugs_master[dfm_drugs_master.EntryName.isin([i+'_human' for i in list_of_proteins_lB])]
        LB = LB_subset[(LB_subset.LB_fam=='interacting')].sort_values(by=['Allele Frequency'], ascending=False)[subcols].reset_index().drop_duplicates(subset=['EntryName','NMaa', 'SequenceNumber','Allele Frequency'], keep='last')
        LB.to_excel(writer, sheet_name='LB', index=False)

        ## Signalling mutations
        SB_subset = dfm_drugs_master[dfm_drugs_master.EntryName.isin([i+'_human' for i in list_of_proteins_SB])]
        SB = SB_subset[(SB_subset.ArrestinInteraction=='putative') | (SB_subset.GProteinInteraction=='putative')].sort_values(by=['Allele Frequency'], ascending=False)[subcols].reset_index().drop_duplicates(subset=['EntryName','NMaa', 'SequenceNumber','Allele Frequency'], keep='last')
        SB.to_excel(writer, sheet_name='SB', index=False)

        ## PTM mutations
        PTM_subset = dfm_drugs_master[dfm_drugs_master.EntryName.isin([i+'_human' for i in list_of_proteins_SB+list_of_proteins_lB])]
        PTM = PTM_subset[(PTM_subset.PTMsite=='yes')].sort_values(by=['Allele Frequency'], ascending=False)[subcols].reset_index().drop_duplicates(subset=['EntryName','NMaa', 'SequenceNumber','Allele Frequency'], keep='last')
        PTM.to_excel(writer, sheet_name='PTM', index=False)
