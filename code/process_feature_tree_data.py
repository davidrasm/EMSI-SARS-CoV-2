"""
Created on Wed Aug 19 18:55:30 2020

Process SARS-CoV-2 sequences, reconstruct ancestral states and return fitness features for model trainng

@author: david
"""

import pandas as pd
import re
import dendropy
from Bio import SeqIO
import AncestralReconstruction
import numpy as np
from ete3 import Tree
from pathlib import Path

def fasta2df(fasta_file):
    
    "Get sequences from fasta"
    seq_dic = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_dic[record.id] = [i for i in record.seq]
    df = pd.DataFrame.from_dict(seq_dic, orient='index')
    
    return df

def filter_sites(df,freq_cutoff = 0.0):
    
    "Code mostly taken from filter_invariants.py function:"
    "https://github.com/btmartin721/raxml_ascbias/blob/master/ascbias.py"
    
    #bases = ["A","G","C","T"]
    alphabet = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

    invariant_lst = list()

    # Loop through each dataframe column
    for i in df.columns:

        # Gets unique values at each column and saves to list
        column_unique = df[i].unique().tolist()

        # Intersects column_unique with bases list
        intersect = [value for value in alphabet if value in column_unique]

        # If column contains only ambigous or IUPAC characters
        # Save the column index for dropping later
        if not any(value for value in alphabet if value in column_unique):
            invariant_lst.append(i)
        
        # Get frequency of each variant
        freqs = df[i].value_counts(normalize=True)
        if len(freqs.index) > 1:
            major_var = freqs.index[1]
            major_var_freq = freqs[1]
            if major_var == 'X' and len(freqs.index) > 2: # was major_var == '-'
                major_var = freqs.index[2]
                major_var_freq = freqs[2]
            if major_var == '-': # added to catch case were major variant is '-" in translated seqs
                major_var_freq = 0
        else:
            major_var_freq = 0
        
        # If site is invariant (only A, C, G, or T); ignores N's and "-"
        # OR frequency of major variant is less than frequency threshold cutoff
        if len(intersect) == 1 or major_var_freq < freq_cutoff:

            # Saves column indexes to list
            invariant_lst.append(i)

    # Drops invariant sites from dataframe
    df.drop(invariant_lst, axis=1, inplace=True)

    return df

def filter_variants(df,freq_cutoff):
    
    drop_lst = list()

    # Loop through each dataframe column
    for i in df.columns:

        # Get frequency of each variant
        freqs = df[i].value_counts(normalize=True)
        if len(freqs.index) > 1:
            var_freq = freqs[1]
        else:
            var_freq = 0 # Presumably this dosn't happen but just in case
        
        # If frequency of variant is less than frequency threshold cutoff
        if var_freq < freq_cutoff:
            drop_lst.append(i)

    # Drop sites from dataframe
    df.drop(drop_lst, axis=1, inplace=True)

    return df

def replace_missing(df):
    
    "Replace missing characters with consensus type"
    for i in df.columns:
        consensus = df[i].value_counts().index[0]
        if consensus == '-':
            print("WARNING: missing value '-' is consensus type at: ",str(i))
        df[i].replace(['-'], consensus, inplace=True)
    return df

def drop_consensus_dummies(df):
    
    "Drop dummy variable columns for consensus type"
    drop_lst = list()
    for i in df.columns:
        if '_' in i:
            var_str = i.split('_')[1]
        else:
            var_str = i            
        splits = re.split(r'(\d+)', var_str)
        consensus = splits[0]
        var = splits[-1]
        if consensus == var: # If variant is consensus
            drop_lst.append(i)
    df.drop(drop_lst, axis=1, inplace=True)
    
    return df

def drop_root_state_dummies(df):
    
    "Drop dummy variable columns for consensus type"
    drop_lst = list()
    for i in df.columns:
        if df.loc['root'][i]: # If dummy represents root state
            drop_lst.append(i)
    df.drop(drop_lst, axis=1, inplace=True)
    
    return df

def drop_X_variants(df):
    
    "Drop dummy variable columns for consensus type"
    drop_lst = list()
    for i in df.columns:
        if i[-1] == 'X': # if undetermined
            drop_lst.append(i)
    df.drop(drop_lst, axis=1, inplace=True)
    
    return df

def rename_sites_defunct(df,prefix=''):
    
    for i in df.columns:
        consensus = df[i].value_counts().index[0]
        new_label = prefix + consensus+str(i+1) # index from 1
        df.rename(columns={i: new_label},inplace=True)
        
    return df

def rename_sites(df,prefix='',method='consensus'):
    
    for i in df.columns:
        if method == 'consensus':
            consensus = df[i].value_counts().index[0]
            new_label = prefix + consensus+str(i+1) # index from 1
        elif method == 'root':
            root_state = df.loc['root'][i]
            splits = i.split('_')
            if len(splits) > 1:
                prefix = splits[0] + '_'
                variant = splits[1]
            else:
                prefix = ''
                variant = i
            new_label = prefix + root_state + variant[1:]
        else:
            print("Renaming method not recognized")
            new_label = i        
        df.rename(columns={i: new_label},inplace=True)
        
    return df

def get_geographic_locs(df,column_name,mapToRegion=False):
    
    path = './covid-analysis/maps/'
    code_map_file = path + 'code2StateMap.csv'
    state_map_file = path + 'state2HHSRegionMap.csv'
    
    code_df = pd.read_csv(code_map_file, index_col=0)
    codeMap = {index:row['State'] for index, row in code_df.iterrows()}
    
    state_df = pd.read_csv(state_map_file, index_col=0)
    stateMap = {index:row['Region'] for index, row in state_df.iterrows()}
    
    drop_indexes = []
    locations = []
    #for idx in df.index:
    for index, row in df.iterrows():
        
        gloc = row['virus_name'].split('/')[2] # if taking from GISAID metadata file
        #gloc = idx.split('/')[2] # after 2nd backslash
        gloc = gloc[:2] # first two chars should be state
        
        if gloc in codeMap:
            state_loc = codeMap[gloc]
            if mapToRegion:
                mapped_loc = stateMap[state_loc]
            else:
                mapped_loc = state_loc
        else:
            mapped_loc = 0 # 0 represents unknown
            print("Unrecognized geographic location: " + gloc)
        
        locations.append(mapped_loc)
    
    df[column_name] = locations
    df.drop(drop_indexes,inplace=True) # drop unknowns
        
    return df


def assign_pango_lineage_from_gisaid(df,column_name,pango_df):
    
    "Assign lineages"
    lines = []
    for index, row in df.iterrows():
        if index in pango_df.index:
            lines.append(pango_df.loc[index]['covv_lineage'])
        else:
            lines.append('None')
    df[column_name] = lines
    df[column_name].replace([np.nan], 'None', inplace=True)
    
    "Reduce # of lineages by iterativly mapping rare lineages to parent lineage until parent meets inclusion criteria"
    line_counts = df[column_name].value_counts()
    lineages_to_retain = line_counts[line_counts>=35].index # using 300 for paper results, was 100
    lineages_to_retain = lineages_to_retain.insert(0,'B.1.526')
    lineages_to_retain = lineages_to_retain.insert(0,'B.1.427')
    lineMap = {}
    for index, value in line_counts.items():
        if index in lineages_to_retain:
            lineMap[index] = index
        else:
            parent = index
            while parent not in lineages_to_retain:
                splits = parent.split('.')
                if len(splits) == 1:
                    parent = 'LowFreq' # splits[0]
                    break
                parent = '.'.join(splits[0:-1]) # returns parent e.g. B.1.1.7 --> B.1.1
            lineMap[index] = parent
    
    "Re-assign lineages"
    lines = []
    for index, row in df.iterrows():
        lines.append(lineMap[row[column_name]])
    df[column_name] = lines
        
    return df



def get_nsp6_del_state(df,align_file,column_name):
    
    "This is specifically hardcoded the nsp6 del9 deletion"
    seq_dict = SeqIO.to_dict(SeqIO.parse(align_file, "fasta"))
    deletions = []
    for index, row in df.iterrows():
        rec = seq_dict.get(index)
        seq_str = str(rec.seq[11287:11296]) # genomic start pos is shifted -1 for zero-based indexing
        del_str = '-'*9
        if seq_str == del_str:
            deletions.append(1)
        else:
            deletions.append(0)
            
    df[column_name] = deletions
    
    return df

def get_ORF9_del_state(df,align_file,column_name):
    
    "This is specifically hardcoded the ORF9 TRS -3 deletion"
    seq_dict = SeqIO.to_dict(SeqIO.parse(align_file, "fasta"))
    deletions = []
    for index, row in df.iterrows():
        rec = seq_dict.get(index)
        seq_str = str(rec.seq[28270]) # genomic start pos is shifted -1 for zero-based indexing
        del_str = '-'
        if seq_str == del_str:
            deletions.append(1)
        else:
            deletions.append(0)
            
    df[column_name] = deletions
    
    return df


def get_features(align_file,as_string=True):
    
    df = pd.read_csv(align_file, index_col=0)
    features = [column for column in df.columns]
    if as_string:
        features = ' '.join(features)
    return features


analysis_dir = Path("/Users/david/Documents/GitHub/phyloTF2/covid-analysis")
tree_dir = analysis_dir / "phylogenies" / 'GISAID-hCoV-19-phylogeny-2021-03-08'
base_dir = Path(__file__).parent.parent / "data"

meta_file = base_dir / "hcov_USA_post2020-09-01_EMSI_subsampled_metadata.csv"
tree_file = str(tree_dir / "hcov_march2021_USA_post2020-09-01_dated.tre")
subtree_file = str(base_dir / "hcov_USA_post2020-09-01_EMSI_mini.tre")
labeled_tree_file = str(tree_dir / "hcov_USA_post2020-09-01_EMSI_mini_labeled.tre") # tree with internal labels

align_file = base_dir / 'hcov_USA_post2020-09-01_EMSI_subsampled_aligned.fasta'
unencoded_csv = base_dir / "hcov_USA_post2020-09-01_unencodedAncStates_dels+pangoline.csv"

"Get meta data file"
merged_df = pd.read_csv(meta_file,sep=",",index_col='accession_id')
filtered_taxa = merged_df.index.tolist()

"Extract tree of desired samples using dendropy"
#subtree_file = str(tree_dir / "hcov_USA_post2020-09-01_EMSI_mini.tre")
#filtered_taxa_spaces = [t.replace('_',' ') for t in filtered_taxa] # need to replace underscores with spaces to match dendropy taxon labels
#taxa = dendropy.TaxonNamespace()
#tree = dendropy.Tree.get(file=open(tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa) 
#taxa_to_retain = set([taxon for taxon in tree.taxon_namespace if taxon.label in filtered_taxa_spaces])
#filtered_tree = tree.extract_tree_with_taxa(taxa=taxa_to_retain)
#filtered_tree.write(path=subtree_file,schema='newick',suppress_annotations=True,suppress_rooting=True)

"Get geographic locations"
#merged_df = get_geographic_locs(merged_df,'STATE',mapToRegion=False) # get states
#merged_df = get_geographic_locs(merged_df,'REGION',mapToRegion=True) # get states

"Work around since pangolin lineage assignments are not working"
desktop_dir = Path('/Users/david/Desktop/')
gisaid_metadata_file = desktop_dir / "GISAID_metadata_2021-03-10_USA.csv"
pangolin_df = pd.read_csv(gisaid_metadata_file,sep=",",index_col='covv_accession_id')
merged_df = assign_pango_lineage_from_gisaid(merged_df,'PANGOLINE',pangolin_df)
pangolin_lineage_counts = merged_df['PANGOLINE'].value_counts()
print(pangolin_lineage_counts)

"Add deletion features"
"Encode nsp6 deletion as binary variable in meta data"
del_label = 'nsp6_Delta9'
merged_df = get_nsp6_del_state(merged_df,align_file,del_label)

"Encode ORF9 TRS deletion as binary variable in meta data"
del_label = 'ORF9_TRS-3'
merged_df = get_ORF9_del_state(merged_df,align_file,del_label)

# "Get features to reconstruct"
features = ['PANGOLINE', 'nsp6_Delta9', 'ORF9_TRS-3']
# #features = merged_df.columns.to_list()

# "First assign internal labels to tree and then create new df for ancestral features"
tree = Tree(subtree_file, format=1)
tree = AncestralReconstruction.label_internal_nodes(tree)
node_labels = [node.name for node in tree.traverse("preorder")]
anc_df = pd.DataFrame(index=node_labels)

"Run MP ancestral state reconstruction for all features"
reconstruct = True # set false if already reconstructed
if reconstruct:
    for f in features:
        print('Reconstructing: ' + f)
        feature_dic = {}
        for index, row in merged_df.iterrows():
            feature_dic[index] = row[f]
        tree, anc_dic = AncestralReconstruction.reconstruct_MP(tree,feature_dic)
        anc_df[f] = pd.Series(anc_dic)
    tree.write(format=1, outfile=labeled_tree_file)

# "Save unencoded ancestral states df"
anc_df.to_csv(unencoded_csv,index_label='node')
# #anc_df = pd.read_csv(unencoded_csv,sep=",",index_col='node')


    
    