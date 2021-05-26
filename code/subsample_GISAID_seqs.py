"""
Created on Fri Feb 12 06:37:59 2021

Runs bioinformatic pipeline for global trees/seqs from GISAID

@author: david
"""
from pathlib import Path
import pandas as pd
import numpy as np
import dendropy
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
import time
import subprocess
import sys
import re


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
    for index, row in df.iterrows():
        
        gloc = row['virus_name'].split('/')[2] # after 2nd backslash
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

def fasta2df(fasta_file):
    
    "Get sequences from fasta"
    seq_dic = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_dic[record.id] = [i for i in record.seq]
    df = pd.DataFrame.from_dict(seq_dic, orient='index')
    
    return df

def date_tree(tree_file,tip_date_file,rate_file,verbose = False):
    
    cluster = False
    if cluster:
        cmd = '~/lsd -i ' + tree_file + ' -d ' + tip_date_file + ' -c -w ' + rate_file
    else:
        cmd = 'lsd -i ' + tree_file + ' -d ' + tip_date_file + ' -c -w ' + rate_file
    try:
        output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
        if verbose:
            sys.stdout.write(output)
    except subprocess.CalledProcessError as exc:
        print(exc.output)
        print('Execution of "%s" failed!\n' % cmd)
        sys.exit(1)
        
    "Parse newick tree and remove date annotation"
    f = open(tree_file + '.result.date.nexus')
    line = f.readline()
    while line:
        if "tree 1 =" in line:
            tree = line.split()[3]
        line = f.readline()
    f.close()
    tree = re.sub("[\[].*?[\]]", "", tree) # remove bracketed dates
    out_tree = tree_file.split('.')[0] + '_dated.tre'
    nwk=open(out_tree,"w")
    nwk.write(tree + '\n')
    nwk.close()    

def write_tip_dates(tip_date_dict,date_file):
    
    txt=open(date_file,"w")
    txt.write(str(len(tip_date_dict)) + '\n')
    for k,v in tip_date_dict.items():
        txt.write(k + '\t' + str(v) + '\n')
    txt.close()
    
def rename_seqs(df,records,del_label):
    
    "Iterate through list of records, renaming each as we go"
    new_records = []
    for rec in records:
        date = df.loc[rec.name]['collection_date']
        del_state = df.loc[rec.name][del_label]
        if del_state:
            del_state = 'Present'
        else:
            del_state = 'Absent'
        rec.id = rec.name + '_' + del_state + '_' + date
        rec.description = '' # set description blank
        new_records.append(rec)

    return new_records

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

def subsample_align(sampled_taxa,aln_file,new_aln_file):
    
    "Get sequences in tree from GISAID fasta file"
    seq_dict = SeqIO.to_dict(SeqIO.parse(aln_file, "fasta")) # one line alternative
    seq_records = []
    missing_taxa = []
    for tx in sampled_taxa:
        rec = seq_dict.get(tx)
        if rec is not None:
            seq_records.append(rec)
        else:
            missing_taxa.append(tx)
            print('WARNING: ' + tx + " is not in GISAID sequence file")
    SeqIO.write(seq_records, new_aln_file, "fasta")

"""
    Set up directories and files
"""

"GISAID input files"
# base_dir = Path(__file__).home() / 'Documents' / 'GitHub' / 'phyloTF2' / "covid-analysis"
# tree_dir = base_dir / "phylogenies" / 'GISAID-hCoV-19-phylogeny-2021-03-08'
# align_dir = Path.home() / 'Desktop' / 'msa_0314'
# tree_file = str(tree_dir / "hcov_march2021_USA_post2020-09-01.tree")

"Subsample full sequence data to get smaller EMSI data set"
# meta_file = str(tree_dir / 'hcov_USA_post2020-09-01_EMSI_metadata.csv')
# aln_fasta_file = str(align_dir / "hcov_march2021_USA_post2020-09-01_aligned.fasta")
# sub_meta_file = str(tree_dir / "hcov_USA_post2020-09-01_EMSI_subsampled_metadata.csv")
# sub_aln_file = str(align_dir / "hcov_USA_post2020-09-01_EMSI_subsampled_aligned.fasta")
# meta_df = pd.read_csv(meta_file,sep=",",index_col='accession_id')
# sample_count = 2000
# sub_df = meta_df.sample(n=sample_count, axis=0)
# sub_df.to_csv(sub_meta_file,index=True)
# sampled_taxa = sub_df.index.tolist()
# subsample_align(sampled_taxa,aln_fasta_file,sub_aln_file)

"Subsampled data set file names"
base_dir = Path(__file__).parent.parent / "data"
meta_file = str(base_dir / "hcov_USA_post2020-09-01_EMSI_subsampled_metadata.csv")
aln_fasta_file = str(base_dir / "hcov_USA_post2020-09-01_EMSI_subsampled_aligned.fasta")

"""
    Create new sub_sampled data set
"""

"Get metadata for GISAID global tree"
meta_df = pd.read_csv(meta_file,sep=",",index_col='accession_id')

"Encode nsp6 deletion as binary variable in meta data"
#del_label = 'nsp6_Delta9'
#meta_df = get_nsp6_del_state(meta_df,aln_fasta_file,del_label)

"Encode ORF9 TRS deletion as binary variable in meta data"
del_label = 'ORF9_TRS-3'
meta_df = get_ORF9_del_state(meta_df,aln_fasta_file,del_label)

"Split df in presence/absence of deletion"
sub_df = meta_df[meta_df[del_label] == 0]
del_sub_df = meta_df[meta_df[del_label] == 1]

"Subsample seqs without deletion"
sample_count = len(del_sub_df.index)
sub_df = sub_df.sample(n=sample_count, axis=0)

"Merge back into one dataframe"
sub_df = sub_df.append(del_sub_df)

"Get subsampled alignment for samples in sub_df"
sampled_taxa = sub_df.index.tolist()
#new_aln_file = str(base_dir / "hcov_USA_post2020-09-01_EMSI_nsp6_Delta9_subsampled_aligned.fasta") 
new_aln_file = str(base_dir / "hcov_USA_post2020-09-01_EMSI_ORF9_TRS-3_subsampled_aligned.fasta") 
subsample_align(sampled_taxa,aln_fasta_file,new_aln_file)

"Rename records in fasta file by accession id to match taxa labels in tree"
records = SeqIO.parse(new_aln_file, "fasta")
records = rename_seqs(sub_df,records,del_label) # reanme by accession id
SeqIO.write(records,new_aln_file, "fasta")

"Extract tree of desired samples using dendropy"
# filtered_taxa = set(filtered_taxa) # index is 'accession_id'
# filtered_taxa = filtered_taxa.difference(set(missing_taxa)) # remove missing taxa
# filtered_taxa_spaces = [t.replace('_',' ') for t in filtered_taxa] # need to replace underscores with spaces to match dendropy taxon labels
# taxa = dendropy.TaxonNamespace()
# global_tree = dendropy.Tree.get(file=open(tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa) 
# if not global_tree.is_rooted:
#     print('WARNING: Global tree is not rooted') # Should be rooted, this is just to check
# taxa_to_retain = set([taxon for taxon in global_tree.taxon_namespace if taxon.label in filtered_taxa_spaces])
# filtered_tree = global_tree.extract_tree_with_taxa(taxa=taxa_to_retain)
# filtered_tree.write(path=filtered_tree_file,schema='newick',suppress_annotations=True,suppress_rooting=True)

"Get tip dates file for dating with lsd and date tree"
# tip_date_dict = {}
# for sample, row in meta_df.iterrows():
#  	tip_date_dict[sample] = date2FloatYear(row['collection_date'])
# write_tip_dates(tip_date_dict,tip_date_file)
# date_tree(filtered_tree_file,tip_date_file,rate_file)



 