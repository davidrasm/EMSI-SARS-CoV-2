"""
Created on Fri Feb 12 06:37:59 2021

Runs bioinformatic pipeline for global trees/seqs from GISAID

@author: david
"""
from pathlib import Path
import pandas as pd
import numpy as np
from Bio import AlignIO
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns

def make_features_map():
    
    features = ['nsp01','nsp02','nsp03','nsp04'
		,'nsp05','nsp06','nsp07','nsp08'
		,'nsp09','nsp10','nsp11','nsp12'
		,'nsp13','nsp14','nsp15','nsp16'
		,'Spike','ORF03a','ORF04','ORF05','ORF06'
		,'ORF07a','ORF07b','ORF08','ORF09','ORF10']
    start_locs = [266,806,2720,8555,10055,10973,11843,12092,12686,13025,13442,13442,16237,18040,19621
    			,20659,21563,25393,26245,26523,27202,27394,27756,27894,28274,29558]
    end_locs = [805,2719,8554,10054,10972,11842,12092,12685,13024,13441,13483,16236,18039,19620,20658
    			,21555,25384,26220,26472,27192,27387,27759,27887,28259,29533,29674]
    feature_data = {'Feature':features,'StartPosition': start_locs, 'EndPosition': end_locs}   
    feature_df = pd.DataFrame(feature_data)
    feature_df.to_csv('hcov_genomic_features.csv',index=False)

def add_feature(ax, label, start, end, offset_cntr):
    
    "Add gene or feature annotation to plot ax from start to end"
    left = start
    bottom = (-0.2 - (offset_cntr*(0.1))) * ax.get_ylim()[1]
    width = end - start
    height = 0.08 * ax.get_ylim()[1]
    p = patches.Rectangle(
        (left, bottom), width, height,
        fill=True, color='gray', transform=ax.transData, clip_on=False, alpha=0.75
        )
    ax.add_patch(p)
    ax.text(left+(width/2), bottom+(height/2), label,
        horizontalalignment='center',
        verticalalignment='center',
        transform=ax.transData)

"""
    Set up directories and files
"""

"Set directories"
base_dir = Path(__file__).parent.parent / "data"

"Subsampled data set file names"
aln_fasta_file = str(base_dir / "hcov_USA_post2020-09-01_EMSI_subsampled_aligned.fasta")

feature_file = base_dir / 'hcov_genomic_features.csv'
feature_df = pd.read_csv(feature_file,sep=",",index_col='Feature')

"Count deletions at each site"
align = AlignIO.read(aln_fasta_file, "fasta")
n_sites = align.get_alignment_length()
n_samples = len(align[:,0])
deletions = []
for i in range(n_sites):
    print(i)
    deletions.append(align[:,i].count('-'))
site_del_freq = np.array(deletions) / n_samples 
site_del_freq[:200] = np.nan # ignore first 200 sites
site_del_freq[-200:] = np.nan # ignore last 200 sites

"Plot deletion freqs across genome"
sns.set(style="darkgrid")
fig, ax = plt.subplots(1, 1)
ax.plot(site_del_freq)
#ax.set_xlabel('Position')
ax.set_ylabel('Deletion Frequency')
ax.grid(True)
cntr = 0
for feature, row in feature_df.iterrows():
    add_feature(ax,feature, row['StartPosition'], row['EndPosition'],cntr) # add genomic feature annotation
    cntr += 1
    if cntr > 3: cntr = 0
fig.set_size_inches(10, 6)
fig.tight_layout()
plt.show()

fig.savefig('hcov_USA_post2020-09-01_site_del_freqs.png', dpi=200)



 