"""

Process SARS-CoV-2 sequences, reconstruct ancestral states and return anc state csv

@author: david
"""
import pandas as pd
from Bio import SeqIO
import AncestralReconstruction
from ete3 import Tree
from pathlib import Path

def fasta2df(fasta_file):
    
    "Get sequences from fasta"
    seq_dic = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_dic[record.id] = [i for i in record.seq]
    df = pd.DataFrame.from_dict(seq_dic, orient='index')
    
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

def plot_tree(tree,mtt_file,fig_file):
    
    import balticmod as bt
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    import TreeUtils
    import seaborn as sns
    
    "Write tree with fit vals to multi-type Newick file"
    absolute_time = 2020.67
    TreeUtils.write_MTT_newick(tree,mtt_file)
    
    sns.set(style="darkgrid")
    myTree=bt.loadNewick(mtt_file,absoluteTime=False)
    myTree.traverse_tree() ## required to set heights
    myTree.setAbsoluteTime(absolute_time) ## set absolute time of all branches by specifying date of most recent tip
    myTree.treeStats() ## report stats about tree
    
    fig,ax = plt.subplots(figsize=(20,20),facecolor='w')

    x_attr=lambda k: k.absoluteTime ## x coordinate of branches will be absoluteTime attribute
    c_func=lambda k: 'darkorange' if k.traits['type']=='1' else 'steelblue' ## colour of branches
    s_func=lambda k: 50-30*k.height/myTree.treeHeight ## size of tips
    z_func=lambda k: 100
    
    cu_func=lambda k: 'k' ## for plotting a black outline of tip circles
    su_func=lambda k: 2*(50-30*k.height/myTree.treeHeight) ## black outline in twice as big as tip circle 
    zu_func=lambda k: 99
    myTree.plotTree(ax,x_attr=x_attr,colour_function=c_func) ## plot branches
    myTree.plotPoints(ax,x_attr=x_attr,size_function=s_func,colour_function=c_func,zorder_function=z_func) ## plot circles at tips
    myTree.plotPoints(ax,x_attr=x_attr,size_function=su_func,colour_function=cu_func,zorder_function=zu_func) ## plot circles under tips (to give an outline)
    
    "Add legend the hard way"
    import matplotlib.patches as mpatches
    blue_patch = mpatches.Patch(color='steelblue', label = 'WT')
    red_patch = mpatches.Patch(color='darkorange', label = 'DEL')
    handles = [blue_patch,red_patch]
    
    "For Regions"
    ax.legend(handles=handles,prop={'size': 24}) #loc='upper left'
    
    "Add month labels as xticks"
    # step_freq = 1/12
    # xticks = np.arange(2020,absolute_time,step_freq)
    # ax.set_xticks(xticks)
    # labels = ['Jan','Feb','Mar','Apr','May','June','July','Aug','Sep']
    # ax.set_xticklabels(labels, fontsize=24) #rotation='vertical'
    # ax.set_xlabel('Time', fontsize=24)
    
    ax.set_ylim(-5,myTree.ySpan+5)
    plt.savefig(fig_file, dpi=300)

plot = True # plot tree colored by anc states?
base_dir = Path(__file__).parent.parent / "data"

"Input files"
meta_file = base_dir / "hcov_USA_post2020-09-01_EMSI_subsampled_metadata.csv"
tree_file = str(base_dir / "hcov_USA_post2020-09-01_EMSI_mini.tre")
align_file = base_dir / 'hcov_USA_post2020-09-01_EMSI_subsampled_aligned.fasta'

"Output files"
labeled_tree_file = str(base_dir / "hcov_USA_post2020-09-01_labeled_forCatherine.tre") # tree with internal labels
anc_state_csv = base_dir / "hcov_USA_post2020-09-01_anc_states_forCatherine.csv"

"Get meta data file"
merged_df = pd.read_csv(meta_file,sep=",",index_col='accession_id')

"Add deletion features"
"Encode nsp6 deletion as binary variable in meta data"
del_label = 'nsp6_Delta9'
merged_df = get_nsp6_del_state(merged_df,align_file,del_label)

"Encode ORF9 TRS deletion as binary variable in meta data"
del_label = 'ORF9_TRS-3'
merged_df = get_ORF9_del_state(merged_df,align_file,del_label)

"Get features to reconstruct"
features = ['nsp6_Delta9', 'ORF9_TRS-3']
#features = merged_df.columns.to_list()

"First assign internal labels to tree and then create new df for ancestral features"
tree = Tree(tree_file, format=1)
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
        if plot:
            mtt_file = "hcov_USA_post2020-09-01_MPreconstruct" + f + ".tre"
            mtt_file = str(base_dir / mtt_file)
            fig_file = "hcov_USA_post2020-09-01_MPreconstruct" + f + ".png"
            fig_file = str(base_dir / fig_file)
            plot_tree(tree,mtt_file,fig_file)
    tree.write(format=1, outfile=labeled_tree_file)

"Save unencoded ancestral states df"
anc_df.to_csv(anc_state_csv,index_label='node')
#anc_df = pd.read_csv(anc_state_csv,sep=",",index_col='node')


    
    