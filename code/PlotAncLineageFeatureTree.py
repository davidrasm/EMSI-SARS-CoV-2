"""
Created on Thu Mar 28 12:14:01 2019

Plot tree using baltic with lineages colored according to their fitness

@author: david
"""
import balticmod as bt
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
import TreeUtils
from ete3 import Tree
import seaborn as sns
import random
import json
from pathlib import Path

def add_line_feature(tree,df,feature):
    
    for node in tree.traverse("postorder"):
        print(node.name)
        if node.name == '':
            node.name = 'root'
        node.add_features(state = df.loc[node.name][feature])
    
    return tree

def sample_tips(tree,number=None,fraction=None):

    "Thin tree for plotting by randomly pruning samples"
    #sample_tree = copy.deepcopy(tree)  # Create sample tree
    names = [node.name for node in tree.traverse() if node.is_leaf()]
    leafs = len(names)
    
    sample_count = 0
    if not number and not fraction:
        print("Need to specify a number or fraction to subsample")
    elif number:
        if number > leafs:
            print("Sample size greater than number of leafs in tree")
        else:
            sample_count = number    
    elif fraction:
            sample_count = round(fraction * leafs)
                
    "Return subsample"   
    return random.sample(names, sample_count)   


def thin_tree(tree,sample):
    
    problem_sample = 'hCoV-19/USA/WA-UW-1731/2020|EPI_ISL_424231|2020-03-21'
    if problem_sample in sample:
        print("Removing sample")
        sample.remove(problem_sample)

    print("Pruning tree")
    tree.prune(sample, preserve_branch_length=True)

    return tree

def transform_fvals(fit_vals, tr):
    
    "Rescale fitness values"
    max_fit = np.max(fit_vals)
    min_fit = np.min(fit_vals)
    new_max_fit = max_fit - min_fit # rescale
    fit_vals = np.array([])
    for k in tr.Objects: ## iterate over a flat list of branches
        k.traits['type'] = (float(k.traits['type']) - min_fit) / new_max_fit
        fit_vals = np.append(fit_vals, float(k.traits['type']))
        
    return fit_vals, tr, min_fit, max_fit

"Load in ancestral features as df before encoding"
base_dir = Path(__file__).parent.parent / "data"

features_file = str(base_dir / 'hcov_USA_post2020-09-01_unencodedAncStates_dels+pangoline.csv')
tree_file = str(base_dir / "hcov_USA_post2020-09-01_EMSI_mini_labeled.tre") # tree with internal labels

anc_df = pd.read_csv(features_file,sep=",",index_col='node')

features = ['nsp6_Delta9', 'ORF9_TRS-3']

"Import tree and index branches"
tree = Tree(tree_file, format=1)
absolute_time = 2021.1616 # absolute time of last sample ('2021-02-01')

"Reduce PANGOLINE lineages"
pango_map = {'LowFreq':'LowFreq',
             'B.1':'B.1',
             'B.1.2':'B.1.2',
             'B.1.234':'B.1',
             'B.1.243':'B.1',
             'B.1.427':'B.1.42(7/9)',
             'B.1.429':'B.1.42(7/9)',
             'B.1.526':'B.1.526',
             'B.1.1.7':'B.1.1.7',
             'B.1.1.222':'B.1.1.222'}
anc_df['PANGOLINE'] = anc_df['PANGOLINE'].apply(lambda x: pango_map[x])

"Map anc states to integers"
feature = 'PANGOLINE'
states = anc_df[feature].unique()
num_states = len(states)
state2int = {state:index for index, state in enumerate(states)}
if feature == 'REGION':
    labels = ['REGION_' + str(i) for i in states]
else:
    labels = states
int2Label =  {index:label for index, label in enumerate(labels)}
anc_df[feature] = anc_df[feature].apply(lambda x: state2int[x])

"If thinning by subsampling leaf taxa"
resample = True
if resample:
    sample = sample_tips(tree,number=800)
    with open("sampled_tips_del_tree.json", "w") as fp:
        json.dump(sample, fp)
else:
    with open("sampled_tips_del_tree.json", "r") as fp:
        sample = json.load(fp)

"Add ancestral features"
tree = add_line_feature(tree,anc_df,feature)
tree = thin_tree(tree,sample)
tree, tree_times = TreeUtils.add_tree_times(tree)
final_time = max(tree_times)
root_time = absolute_time - final_time

"Write tree with fit vals to multi-type Newick file"
mtt_file = 'hcov_PANGO_post2020-09-01_del_features.tre'
fig_file = 'hcov_PANGO_post2020-09-01_del_features.png'
TreeUtils.write_MTT_newick(tree,mtt_file)

"""
    Plot tree with lineages colored by fitness
"""
sns.set(style="dark")
myTree=bt.loadNewick(mtt_file,absoluteTime=False)
myTree.traverse_tree() ## required to set heights
myTree.setAbsoluteTime(absolute_time) ## set absolute time of all branches by specifying date of most recent tip
myTree.treeStats() ## report stats about tree

cmap = sns.color_palette("husl", 10) #mpl.cm.get_cmap('tab10', num_states)
location_cmap = sns.color_palette("muted", 10) #mpl.cm.get_cmap('tab10', 10)

fig,ax = plt.subplots(figsize=(20,20),facecolor='w')

gs = GridSpec(1,2,width_ratios=[6,1],wspace=0.0)
ax_tree = plt.subplot(gs[0])
ax_genome = plt.subplot(gs[1],sharey=ax_tree)

x_attr=lambda k: k.absoluteTime ## x coordinate of branches will be absoluteTime attribute
#c_func=lambda k: cmap(int(k.traits['type'])-1) if int(k.traits['type'])>1 else 'black' # for regions, set region 0 to black
c_func=lambda k: cmap[int(k.traits['type'])] # for regions, set region 0 to black
s_func=lambda k: 50-30*k.height/myTree.treeHeight ## size of tips
z_func=lambda k: 100

cu_func=lambda k: 'k' ## for plotting a black outline of tip circles
su_func=lambda k: 2*(50-30*k.height/myTree.treeHeight) ## black outline in twice as big as tip circle 
zu_func=lambda k: 99
myTree.plotTree(ax_tree,x_attr=x_attr,colour_function=c_func) ## plot branches
myTree.plotPoints(ax_tree,x_attr=x_attr,size_function=s_func,colour_function=c_func,zorder_function=z_func) ## plot circles at tips
myTree.plotPoints(ax_tree,x_attr=x_attr,size_function=su_func,colour_function=cu_func,zorder_function=zu_func) ## plot circles under tips (to give an outline)

"Add genome annotations"
for k in myTree.Objects: ## iterate over branches
    if k.branchType=='leaf':
        for f in range(len(features)): ## iterate over trait keys
            ftype = anc_df.loc[k.numName][features[f]]
            c='midnightblue' if ftype==1 else 'darkgrey' # block mapper
            #c='darkorange' if k.traits[features[f]]=='1' else 'steelblue' # block mapper
            lineage=plt.Rectangle((f,k.y-0.5),1,1,facecolor=c,edgecolor='none') ## rectangle with height and width 1, at y position of tip and at the index of the key
            ax_genome.add_patch(lineage) ## add coloured rectangle to plot
ax_genome.set_xticks(np.arange(0.5,len(features)+0.5))
#clean_feature_names = [n.replace('_',':') for n in features]
clean_feature_names = [r'nsp6:$\Delta$9','ORF9:TRS-3']
ax_genome.set_xticklabels(clean_feature_names,rotation=90)
[ax_genome.axvline(x,color='w') for x in range(len(features))]

"Add legend for lineages"
import matplotlib.patches as mpatches
handles = [mpatches.Patch(color=cmap[i], label=int2Label[i]) for i in range(num_states)]
legend1 = ax_tree.legend(handles=handles,prop={'size': 24},loc='upper left',bbox_to_anchor=(-0.36, 1.)) #was 1.32

"Add month labels as xticks"
step_freq = 1/12
xticks = np.arange(2020.25,absolute_time,step_freq)
ax_tree.set_xticks(xticks)
labels = ['Apr','May','June','July','Aug','Sep','Oct','Nov','Dec','Jan','Feb','March']
ax_tree.set_xticklabels(labels, fontsize=24) #rotation='vertical'
ax_tree.set_xlabel('Month', fontsize=24)

"Or set times"
ats=[k.absoluteTime for k in myTree.Objects]
fr=0.05
ax_tree.set_xlim(min(ats)-fr,max(ats)+fr)

ax_genome.set_xlim(0,len(features))
ax_tree.set_ylim(-5,myTree.ySpan+5)

"Turn axis spines invisible"
[ax_tree.spines[loc].set_visible(False) for loc in ['top','right','left']] ## no axes
[ax_genome.spines[loc].set_visible(False) for loc in ['top','right','left','bottom']] ## no axes

ax_tree.tick_params(axis='x',size=24) ## no labels
ax_genome.tick_params(size=0,labelsize=24)
ax_tree.set_yticklabels([])
ax_genome.xaxis.set_ticks_position('top')
ax_tree.grid(axis='x')

#plt.show()
plt.subplots_adjust(left=0.22)
plt.savefig(fig_file, dpi=300)
