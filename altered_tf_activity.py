# -*- coding: utf-8 -*-
"""
Created by Gaurang Mahajan on June 21 2021.

"""

import streamlit as st
import numpy as np
import seaborn as sns
import matplotlib
from matplotlib.figure import Figure
import pandas as pd
from collections import Counter,OrderedDict
from scipy import stats
from operator import itemgetter
import networkx as nx

####################################
network_file = r'.\RegulonDB_pairs.txt'
geneexp_file = r'.\ph8p7_degenes_fdr1e-3_fc2.txt'
####################################

def targets(tf,g):
    return [s for s in g.successors(tf)]

def find_best_tf(tfs,bkg,sig,trn):
    tflist=[]
    for tf in tfs:
        a=len(np.intersect1d(targets(tf,trn),sig))
        b=len(sig)-a
        c=len(np.intersect1d(targets(tf,trn),bkg))-a
        d=len(bkg)-a-b-c
        o,p=stats.fisher_exact([[a,b],[c,d]],alternative='greater')
        tflist.append((tf,p))
    return sorted(tflist,key=itemgetter(1))[0][0]

def combined_pval(tfsol,allgenes,siggenes,trn):
    l = [set(targets(tf,trn)) for tf in tfsol]
    targetset = list(frozenset().union(*l))
    a=len(np.intersect1d(targetset,siggenes))
    b=len(siggenes)-a
    c=len(np.intersect1d(targetset,allgenes))-a
    d=len(allgenes)-a-b-c
    o,p=stats.fisher_exact([[a,b],[c,d]],alternative='greater')   
    return np.log10(p)
    
#st.write("""Hello *World*!!!  """)

plot_types = ("line","scatter")

with st.sidebar:
    st.header("Input parameters")
    with st.form(key="grid_reset"):
        n_cutoff = st.slider("Cutoff on TF network degree", 2, 20)
        p_cutoff = st.number_input("-Log P-val cutoff\n on indiv. TF enrichment", 2, 6)
        chart_type = st.selectbox("Choose chart type", plot_types)
        st.form_submit_button(label="Score TFs on input dataset")
    
    with st.beta_expander("About this app"):
        st.markdown("Computes differential TF activity, given a list of gene expression changes and \
        precompiled TF-gene regulatory network (implements method proposed in \
        [Mahajan & Mande, 2015](https://doi.org/10.1371/journal.pone.0142147))")

with st.beta_container():
    st.header("**TF scores on input DEG dataset**")


e1=0
e2=0
####################################
g = nx.DiGraph()
try:
    f = open(network_file,'r')
    for line in f.readlines():
        line=line.strip('\n').split('\t')
        g.add_edge(line[0],line[1])
    f.close()
except IOError as e1:
    with st.beta_container():
        st.write("I/O error({0}): {1}".format(e1.errno, e1))

regs = [r for r in g.nodes() if g.out_degree(r) >= n_cutoff]
with st.beta_container(): 
    st.write('Pre-compiled regulatory network has',g.number_of_nodes(),'genes,',\
      g.number_of_edges(),'directed links, and',len(regs),'TFs')
####################################

####################################
siggenes=[]
allgenes=[]
try:
    f = open(geneexp_file,'r')
    for line in f.readlines():
        line=line.strip('\n').split('\t')
        allgenes.extend(line[0].split(';'))
        if int(line[1])==1: siggenes.extend(line[0].split(';')) 
    f.close()
except IOError as e2:
    with st.beta_container():
        st.write("I/O error({0}): {1}".format(e2.errno, e2))
    
siggenes=list(set(siggenes))
allgenes=list(set(allgenes))
tflist = [tf for tf in regs if len(np.intersect1d(siggenes,targets(tf,g)))>0]
with st.beta_container():
    st.write(len(siggenes),'genes out of a total of', len(allgenes),'are labeled as significantly diff. expressed (DEG).')
    st.write(len(tflist),'TFs have non-zero overlap with DEGs.')
####################################

if e1!=0 or e2!=0 or len(tflist)==0:
    with st.beta_container(): st.write("\nError: Can't run test.")
    
if ((e1==0) and (e2==0) and (len(tflist)>0)):

    ####################################
    sigtfs={}
    for tf in tflist[:]:
       a = len(np.intersect1d(siggenes,targets(tf,g)))
       b = len(siggenes)-a
       c = len(np.intersect1d(allgenes,targets(tf,g)))-a
       d = len(allgenes)-a-b-c
       o,p = stats.fisher_exact([[a,b],[c,d]],alternative='greater')
       sigtfs[tf] = np.min([1,p*len(tflist)])
    with st.beta_container(): 
        st.write(len([tf for tf in sigtfs if sigtfs[tf] < 10**(-p_cutoff)]),\
        'TFs are individually associated with DEGs at p<',10**(-p_cutoff),'level:')
        
    df1 = pd.DataFrame(sorted([[tf,len(targets(tf,g)),len(np.intersect1d(siggenes,targets(tf,g))),sigtfs[tf]] for tf in sigtfs \
    if sigtfs[tf] < 10**(-p_cutoff)],key=itemgetter(-1)),\
    columns=['TF','Total # of targets','# DEG targets','Indiv. p-val\n(FET; adjusted)'])
    with st.beta_container(): df1
    ####################################

    with st.beta_container(): 
        st.write(""" TF subnetwork output: """)
        
    ####################################
    tfset=tflist[:]
    tf = find_best_tf(tfset,allgenes,siggenes,g)
    tfsol_old=[tf]
    tfset.remove(tf)
    tfset_rem=tfset[:]
    old_pval=combined_pval(tfsol_old,allgenes,siggenes,g)
    if len(tfset)>0:
       while len(tfset_rem[:])>0: 
             tflist=[]
             for tf in tfset_rem:
                 tfsol_temp = tfsol_old[:]
                 tfsol_temp.append(tf)
                 tflist.append((tf,combined_pval(tfsol_temp,allgenes,siggenes,g)))
             (tf,p)=sorted(tflist,key=itemgetter(1))[0]
             if p >= old_pval:
                break
             else:  
                tfsol_new=tfsol_old[:]
                tfsol_new.append(tf)
                old_pval = p
                tfsol_old = tfsol_new[:]
                tfset_rem = [r for r in tfset_rem[:] if r!=tf]
    tfsol_approx = tfsol_old[:]
    pval = combined_pval(tfsol_approx,allgenes,siggenes,g)
    
    df2 = pd.DataFrame(sorted([[tf,len(targets(tf,g)),len(np.intersect1d(siggenes,targets(tf,g))),sigtfs[tf]] for tf in tfsol_approx],\
    key=itemgetter(-1)),columns=['TF','Total # of targets','# DEG targets','Indiv. p-val\n(FET; adjusted)'])
    with st.beta_container(): 
        st.write('Optimal (log10) combined p-value =',round(pval,1))
        df2
