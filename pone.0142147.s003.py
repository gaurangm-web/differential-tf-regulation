#!/usr/bin/env python

## Usage: python S1_File.py
##
## Location of network to be specified by "network_file"
## Location of DEG file to be specified by "DEG_file"
## network_file: Two-column tab-delimited text file listing regulatory links (TF --> target gene).
## DEG_file: Two-column tab-delimited text file listing all genes with expression info; binary label (0/1)
##           specifies whether the gene is significantly differentially expressed (1) or not (0).
## Note: The gene symbols in the two files must be consistent.
##
## The script finds individually over-represented TFs at Bonferroni p < 0.05,
## and runs approximate methods A, B and C before listing out the optimal TF subset identified by method C.
##
## Note: The script makes use of the python package NetworkX


####################################
network_file = r'network_file.txt'
DEG_file = r'DEG_file.txt'
####################################


import numpy as np
from collections import Counter,OrderedDict
from scipy import stats
from operator import itemgetter
import networkx as nx

def find_best_tf(tfs,bkg,sig,trn):
    tflist=[]
    for tf in tfs:
        a=len(np.intersect1d(trn.successors(tf),sig))
        b=len(sig)-a
        c=len(np.intersect1d(trn.successors(tf),bkg))-a
        d=len(bkg)-a-b-c
        o,p=stats.fisher_exact([[a,b],[c,d]],alternative='greater')
        tflist.append((tf,p))
    return sorted(tflist,key=itemgetter(1))[0][0]

def combined_pval(tfsol,allgenes,siggenes,trn):
    l = [set(trn.successors(tf)) for tf in tfsol]
    targetset = list(frozenset().union(*l))
    a=len(np.intersect1d(targetset,siggenes))
    b=len(siggenes)-a
    c=len(np.intersect1d(targetset,allgenes))-a
    d=len(allgenes)-a-b-c
    o,p=stats.fisher_exact([[a,b],[c,d]],alternative='greater')   
    return np.log10(p)


e1=0
e2=0
####################################
print 'Reading regulatory network information...'
g = nx.DiGraph()
try:
    f = open(network_file,'r')
    for line in f.readlines():
        line=line.strip('\n').split('\t')
        g.add_edge(line[0],line[1])
    f.close()
except IOError as e1:
    print "I/O error({0}): {1}".format(e1.errno, e1)

regs = [r for r in g.nodes() if g.out_degree(r)>0]
print 'Regulatory network has',g.number_of_nodes(),'genes,',\
      g.number_of_edges(),'links and',len(regs),'TFs.'
print 'TFs:\n',','.join(regs)
####################################

####################################
print '\nReading differential expression information...'
siggenes=[]
allgenes=[]
try:
    f = open(DEG_file,'r')
    for line in f.readlines():
        line=line.strip('\n').split('\t')
        allgenes.extend(line[0].split(';'))
        if int(line[1])==1: siggenes.extend(line[0].split(';')) 
    f.close()
except IOError as e2:
    print "I/O error({0}): {1}".format(e2.errno, e2)
    
siggenes=list(set(siggenes))
allgenes=list(set(allgenes))
print len(siggenes),'genes out of a total of', len(allgenes),'are DEGs.'
tflist = [tf for tf in regs if len(np.intersect1d(siggenes,g.successors(tf)))>0]
print '\n',len(tflist),'TFs have non-zero overlap with DEGs.'
####################################



if e1!=0 or e2!=0 or len(tflist)==0:

    print "\nError: Can't run test."

else:
    
    ####################################
    print '\nTesting target sets of individual TFs for association with DEGs...'
    sigtfs={}
    for tf in tflist[:]:
       a = len(np.intersect1d(siggenes,g.successors(tf)))
       b = len(siggenes)-a
       c = len(np.intersect1d(allgenes,g.successors(tf)))-a
       d = len(allgenes)-a-b-c
       o,p = stats.fisher_exact([[a,b],[c,d]],alternative='greater')
       sigtfs[tf] = np.min([1,p*len(tflist)])
    print len([tf for tf in sigtfs if sigtfs[tf]<0.05]),\
          'TFs are individually associated with DEGs at p<0.05 level'
    ####################################

    ####################################
    print '\nSequential search for TF combination (Method A)...'
    tfset=tflist[:]
    orderedlist=sorted([(tf,combined_pval([tf],allgenes,siggenes,g)) for tf in tfset],key=itemgetter(1))
    tfsol=[orderedlist[0][0]]
    old_pval=combined_pval(tfsol,allgenes,siggenes,g)
    tfsol_temp=tfsol[:]
    if len(orderedlist)>1:
      for tf,p in orderedlist[1:]:
          tfsol_temp.append(tf)
          if combined_pval(tfsol_temp,allgenes,siggenes,g)>old_pval:
             break
          else:
             tfsol.append(tf) 
             old_pval=combined_pval(tfsol[:],allgenes,siggenes,g)
    tfsol_approx1 = tfsol[:]
    pval1=combined_pval(tfsol_approx1,allgenes,siggenes,g)
    print 'Log combined p-value =',pval1
    print '# TFs =',len(tfsol_approx1)
    #####################################

    #####################################
    print '\nSequential search for TF combination (Method B)...'
    tfset=tflist[:]
    tf = find_best_tf(tfset,allgenes,siggenes,g)
    tfsol_old=[tf]
    tfset.remove(tf)
    tfset_rem=tfset[:]
    old_pval=combined_pval(tfsol_old,allgenes,siggenes,g)
    removed_nodes = np.intersect1d(siggenes,g.successors(tf))
    bkg = [gene for gene in allgenes if gene not in removed_nodes]
    sig = [gene for gene in siggenes if gene not in removed_nodes]
    if len(tfset)>0:
       while len(tfset_rem[:])>0:
             tf = find_best_tf(tfset_rem,bkg,sig,g)
             tfsol_new=tfsol_old[:]
             tfsol_new.append(tf)
             new_pval = combined_pval(tfsol_new,allgenes,siggenes,g)   
             if new_pval > old_pval:
                break
             else:
                old_pval = new_pval 
                tfsol_old = tfsol_new[:]
                tfset_rem = [r for r in tfset_rem[:] if r!=tf]
                removed_nodes = np.intersect1d(sig,g.successors(tf))
                bkg = [gene for gene in bkg[:] if gene not in removed_nodes]
                sig = [gene for gene in sig[:] if gene not in removed_nodes]
    tfsol_approx2=tfsol_old[:]
    pval2=combined_pval(tfsol_approx2,allgenes,siggenes,g)
    print 'Log combined p-value =',pval2
    print '# TFs =',len(tfsol_approx2)
    ######################################

    ######################################
    print '\nSequential search for TF combination (Method C)...'
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
    tfsol_approx3 = tfsol_old[:]
    pval3=combined_pval(tfsol_approx3,allgenes,siggenes,g)
    print 'Log combined p-value =',pval3
    print '# TFs =',len(tfsol_approx3)
    #######################################

    #######################################
    print 'TFs identified by Method C:\n---------------------------------------\n TF\tIndividual p-value (corrected)\n---------------------------------------'
    for tf in tfsol_approx3: print tf,'\t',sigtfs[tf]
    print '---------------------------------------'
    #######################################
