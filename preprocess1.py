# -*- coding: utf-8 -*-
import pandas as pd
import vcf,argparse,os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from time import time
import concurrent.futures
import copy as c
import random

# mudule of single chrom region check
class check(object):
    def __init__(self, bed=None):
        self.bed = bed.sort_values(by='start')
        self.bed['length']=self.bed['end']-self.bed['start']
        self.chr_now = bed['chrom'].iloc[0]



def addtitle(path):
    with open(path) as f:
        lines = f.readlines()
        lines.insert(0,title)
    with open(path,'w+') as w:
        w.writelines(lines)

def rawsplit(data):
    l = [5000,50000,500000,5000000]
    data['center'] = (data['end']+data['start'])*0.5
    data['center'] = data['center'].apply(int)
    #type(data['center'][0])
    ns,ne = [],[]
    for i in range(len(data)):
        il = i % 4
        ns.append(int(data['center'][i]-l[il]))
        ne.append(int(data['center'][i]+l[il]))
    data['ns'] = pd.DataFrame(ns)
    data['ne'] = pd.DataFrame(ne)
    return data
   
if __name__ == '__main__': 
    '''
    #gene = pd.read_csv('../gene_region.list',sep='\t',header=None,names=['name','chrom','start','end'])
    #bed = pd.read_csv('chr1_10k.bed',skiprows=1, sep = '\t')
    rawdata=pd.read_csv('gene_region.list',sep='\t',header=None,names=['name','chrom','start','end'])
    ddd1 = rawsplit(rawdata)
    ddd2 = c.deepcopy(rawdata)
    ddd2['start'] = ddd2['start'] + random.randint(5000000,10000000)
    ddd2['end'] = ddd2['end'] + random.randint(5000000,10000000)
    ddd2 = rawsplit(ddd2) 
    ddd=pd.concat((ddd1,ddd2),ignore_index=True)
    ddd.to_csv('No_4.gene.list',index=False,sep='\t')
    '''
    title = '#HG00403:P:L HG00404:M:R\n'
    data = pd.read_csv('No_4_gene.csv',sep=',')
    gene = list(data.groupby('chrom'))
    for i in gene:
        path = '../bedfile/'+i[0]+'.bed'
        i[1].to_csv(path,sep='\t',index=False)
        addtitle(path)

        
    
    
    
    
    
    
    
    
    
    
    
    

