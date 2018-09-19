# -*- coding: utf-8 -*-
import pandas as pd

class check(object):
    def __init__(self, bed=None):
        self.bed = bed
        print(bed)
    
    def pa(self,ag):
        print(self.name)
        print(ag)
        

    
if __name__ == '__main__': 
    gene = pd.read_csv('gene_region.list',sep='\t',header=None,names=['name','chrom','start','end'])
    s = check('aaa')
