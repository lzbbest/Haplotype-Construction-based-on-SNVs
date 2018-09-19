# -*- coding: utf-8 -*-
import pandas as pd
import vcf,argparse,os,re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from time import time
from glob import glob


if __name__ == '__main__': 
    #gene = pd.read_csv('../gene_region.list',sep='\t',header=None,names=['name','chrom','start','end'])
    #bed = pd.read_csv('chr1_10k.bed',skiprows=1, sep = '\t')
    bedfiles = glob('../bedfile/*.bed')
    snv_path = '/home/tanbowen/HRC_Chinese/'
    chrom_path = '/home/BIOINFO_DATABASE/reference/genome_DNA/Homo_sapiens/hg19/splitbychr_hg19/'
    
    out_path = '/home/liuzhibo/out_simulate/No_4'
    out_base = 'user_specified_prefix.'
    for b in bedfiles:
        #######  Input  #######
        parent_path,name = os.path.split(b)
        name = name.replace('.bed','')
        print('Simulating '+name+'......')
        #snv_name = 'CDX_CHB_CHS.' + name + '.vcf.gz'
        #snv_vcf = os.path.join(snv_path,snv_name)
        #print(snv_vcf)
        chrom_name = name + '.fa.gz'
        chrom_seq = os.path.join(chrom_path,chrom_name)
        ####### Output  #######
        report = os.path.join(out_path,name+'_No4_Report.csv')
        outfile0p = os.path.join(out_path,'No_4.paternal.'+name+'.normal.fasta')
        outfile0m = os.path.join(out_path,'No_4.maternal.'+name+'.normal.fasta')
        outfile1 = os.path.join(out_path,'No_4.paternal.'+name+'.CNV.fasta')
        outfile2 = os.path.join(out_path,'No_4.maternal.'+name+'.CNV.fasta')
        
        #command = 'python simulate11.py -v '+snv_vcf+' -s '+chrom_seq+' -o '+outfile0p+' -l '+outfile0m+' -p '+outfile1+' -m '+outfile2+ ' -b '+b+' -r '+report
        snv_vcf = '/home/tanbowen/HRC_Chinese/CDX_CHB_CHS.vcf.gz'
        command = 'python simulate11.py -v '+snv_vcf+' -s '+chrom_seq+' -o '+outfile0p+' -l '+outfile0m+' -p '+outfile1+' -m '+outfile2+ ' -b '+b+' -r '+report
        os.system(command)

        print('Simulating End...')
        print('\n\n')
    
    
#nohup python simulate11.py -v ../CDX_CHB_CHS.chr1.vcf.gz -s /home/BIOINFO_DATABASE/reference/genome_DNA/Homo_sapiens/hg19/splitbychr_hg19/chr1.fa.gz -o testfile.fasta -b chr1.bed -r report.csv
    
    
    
    
    
    
    
    
    
    

