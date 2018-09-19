# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import vcf,argparse,os
import pysam
import copy as c
from Bio import bgzf
from ACO4 import MILPsolver

# judge male of famale
def is_male():
    try:
        temp_r = vcf_reader.fetch('chrY')
        return True
    except:
        return False

# extract segment information from 'CNV' or other file
def seg_extract(input_cnv):
    ########   different file use different reading format
    #name = ['chrom','start','end','pho']
    #cnv = pd.read_csv(input_cnv, skiprows=1,usecols=[1,2,3,14],names=name)
    cnv = pd.read_csv(input_cnv, skiprows=0,sep='\t',usecols=[0,1,2,3,4,5])
    #######
    if is_male():
        cnv = cnv[(cnv.chrom != 'chrX') & (cnv.chrom !='chrY')]
    #normal = cnv[(cnv.P_Cn==1)&(cnv.M_Cn==1)]
    return cnv

# obtain "avg_depth_tumor" and C_i in each seg
class coverage(object):
    def __init__(self,filename,cnv):
        self.bamfile = pysam.AlignmentFile(filename,'rb')
        self.cnvfile = cnv
    
    def getDepth(self):
        avg_depth_tumor = []
        for i in range(len(self.cnvfile)):
            s = self.cnvfile['start'].iloc[i]
            e = self.cnvfile['end'].iloc[i]
            chrom = self.cnvfile['chrom'].iloc[i]
            d = self.bamfile.count_coverage(chrom,s,e)
            t_coverage = sum(d[0]+d[1]+d[2]+d[3])/(e-s)
            avg_depth_tumor.append(t_coverage)
            print(t_coverage)
        return np.array(avg_depth_tumor)
    
    # need to get ploidy of normal segs(P_cn=1 & M_cn=1)
    def getCi(self,avg_depth_tumor):
        normal = self.cnvfile[(self.cnvfile.P_Cn==1)&(self.cnvfile.M_Cn==1)]
        ix = list(normal.index)
        normal_d = avg_depth_tumor[ix]
        d_hat = np.mean(normal_d)
        print(d_hat)
        ci = 2*avg_depth_tumor/d_hat
        return np.array(ci)



#read & write "INFO" and extract DP4 in one segment
def rwDP4(chr_now,start,end):
    print('Processing '+chr_now+' start: '+str(start)+'  end: '+str(end))
    M,m,tl = 0.0,0.0,[]
    n_SNP,n_het,psnp,phet = 0,0,0,0
    
    temp_r = vcf_reader.fetch(chr_now,start,end)
        
    for record in temp_r:
        t = record.INFO['DP4'] #type(t)=list
        
        if not record.is_snp:
            continue
        elif not record.genotype(sample_name).is_het:
            n_SNP += 1
            psnp += 1
            vcf_writer.write_record(record)
            continue
        elif t[2]<2 or t[3]<2 or t[1]<2 or t[0]<2:
            n_SNP += 1
            n_het += 1
            continue
        else:
            n_SNP += 1
            n_het += 1
            psnp += 1
            phet += 1
            vcf_writer.write_record(record)

            M = t[0] + t[1]
            m = t[2] + t[3]
            if m > M and M !=0:
                tl.append(m/M)
            else:
                tl.append(M/m)
        #print(t)
    return tl,n_SNP,n_het,psnp,phet

# SNVs checking and aggregation
def ksegment(cnv):
    snps,hets,phasedSNP,phasedhet = [],[],[],[]
    dropline,dp4 = [],[]
    for i in range(len(cnv)):
        start = cnv.start[i]
        end = cnv.end[i]
        chr_now = cnv.chrom[i]
        #pho = cnv.pho[i]
        #search corresponding segment in snv
        try:
            ratio,n_SNP,n_het,psnp,phet = rwDP4(chr_now,start,end)
            snps.append(n_SNP)
            hets.append(n_het)
            phasedSNP.append(psnp)
            phasedhet.append(phet)
            ratio = np.array(ratio)
        except:
            snps.append(0)
            hets.append(0)
            phasedSNP.append(0)
            phasedhet.append(0)
            ratio = 0        
        if len(ratio)<3:
            dropline.append(i)
        else:
            dp4.append(ratio)

    snps = pd.DataFrame(snps,columns=['n_SNP'])
    hets = pd.DataFrame(hets,columns=['n_het'])
    phasedSNP = pd.DataFrame(phasedSNP,columns=['n_phased'])
    phasedhet = pd.DataFrame(phasedhet,columns=['n_phased_het'])
    
    #cnv_out = c.deepcopy(cnv.drop('pho',axis=1))
    cnv_out = c.deepcopy(cnv)
    #cnv_out['ploidy'] = rm['maternal_copy'].apply(int) + rM['paternal_copy'].apply(int)
    cnv_out['length'] = cnv_out['end']-cnv_out['start'] + 1
    cnv_out['ploidy'] = None
    cnv_out['paternal_copy'] = None
    cnv_out['maternal_copy'] = None
    cnv_out['avg_depth_tumor'] = None
    cnv_out['n_SNP'] = snps['n_SNP']
    cnv_out['n_het'] = hets['n_het']
    cnv_out['n_phased'] = phasedSNP['n_phased']
    cnv_out['n_phased_het'] = phasedhet['n_phased_het']
    
    cnv_out2 = cnv_out.ix[dropline,[0,1,2]]
    cnv_out2['n_SNP'] = 'not_pass'
    cnv_out1 = cnv_out.drop(dropline)
    cnv_out1 = cnv_out1.reset_index(drop=True)
    
    return cnv_out1,cnv_out2,dp4

# optimization
def optimal(pho,dp4,c_i):
    pho = np.array(pho)
    '''
    print('Optimization')
    print('The dp4 length is :')
    print(len(dp4))
    print(c_i)
    print(len(c_i))
    a=pd.DataFrame(dp4)
    a.to_csv('dp4.csv',index=None,header=None,sep='\t')
    np.savetxt('c_i.txt',c_i)
    #a=a.reshape(-1,1) 9x1 mat
    '''
    solver = MILPsolver(dp4,c_i)
    M,m = solver.optimal()
    return M,m

# segments write to disk
def segwrite(cnv_out1,cnv_out2,out_seg,out_error,M,m,avg_depth_tumor):
    cnv_out1 = cnv_out1.drop('pho',axis=1)
    cnv_out1['paternal_copy'] = pd.DataFrame(M)
    cnv_out1['maternal_copy'] = pd.DataFrame(m)
    cnv_out1['ploidy'] = cnv_out1['maternal_copy'].apply(int) + cnv_out1['paternal_copy'].apply(int)
    cnv_out1['avg_depth_tumor'] = pd.DataFrame(avg_depth_tumor)
    
    cnv_out1.to_csv(out_seg,sep='\t',index=False)
    cnv_out2.to_csv(out_error,sep='\t',index=False)

def cmd():# test
    parser = argparse.ArgumentParser(description = 'Haplotype Data simulation program',epilog = 'Parameter describe ')
    #required parameters
    parser.add_argument('-v','--vcf',type = str,dest='snv',required = True,help = 'Please input your phased VCF file')
    parser.add_argument('-s','--sample',type = str,dest='sample',required = True,help = 'Please input your sample name')
    
if __name__ == '__main__':
    #input_cnv,input_snvT,chrom,start,end,outfile = cmd()
    
    #########  input and output
    global sample_name
    sample_name = 'No_3'#'70T'
    base_path = sample_name
    if not os.path.exists(base_path):
        os.makedirs(base_path)

    bamfile = '/home/liuzhibo/out_simulate/No_3_combined/aln/No_3.CNV.sort.bam'
    #input_cnv = 'cnv/70t_Copynumbers.csv'
    #input_snvT = 'snv/70T.vcf.gz'
    input_cnv = 'cnv/No_3_Copynumbers.csv'
    input_snvT = '/home/liuzhibo/out_simulate/No_3_combined/snv/No_3.CNV.vcf.gz'
    out_vcf = os.path.join(base_path,'SNP.phased.vcf.gz')
    out_seg = os.path.join(base_path,'segments.list')
    out_error = os.path.join(base_path,'unphased_segments.list') 
    #########
    
    global vcf_reader,vcf_writer    
    vcf_reader = vcf.Reader(filename=input_snvT)
    vcf_writer = vcf.Writer(bgzf.open(out_vcf, 'w'), vcf_reader)

    cnv = seg_extract(input_cnv)
    cnv_out1,cnv_out2,dp4 = ksegment(cnv)
    
    cov = coverage(bamfile,cnv_out1)
    avg_depth_tumor = cov.getDepth() # numpy.array
    c_i = cov.getCi(avg_depth_tumor)
    
    M,m = optimal(cnv_out1['pho'],dp4,c_i)
    segwrite(cnv_out1,cnv_out2,out_seg,out_error,M,m,avg_depth_tumor)
   
    vcf_writer.flush()
    vcf_writer.close()


    