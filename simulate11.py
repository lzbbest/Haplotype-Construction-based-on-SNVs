# -*- coding: utf-8 -*-
import pandas as pd
import vcf,argparse,os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from time import time
import concurrent.futures
import copy as c
import pysam,gzip
'''
Testing command：
hg19chr1 = '/home/BIOINFO_DATABASE/reference/genome_DNA/Homo_sapiens/hg19/splitbychr_hg19/chr1.fa.gz'
snvIN = '/home/BIOINFO_DATABASE/dbSNP/dbSNP_150/Homo_sapiens/GRCh38/dbSNP_150_hg38.common.prefixChr.vcf.bgz'
/home/tanbowen/HRC_Chinese/
'''

def readBed(bedfile):
    with open(bedfile) as f:
        info = f.readline()[1:].split()
        info1 = info[0].split(':')
        info2 = info[1].split(':')
    bed = pd.read_csv(bedfile,skiprows=1, sep = '\t')    
    #bed.rename(columns={'end':'length'},inplace=True)
    #bed['length'] = bed['length'] - bed['start']
    bed['length'] = bed['end']-bed['start']
    return bed,info1,info2

def readSeq(seqin): #read seq as a string
    handle = gzip.open(seqin, "rt")
    genome = SeqIO.parse(handle, "fasta")
    
    #genome = SeqIO.parse(seqin, "fasta")
    seq_record =next(genome)
    handle.close()
    return seq_record

def report(bedfile,end_pos,rfile,bed):
    with open(bedfile) as f:
        firstrow = f.readline()
    with open(rfile,'w+') as w:
        w.write(firstrow)   
    #bed = pd.read_csv(bedfile,skiprows=1, sep = '\t')
    r_list = [['chrom','sid','start','end','Xmx','P','M']] # Report file record
    r_list.append([str(bed['chrom'][0]),'1','1',str(bed['start'][0]),'2m1','1','1'])
    for i in range(len(bed)):
        temp1,temp2=[],[]
        # changed seg
        temp1.append(bed['chrom'][i])
        temp1.append(str(2*i+2))
        temp1.append(str(bed['start'][i]+1))
        temp1.append(str(bed['end'][i]))
        temp1.append(str(bed['M_Cn'][i]+bed['P_Cn'][i])+'m'+str(bed['M_Cn'][i]))
        temp1.append(str(bed['P_Cn'][i]))
        temp1.append(str(bed['M_Cn'][i]))
        r_list.append(temp1)
        # normal seg
        temp2.append(str(bed['chrom'][i]))
        temp2.append(str(2*i+3))
        temp2.append(str(bed['end'][i]+1))
        if i == len(bed)-1:
            temp2.append(str(end_pos))
        else:
            temp2.append(str(bed['start'][i+1]))
        temp2.append('2m1')
        temp2.append('1')
        temp2.append('1')
        r_list.append(temp2)
    #out = pd.DataFrame(r_list,columns=['chrom','sid','start','end','Xmx','F','M'])
    with open(rfile,'a+') as w:
        for i in range(len(r_list)):
            w.write('\t'.join(r_list[i])+'\n')

def seq_split(seq_record,length):
    seqlist = [] # seqRecord class list
    last_len = len(seq_record) % length
    l = int((len(seq_record)-last_len)/length)
    for i in range(l):
        head,tail = i*length,i*length+length
        temp_seq = seq_record[head:tail]
        seqlist.append([head,tail,temp_seq])# data structure
    lastseq = seq_record[tail:tail+last_len]
    seqlist.append([tail,tail+last_len,lastseq])
    return seqlist 

def exSNVs(seqclass,snvs,info,chr_now): # Extract corresponding SNVs in dbSNP
    seqrecord = seqclass[2]
    head = seqclass[0]
    tail = seqclass[1]
    mseq = seqrecord.seq.tomutable()
    vcf_reader = vcf.Reader(filename=snvs)
    temp_vcf = vcf_reader.fetch(chr_now,head,tail)
    
    sample_name = info[0]
    if info[2] == 'L':
        L_R = int(0)
    else:  ## R
        L_R = int(1)
    pos,base = 0,''
    
    for record in temp_vcf:#vcf_reader:#
        #n+=1
        #print(n)
        if record.is_indel:
            continue
        else:
            r_a = record.genotype(sample_name)['GT'].split('|')[L_R]
            if r_a == '1': #alt genome
                base = str(record.ALT[0]).upper()
            #else:
                #base = str(record.REF[0]).upper()
                pos = record.POS # int
                ###  multi-process
                try:
                    mseq[pos-head-1] = base
                except (KeyError, TypeError, IndexError, BaseException) as error:
                    pass
                #print('error')
            ###
            #mseq[pos] = base
    seqrecord.seq = mseq.toseq()
    return seqrecord # Seqrecord objects


def addseq(seq_record,bed,outfile,info):
    if info[1] == 'M':
        ploidyList = list(bed['M_Cn'])
        H = info[2]
    else:
        ploidyList = list(bed['P_Cn'])
        H = info[2]
    
    '''
    if hap_flag == 'M':
        ploidyList = list(bed['M_Cn'])
        H = '0'
    else:
        ploidyList = list(bed['P_Cn'])
        H = '1'
    '''
    #sid1 = seq_record[:bed['start'][0]]
    sid1 = cutRegion(seq_record,0,bed['start'][0])###
    
    sid1.description = sid1.description +':'+'1-'+str(bed['start'][0])+' SID=1 CID=1 H='+H
    sid1.id = info[1].lower()+'_'+seq_record.description+'_1_1'
    myrecord = [sid1]
    extendlist = []
    for i in range(len(bed)):
        start = int(bed['start'][i])
        end = start + int(bed['length'][i])
        #temp = seq_record[start:end]
        temp = cutRegion(seq_record,start,end)###
        
        tag = seq_record.description+':'+str(start)+'-'+str(end)+' SID='+str(2*i+2)+' CID=1 H='+H
        temp.description = c.deepcopy(tag)
        temp.id = info[1].lower()+'_'+seq_record.description+'_'+str(2*i+2)+'_1'
        # mutable seg
        copy = ploidyList[i]
        if copy ==1:
            myrecord.append(temp)
        elif copy==0:
            pass
            #print('copy=0')
        else:
            myrecord.append(temp)
            CID = 1
            for j in range(1,int(copy)):
                #print(str(j))
                #temp = seq_record[start:end]
                temp = cutRegion(seq_record,start,end)###
                
                CID += 1
                tag = seq_record.description+':'+str(start)+'-'+str(end)+' SID='+str(2*i+2)+' CID='+str(CID)+' H='+H
                #print(tag)
                temp.description = c.deepcopy(tag)
                temp.id = info[1].lower()+'_'+seq_record.description + '_'+str(2*i+2)+'_'+str(CID)
                extendlist.append(temp)
        
        
        # un-mutable seg
        if i != len(bed)-1:   # non last seg
            #temp = seq_record[end:int(bed['start'][i+1])]
            temp = cutRegion(seq_record,end,int(bed['start'][i+1]))###
            
            tag = seq_record.description+':'+str(end)+'-'+str(bed['start'][i+1])+' SID='+str(2*i+3)+' CID=1 H='+H
            temp.description = tag
            temp.id = info[1].lower()+'_'+seq_record.description + '_'+str(2*i+3)+'_1'
            myrecord.append(temp)
        else:   #last seg
            #temp = seq_record[end:]
            temp = cutRegion(seq_record,end,len(seq_record))###
            
            tag = seq_record.description+':'+str(end)+'-'+str(len(seq_record))+' SID='+str(2*i+3)+' CID=1 H='+H
            temp.description = tag
            temp.id = info[1].lower()+'_'+seq_record.description + '_'+str(2*i+3)+'_1'
            myrecord.append(temp)

    myrecord.extend(extendlist)
    SeqIO.write(myrecord,outfile,'fasta')
    return myrecord

def snp_check(snvs,chr_now,head,tail,zsnp):
    vcf_reader = vcf.Reader(filename=snvs)
    temp_vcf = vcf_reader.fetch(chr_now,head,tail)
    n,flag_snp = 0,False
    for record in temp_vcf:
        n += 1
        if n>=zsnp:
            flag_snp = True
            break
    return flag_snp

def uni_check(unifile,chrom,start,end):
    tbx = pysam.TabixFile(unifile)
    sum_len = 0
    for row in tbx.fetch(chrom, start, end, parser=pysam.asTuple()):
        sum_len += int(row[2])-int(row[1]) + 1
    ratio = sum_len/(end-start+1)
    if ratio>0.9:
        return True
    else:
        return False

def cutRegion(seq_record,s,e): #s,e is start,end
    a=c.deepcopy(non_N)
    #sid1 = seq_record[:bed['start'][0]]
    #temp = seq_record[start:end]
    #non_N = non_N[['start','end']][non_N['chr'] == chr_now]
    #use_region = a[['start','end']][(a.start>s) & (a.end<e)]
    if not a[['start','end']][(a.start<s) & (a.end>e)].empty:
        return seq_record[s:e]
    elif not a[['start','end']][(a.start<s) & (a.end>s)].empty:
        head = a[['start','end']][(a.start<s) & (a.end>s)]
        if not head.empty: head['start'].iloc[0] = s
        mid = a[['start','end']][(a.start>s) & (a.end<e)]
        tail = a[['start','end']][(a.start<e) & (a.end>e)]
        if not tail.empty: tail['end'].iloc[0] = e
        temp_region = head.append(mid)
        temp_region = temp_region.append(tail).reset_index(drop=True)
               
        temp_seq = ''
        for i in range(len(temp_region)):
            temp_seq += seq_record[temp_region['start'][i]:temp_region['end'][i]] 
        return temp_seq
    elif a[((a.start<s) & (a.end>s))&((a.start<e) & (a.end>e))].empty and a[(a.start>s) & (a.start<e)].empty:
        return seq_record[0:0]
    else:
        mid = a[['start','end']][(a.start>s) & (a.end<e)]
        tail = a[['start','end']][(a.start<e) & (a.end>e)]
        if not tail.empty: tail['end'].iloc[0] = e
        temp_region = mid.append(tail).reset_index(drop=True)
        
        temp_seq = ''
        for i in range(len(temp_region)):
            temp_seq += seq_record[temp_region['start'][i]:temp_region['end'][i]] 
        return temp_seq

def judge(chrom,start,end):
    pd1 = Ana[['p_endPos','q_startPos']][Ana['#chr'] == chrom]
    if start > pd1['p_endPos'].iloc[0] and end < pd1['q_startPos'].iloc[0]:
        return False
    '''
    pd2 = non_N[['start','end']][non_N['chr'] == chrom]
    if end < pd1['p_endPos'].iloc[0] or start > pd1['q_startPos'].iloc[0]:
        for i in pd2.index:
            if not (end<pd2['start'][i] or start>pd2['end'][i]):
                return False
            else:
                return True
    else:
        return False
    '''
def intervalcheck(bed):
    for i in range(len(bed)-1):
        if bed['end'].iloc[i]>bed['start'].iloc[i+1]:
            if bed['length'].iloc[i]>bed['length'].iloc[i+1]:
                bed['end'].iloc[i] = bed['start'].iloc[i+1]-1
            else:
                bed['start'].iloc[i+1] = bed['end'].iloc[i]+1
    return bed

def bed_check(bed,snvs,zsnp,xlength,unifile):
    bed = intervalcheck(bed)
    dropline=[]
    for i in range(len(bed)):            
        chrom,start = bed['chrom'][i],bed['start'][i]
        end = start + bed['length'][i]
        
        flag_snp = snp_check(snvs,chrom,start,end,zsnp)
        if flag_snp == False: 
            dropline.append(i)
            print('Warning: The segment '+ str(i+1)+' does not have sufficient SNPs')
        if bed['length'][i] < xlength:
            print(bed['length'][i])
            dropline.append(i)
            print('Warning: The segment '+ str(i+1)+' does not have sufficient length')        
        flag_uni = uni_check(unifile,chrom,start,end)
        if flag_uni == False:
            dropline.append(i)
            print('Warning: The segment '+ str(i+1)+' does not have sufficient high confidence region')        
        tag = judge(chr_now,start,end)
        if tag==False: print('Warning: Region error in  '+str(chrom)+' start:'+str(start))
        
    dropline = list(set(dropline))
    new_bed = bed.drop(dropline)
    new_bed = new_bed.reset_index(drop=True)# rebuile index from 0
    return new_bed
    
def cmd():# test
    parser = argparse.ArgumentParser(description = 'Haplotype Data simulation program',epilog = 'Parameter describe ')
    #required parameters
    parser.add_argument('-v','--vcf',type = str,dest='snv',required = True,help = 'Please input your phased VCF file')
    parser.add_argument('-s', '--sequence', type = str,dest='seqin', required = True, help = 'Please input your ref Fasta file')
    parser.add_argument('-o','--origin1',type = str,dest='outfile0p',required = True, help = 'Please input your reference first mutable output file')
    parser.add_argument('-l','--origin2',type = str,dest='outfile0m',required = True, help = 'Please input your reference second mutable output file')
    parser.add_argument('-p','--paternal',type = str,dest='outfile1',required = True, help = 'Please input your first mutable output file')
    parser.add_argument('-m','--maternal',type = str,dest='outfile2',required = True, help = 'Please input your second mutable output file')
    #parser.add_argument('-c','--centromere',type = str,dest='anafile',required = True,help = 'Please input your centromere region file') 
    #parser.add_argument('-n','--nonN',type = str,dest='nonNfile',required = True,help = 'Please input your non_N region file')
    parser.add_argument('-b','--bed',type = str,dest='bedfile',required = True,help = 'Please input your simulation imformation file')
    #parser.add_argument('-f','--forks',type = int,dest='n_jobs',required = True,help = 'Please input how many processor you want to use')
    parser.add_argument('-r','--rfile',type = str,dest='rfile',required = True,help = 'Please input your report file')
    
    #parser.add_argument('-z','--zsnp',type = str,dest='zsnp',required = False,default='5',help = 'The minimum number of snps in a segment')
    #parser.add_argument('-x','--xlength',type = str,dest='xlength',required = False,default='100000',help = 'The minimal length of a segment')
    #optional parameters

    args = parser.parse_args() #
    return args

if __name__ == '__main__':
    
    begin = time()
    #print(str(begin))
    
#    snv = '../CDX_CHB_CHS.chr1.vcf.gz'
    #snv = '200r.vcf.gz'
#    seqin = '../chr1.fa'

#    outfile0p = 'outp_ref_chr1_10k.fa'
#    outfile0m = 'outm_ref_chr1_10k.fa'
#    outfile1 = 'out1_chr1_10k.fa'
#    outfile2 = 'out2_chr1_10k.fa'

    anafile = '../hg19.chrForAna.info.tsv'
    nonNfile = '../nonN_region'
    n_jobs = 25
#    rfile = 'Report.csv'
    zsnp = 5
    xlength = 9997
    unifile = '../wgEncodeCrgMapabilityAlign100mer.ScoreEq1.merged.bed.bgz'
    
    args = cmd()
    snv = args.snv
    seqin = args.seqin
    outfile0p = args.outfile0p
    outfile0m = args.outfile0m
    outfile1 = args.outfile1
    outfile2 = args.outfile2
    bedfile = args.bedfile
    rfile = args.rfile

    seg_length = 500000  # size const

       
    print('Program started......')

    Ana = pd.read_csv(anafile, sep = '\t',usecols=[0,4,5])
    non_N = pd.read_csv('../nonN_region', sep = '\t',skiprows=1,names=['chr','start','end'])
    bed,info1,info2 = readBed(bedfile)
    global chr_now
    chr_now = bed['chrom'][0]
    non_N = non_N[['start','end']][non_N['chr'] == chr_now]

    bed = bed_check(bed,snv,zsnp,xlength,unifile)
    print(chr_now+' Region checking has been completed......')
       
    seq_record = readSeq(seqin)
    report(bedfile,len(seq_record),rfile,bed)#    #seqlist = seq_split(seq_record,seg_length)
    seqlist = seq_split(seq_record,seg_length)
    len_seqlist = len(seqlist)
    print('5%=>===================100%'+'    Time consume: '+str(time()-begin))
  
    #vcf_reader = vcf.Reader(filename=snv)
    # sequence 1

    seq_record_seq = ''
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_jobs) as executor:
        for x in executor.map(exSNVs,seqlist,[snv]*len_seqlist,[info1]*len_seqlist,[chr_now]*len_seqlist):
            #seq_record_list.append(x) # x is SeqRecord
            seq_record_seq += x.seq
    seq1 = SeqRecord(seq_record_seq, id=seq_record.id,description=seq_record.description,name=seq_record.name)
    
    seq11 = c.deepcopy(seq1)
    seq11.id = info1[1].lower()+'_'+seq11.id
    seq11.description = seq11.description+' H = '+info1[2]
    SeqIO.write(seq11,outfile0p,'fasta')
    print('40%========>============100%'+'    Time consume: '+str(time()-begin))

    #seq1 = seq_record
    #hap_flag = info1[1]
    myrecord = addseq(seq1,bed,outfile1,info1)
    print('50%==========>==========100%'+'    Time consume: '+str(time()-begin))
    del seq1,seq11,myrecord

    #sequence 2
    seq_record_seq = ''
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_jobs) as executor:
        for x in executor.map(exSNVs,seqlist,[snv]*len_seqlist,[info2]*len_seqlist,[chr_now]*len_seqlist):
            #seq_record_list.append(x) # x is SeqRecord
            seq_record_seq += x.seq
    seq2 = SeqRecord(seq_record_seq, id=seq_record.id,description=seq_record.description,name=seq_record.name)

    seq22 = c.deepcopy(seq2)
    seq22.id = info2[1].lower()+'_'+seq22.id
    seq22.description = seq22.description+' H = '+info2[2]    
    
    SeqIO.write(seq22,outfile0m,'fasta')
    print('90%==================>==100%'+'    Time consume: '+str(time()-begin))
    #hap_flag = info2[1]
    myrecord = addseq(seq2,bed,outfile2,info2)
    del seq2,seq22,myrecord

    print('100%====================>100%')
    end = time()
    print('Program completed     Runing time: '+str(end-begin))
