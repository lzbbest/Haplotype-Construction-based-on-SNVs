#! /bin/sh

for File in /home/liuzhibo/out_simulate/No_2/*.fasta
do
    prefix=${File%%.fasta*}
    one=$prefix'_1.fastq.gz'
    two=$prefix'_2.fastq.gz'
    #$bioto_wgsim_0_3_2 -N 1000000 -1 100 -2 100 -r 0 -R 0 -X 0 $File $one $two
    /home/BIOINFO_TOOLS/mutation_tools/SamTools/SamTools-1.3/bin/wgsim -N 1000000 -1 100 -2 100 -r 0 -R 0 -X 0 $File $one $two
	echo $one
    echo $two
    echo $File
done







