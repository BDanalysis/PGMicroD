#coding:utf-8

import pysam

def ExtractRef():
    S = set()
    file_in = open('Data/Ref_Align/example.reads.fq.sam', 'r')
    for line in file_in:
        if line[0] == '@':
            continue
        else:
            linex = line.strip().split()
            if linex[2] != '*':
                S.add(linex[2])
    file_in.close()


    # 将比对上的ref作为新ref库，存为TotalRef.fasta
    file=pysam.FastaFile('Data/Ref_Align/TotalRef1.fasta')
    file_out=open('Data/Ref_Align/TotalRef.fasta','w')
    for s in S:
        ref=file.fetch(s)
        file_out.write('>'+s+'\n')
        file_out.write(ref+'\n')
    file_out.close()
    file.close()

ExtractRef()