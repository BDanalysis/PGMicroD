# coding:utf-8
import pysam


Dict={} 


def RefFre():
    N = 0
   
    file = open('Data/example.reads1.fq', 'r')
    for line in file:
        N = N + 1
    file.close()
    N = N / 2

    file_out = open('Data/Result/RefFre.txt', 'w')
    file_refno = open('Data/Ref_Align/exist_Species.fasta', 'r')

    ref = ''  

    line = file_refno.readline()
    while line != '':
        if line[0] == '>':
            i = line.strip().split()[0][1:] 

            name = (line.strip().split()[1]).split(';')[-1]  

            line = file_refno.readline()
            ref = line.strip()

            sum = 0  

            file_in = pysam.AlignmentFile('Data/Ref_Align/pre_Density.bam', 'rb')
            for r in file_in.fetch(i, 0, len(ref)):
                sum += 1
            file_in.close()

            Dict[name]=float(sum)/N

        line = file_refno.readline()

    file_refno.close()



    dic1SortList = sorted(Dict.items(), key=lambda x: x[1], reverse=True)

    for x in dic1SortList:
        file_out.write("%-50s%10s" % (x[0], str(x[1])) + '\n')

    file_out.close()



RefFre()  
