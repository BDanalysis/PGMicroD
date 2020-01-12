#coding:utf-8
import pysam

import math
import os


def singleHVR(i):
    HVRx = {}

    file = open('Data/HVRFile/V' + str(i) + '.fuzznuc', 'r')
    line = ' '
    refNo = -1
    while line != '':
        line = file.readline()
        linex = line.strip().split()
        if 'Sequence:' in linex:
            refNo = linex[linex.index('Sequence:') + 1]
            HVRx[refNo] = []

        if 'Start' in linex and 'End' in linex:
            line = file.readline()
            liney = line.strip().split()
            HVRx[refNo].append([eval(liney[0]), eval(liney[1])])
    file.close()
    return HVRx


def init():
    N = 0
    RefNo = []
    HVR = {}  

    file = open('Data/example.reads1.fq', 'r')
    for line in file:
        N = N + 1
    file.close()
    N = N / 2  


    file = open('Data/Ref_Align/TotalRef.fasta', 'r')
    for line in file:
        if line[0] == '>':
            num = line.strip().split()[0]
            RefNo.append(num[1:len(num)])
    file.close()


    for i in range(len(RefNo)):
        HVR[RefNo[i]] = []


    for i in range(2, 9, 1):
        HVRx = singleHVR(i)

        for (k, v) in HVRx.items():
            if k in HVR.keys():
                for (x0, x1) in HVRx[k]:
         
                    l = HVR[k]
                    if len(l) > 0:
                        (y0, y1) = l[-1]
                        if y1 < x0:
                            HVR[k].append([x0, x1])
                    else:
                        HVR[k].append([x0, x1])

    file_HVR = open('Data/Result/HVR.txt', 'w')
    for (k, v) in HVR.items():
        file_HVR.write(str(k) + ' ')
        for (x0, x1) in v:
            file_HVR.write(str(x0) + ',' + str(x1) + ' ')
        file_HVR.write('\n')
    file_HVR.close()

    return N




def CigerTrans(CIGAR):
    l = []

    num = ''  
    flag = True 
    for x in CIGAR:
        if x >= '0' and x <= '9':
            num += x
            flag = True
        else:
            if flag:
                l.append(eval(num))
                num = ''
            flag = False
            l.append(x)

    return l




def ReadRefTrans(read, ref, CIGAR, score, start):
    l = CigerTrans(CIGAR) 
    
    ref_tran = ''  
    read_tran = ''
    score_tran = ''  

    ref_start = start  
    read_start = 1

    for i in range(0, len(l), 2):
        num = l[i]  
        info = l[i + 1]  

        if info == 'M':

            read_tran += read[read_start - 1: read_start + num - 1]
            score_tran += score[read_start - 1: read_start + num - 1]
            read_start = read_start + num  

            ref_tran += ref[ref_start - 1: ref_start + num - 1]
            ref_start = ref_start + num

        elif info == 'S':

            read_tran += read[read_start - 1: read_start + num - 1]
            score_tran += score[read_start - 1: read_start + num - 1]
            read_start = read_start + num

            ref_tran += 5 * ' '
            # ref_start不变

        elif info == 'I':
            read_tran += read[read_start - 1: read_start + num - 1]
            score_tran += score[read_start - 1: read_start + num - 1]
            read_start = read_start + num

            ref_tran += ' '
           

        elif info == 'D':
            read_tran += ' '
            score_tran += ' '
            

            ref_tran += ref[ref_start - 1: ref_start + num - 1]
            ref_start = ref_start + num

        elif info == 'H':
            continue
        else:
            continue

    return [ref_tran, read_tran, score_tran, l]




def HVRlen(position, start, end):
    sumHVR = 1  
    for (x, y) in position:
        sumHVR += y - x + 1

    hvrLen = 0
    space = []  
    space.append(start)

    for (x0, x1) in position:
        if x0 >= start and x0 <= end:
            space.append(x0)

        if x1 >= start and x1 <= end:
            space.append(x1)

    space.append(end)

    if len(space) > 2:  

        
        flag = False  
        for (x0, x1) in position:
            if x0 <= start and start >= x1:
                flag = True
                break

        if flag == True:  
            for i in range(0, len(space), 2):
                if i + 1 < len(space):
                    hvrLen += space[i + 1] - space[i] + 1
        else:
            for i in range(1, len(space), 2):
                if i + 1 < len(space):
                    hvrLen += space[i + 1] - space[i] + 1

    return float(hvrLen + 1) / sumHVR





def SamLikelihood():
    
    file_nomap = open('Data/Result/NoMapRead.txt', 'w')
    file_oldsam = open('Data/Ref_Align/example.reads.fq.sam', 'r')
    file_sam = open('Data/Ref_Align/example.sam', 'w')
    for sam in file_oldsam:
        if sam[0] == '@':
            file_sam.write(sam)
        else:
            samx = sam.strip().split()
            if samx[2] != '*':
                file_sam.write(sam)
            else:
                file_nomap.write(sam)
    file_nomap.close()
    file_oldsam.close()
    file_sam.close()

    file_sam = open('Data/Ref_Align/example.sam', 'r')
    file_readlike = open('Data/Result/samLike.txt', 'w')
    for sam in file_sam:
        if sam[0] == '@':
            continue

        samx = sam.strip().split()

        
        position = []  
        file_HVR = open('Data/Result/HVR.txt', 'r')
        for refnox in file_HVR:
            refno = refnox.strip().split()[0]

            if refno == samx[2]:
                for HVR in refnox.strip().split()[1:]:
                    [x0, x1] = HVR.split(',')
                    position.append([eval(x0), eval(x1)])
                break
        file_HVR.close()

        
        ref = ''
        file_ref = open('Data/Ref_Align/TotalRef.fasta', 'r')
        line = file_ref.readline()
        while line != '':
            if line[0] == '>':
                linex = line.strip().split()[0]
                refno = linex[1:len(linex)]
                if refno == samx[2]:
                    line = file_ref.readline()
                    ref = line.strip()
                    break
            line = file_ref.readline()
        file_ref.close()

        
        score_phred = samx[10]  
        list = ReadRefTrans(samx[9], ref, samx[5], score_phred, eval(samx[3]))  

        
        w = 0.0  
        align_start = eval(samx[3])  
        info = CigerTrans(samx[5])  
        align_len = 0
        for j in range(0, len(info), 2):
            align_len += info[j]
        align_end = align_start + align_len - 1  

        w = HVRlen(position, align_start, align_end)  

        
        
        ref_tran = list[0]
        read_tran = list[1]
        score_tran = list[2]
        lsigar = list[3]
        multi = 1  

        
        l = [len(ref_tran), len(read_tran), len(score_tran)]

        for j in range(min(l)):  
            if ref_tran[j] != ' ' and read_tran[j] != ' ' and score_tran[j] != ' ':
                q = (math.pow(10, -((ord(score_tran[j]) - 33)) / 10))  

                if ref_tran[j] == read_tran[j]:
                    multi = multi * (1.0 - q)

                else:
                    multi = multi * q / 3.0

        
        gapsum = 1.0  
        for j in range(0, len(lsigar), 2):
            if lsigar[j + 1] in ['S', 'I', 'D', 'H']:
                
                gapsum += math.log(lsigar[j] + 1) / math.log(2)

        if gapsum == 0.0:
            like = multi + w
        else:
            like = (multi + w) / gapsum  

        file_readlike.write(samx[0] + ' ' + str(like) + ' ')
        for i in range(1,len(samx),1):
            file_readlike.write(samx[i]+' ')
        file_readlike.write('\n')

    file_readlike.close()
    file_sam.close()


def LowReadLikelihoodNo(T):

    Recovery_Read = set()
    lowLikeRead = set() 

    file_in = open('Data/Result/samLike.txt', 'r')
    for like in file_in:
        likex = like.strip().split()
        if eval(likex[1]) <= T:
            lowLikeRead.add(likex[0])
            
            C={'D','I','S','H'}
            if len(set(likex[2]) & C)==0:
                Recovery_Read.add(likex[0])
            
            for i in range(12, len(likex), 1):
                if 'XA:Z:' in likex[i]:
                    Recovery_Read.add(likex[0])

    file_in.close()


   
    file_sam = open('Data/Ref_Align/example.sam', 'r')
    file_sam_ = open('Data/Ref_Align/example.reads.sam', 'w')
    for sam in file_sam:
        if sam[0] == '@':
            file_sam_.write(sam)
            continue
        samx = sam.strip().split()
        if samx[0] not in lowLikeRead:
            file_sam_.write(sam)
    file_sam.close()
    file_sam_.close()




    
    file_lowRead = open('Data/Result/Recovery.fastq', 'w')
    for x in Recovery_Read:
        fq1_flag = False

        
        file_fq1 = open('Data/example.reads1.fq', 'r')
        line = ' '
        while line != '':
            line = file_fq1.readline()
            if line.strip() == '@' + x:
                line1 = file_fq1.readline()
                line2 = file_fq1.readline()
                line3 = file_fq1.readline()
                file_lowRead.write('@' + x + '\n' + line1.strip() + '\n' + line2.strip() + '\n' + line3.strip() + '\n')
                fq1_flag = True
                break
        file_fq1.close()

        
        if fq1_flag==False:  
            file_fq2 = open('Data/example.reads2.fq', 'r')
            line = ' '
            while line != '':
                line = file_fq2.readline()
                if line.strip() == '@' + x:
                    line1 = file_fq2.readline()
                    line2 = file_fq2.readline()
                    line3 = file_fq2.readline()
                    file_lowRead.write('@' + x + '\n' + line1.strip() + '\n' + line2.strip() + '\n' + line3.strip() + '\n')
                    break
            file_fq2.close()

    file_lowRead.close()




init()  
SamLikelihood()  
T = 0.44
LowReadLikelihoodNo(T)  
