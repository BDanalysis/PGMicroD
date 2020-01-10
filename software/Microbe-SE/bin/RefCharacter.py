#coding:utf-8
import math
import os
import pysam
from scipy.stats import gamma



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





def CoverageSpecial(refNo):
    ref = ''  
    file = open('Data/Ref_Align/TotalRef.fasta', 'r')
    line = file.readline()
    while line != '':
        if line[0] == '>':
            linex = line.strip().split()[0]
            refno = linex[1:len(linex)]
            if refno == refNo:
                line = file.readline()
                ref = line.strip()
                break
        line = file.readline()
    file.close()

    Cov = [0] * len(ref) 

   
    for i in range(len(Cov)):
        file_bam = pysam.AlignmentFile('Data/Ref_Align/example.reads.bam','rb')
        for r in file_bam.fetch(refNo,i,i+1):
            Cov[i] = Cov[i]+1
        file_bam.close()

    return Cov




def RefCov():
    #file_gamma = open('Data/Result/RefGamma.txt','w') 
    file_cov = open('Data/Result/RefCoverage.txt','w')

    file_noapper = open('Data/Result/NoApperRef.txt','w')

    file_refno = open('Data/Result/HVR.txt', 'r')
    for line in file_refno:
        refno = line.strip().split()[0]

        cov = CoverageSpecial(refno)


        if sum(cov)==0: 
            file_noapper.write(refno+'\n')
            continue


        file_cov.write(refno + ' ')
        for i in range(0,len(cov),1):
            file_cov.write(str(cov[i])+' ')
        file_cov.write('\n')

    file_refno.close()
    #file_gamma.close()
    file_cov.close()
    file_noapper.close()





def RefGap():

    file_gapscore = open('Data/Result/GapScore.txt','w')

    file_cov = open('Data/Result/RefCoverage.txt','r')
    for cov in file_cov:
        Coverage = []
        for k in range(1,len(cov.strip().split()),1):
            Coverage.append(eval(cov.strip().split()[k]))

        Coverage.append(-100)  
        gap = []  

        pre_position = True 
        noreadcount = 0  
        for i in range(len(Coverage)):
            if Coverage[i] == 0:
                if pre_position == True:
                    noreadcount += 1
                else:
                    noreadcount = 1
                    pre_position = True

            else:
                pre_position = False
                if noreadcount != 0:
                    gap.append(noreadcount)
                noreadcount = 0

        #计算gap得分
        gapscore = 0.0
        for k in range(len(gap)):
            #gapscore += math.exp(gap[k])
            if gap[k]!=0:
                #gapscore += 1.0/gap[k]   #gap1
                gapscore += gap[k]     #gap2

        if len(gap) != 0:
            gapscore = gapscore/len(gap)

        file_gapscore.write(cov.strip().split()[0]+' '+str(gapscore)+'\n')

    file_gapscore.close()
    file_cov.close()



def Coverage():
    file_cov = open('Data/Result/RefCoverage.txt', 'r')
    file_out = open('Data/Result/RefCovRate.txt','w')
    for line in file_cov:
        linex = line.strip().split()
        covRate = 0.0
        for i in range(1,len(linex),1):
            if eval(linex[i])!=0:
                covRate += 1
        covRate = covRate/(len(linex)-1)
        file_out.write(linex[0]+' '+str(covRate))
        file_out.write('\n')
    file_out.close()
    file_cov.close()






def HVRScore():
    file_in = open('Data/Result/RefCoverage.txt','r')
    file_out = open('Data/Result/RefHVRScore.txt','w')
    for line in file_in:
        linex = line.strip().split()
        file_out.write(linex[0]+' ')

        cov = [] #记录ref的覆盖度
        for j in range(1,len(linex),1):
            cov.append(eval(linex[j]))


        file_hvr = open('Data/Result/HVR.txt','r')
        hvrlist = []  

        for hvr in file_hvr:
            hvrx = hvr.strip().split()
            if hvrx[0]==linex[0]:
                for i in range(1,len(hvrx),1):
                    hvrxx = hvrx[i].split(',')
                    hvrlist.append(eval(hvrxx[0]))
                    hvrlist.append(eval(hvrxx[1]))
                break
        file_hvr.close()

        hvrscore = 0.0
        for k in range(0,len(hvrlist),2):
            gap = []
           
            Coverage = cov[hvrlist[k]-1 : hvrlist[k+1]]
            Coverage.append(-100)  
            pre_position = True  
            noreadcount = 0  
            for i in range(len(Coverage)):
                if Coverage[i] == 0:
                    if pre_position == True:
                        noreadcount += 1
                    else:
                        noreadcount = 1
                        pre_position = True

                else:
                    pre_position = False
                    if noreadcount != 0:
                        gap.append(noreadcount)
                    noreadcount = 0

            hvrscore += float(sum(cov[hvrlist[k]-1:hvrlist[k+1]]))/sum(cov) * math.exp(-float(sum(gap))/(hvrlist[k+1]-hvrlist[k]+1))
            #hvrscore += math.exp(-float(sum(gap))/(hvrlist[k+1]-hvrlist[k]+1))
        file_out.write(str(hvrscore)+'\n')
    file_out.close()
    file_in.close()





def GetData1():
  
    file = open('Data/PredictData/DataSet.txt', 'w')
    #file_fre = open('Data/Result/RefFreErf.txt','r')
    file_fre = open('Data/Result/RefFre.txt', 'r')
    frex = file_fre.readline()
    file_gamma = open('Data/Result/RefGamma.txt','r')
    gammax = file_gamma.readline()
    file_gapscore = open('Data/Result/GapScore.txt','r')
    gapscorex = file_gapscore.readline()
    file_covrate = open('Data/Result/RefCovRate.txt','r')
    covrate = file_covrate.readline()
    file_hvrscore = open('Data/Result/RefHVRScore.txt','r')
    hvrscore = file_hvrscore.readline()
    file_cluster = open('Data/Result/RefCluster.txt','r')
    cluster = file_cluster.readline()


    while frex!='' and gammax!='' and gapscorex!='' and covrate!='' and hvrscore!='' and cluster!='':
        refno = frex.strip().split()[0]  
        label = ''

        file_fasta = open('Data/Ref_Align/TotalRef.fasta', 'r')
        for fastax in file_fasta:
            if fastax[0]=='>':
                fasta = fastax.strip().split()
                refnofasta = fasta[0][1:len(fasta[0])]
                if refno == refnofasta:
                    label = fasta[len(fastax.strip().split())-1]
                    break
        file_fasta.close()

        if frex.strip().split()[0] != gammax.strip().split()[0]:
            frex = file_fre.readline()


        file.write(refno+' '+frex.strip().split()[1] + ' ' + gammax.strip().split()[1] + ' ' + gammax.strip().split()[2]
                   + ' ' + gapscorex.strip().split()[1] + ' ' +covrate.strip().split()[1]+' '+ hvrscore.strip().split()[1]
                   +' '+cluster.strip().split()[1]+ '\n')

        frex = file_fre.readline()
        gammax = file_gamma.readline()
        gapscorex = file_gapscore.readline()
        covrate = file_covrate.readline()
        hvrscore = file_hvrscore.readline()
        cluster = file_cluster.readline()

    file.close()
    file_fre.close()
    file_gamma.close()
    file_gapscore.close()
    file_covrate.close()
    file_hvrscore.close()
    file_cluster.close()





def GetData():
   
    file_Data = open('Data/PredictData/DataSet.txt', 'w')
    file_covrate = open('Data/Result/RefCovRate.txt','r')
    file_hvrscore = open('Data/Result/RefHVRScore.txt','r')


   
    Ref_No=[]
    dic_gap={}
    dic_hvr={}
    dic_cov={}

    file_gapscore = open('Data/Result/GapScore.txt', 'r')
    for line in file_gapscore:
        Ref_No.append(line.strip().split()[0])
    file_gapscore.close()


    file_gapscore = open('Data/Result/GapScore.txt', 'r')
    for line in file_gapscore:
        linex=line.strip().split()
        dic_gap[linex[0]]=eval(linex[1])

    for line in file_covrate:
        linex=line.strip().split()
        dic_cov[linex[0]]=eval(linex[1])

    for line in file_hvrscore:
        linex=line.strip().split()
        dic_hvr[linex[0]]=eval(linex[1])


   
    for k in Ref_No:
        if k in dic_gap.keys() and k in dic_hvr.keys() and k in dic_cov.keys():
            file_Data.write(str(k)+' '+str(dic_gap[k])+' '+str(dic_cov[k])+' '+str(dic_hvr[k])+' '+'\n')


    file_Data.close()
    file_gapscore.close()
    file_covrate.close()
    file_hvrscore.close()






RefCov()   
RefGap()   
Coverage()  
HVRScore()  
GetData()  
