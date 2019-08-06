#coding:utf-8
import pysam

import math
import os

'''
功能： 计算Ref的可变区范围
参数： 可变区位置i
返回： 所有ref可变区i的字典
'''


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


'''
功能： 初始化计算read数量N，ref名称列表RefNo，ref的可变区HVR
参数： N，RefNo，HVR
返回值： N，RefNo，HVR
'''


def init():
    N = 0
    RefNo = []
    HVR = {}  # ref的可变区

    file = open('Data/example.reads1.fq', 'r')
    for line in file:
        N = N + 1
    file.close()
    N = N / 2  #双端测序的read量

    # 1、生成RefNo.txt
    file = open('Data/Ref_Align/TotalRef.fasta', 'r')
    for line in file:
        if line[0] == '>':
            num = line.strip().split()[0]
            RefNo.append(num[1:len(num)])
    file.close()

    # 2、由RefNo生成HVR的空值的键字典
    for i in range(len(RefNo)):
        HVR[RefNo[i]] = []

    # 3、循环所有可变区2~9，将HVRx合并加入HVR
    for i in range(2, 9, 1):
        HVRx = singleHVR(i)
        # 将HVRx合并到HVR
        for (k, v) in HVRx.items():
            if k in HVR.keys():
                for (x0, x1) in HVRx[k]:
                    # 当后续数字大于前者时加入
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


'''
功能：将CIGAR转为数值和字符组合
参数：CIGAR字符串
返回：数字字符列表
'''


def CigerTrans(CIGAR):
    l = []

    num = ''  # 存放数字
    flag = True  # 上次是数字
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


'''
功能：转化read与ref
参数：原read、原ref、比对信息CIGAR、质量分数score、比对起始位置start
返回值：转换后的：read_tran, ref_tran, score, 比对信息列表l
'''


def ReadRefTrans(read, ref, CIGAR, score, start):
    l = CigerTrans(CIGAR)  # 比对信息
    # 将read与ref的碱基碱基序列转化
    ref_tran = ''  # 转化后的序列
    read_tran = ''
    score_tran = ''  # 将D在read的缺失用空格标识

    ref_start = start  # 当前起始位置,使用逻辑值
    read_start = 1

    for i in range(0, len(l), 2):
        num = l[i]  # 数值
        info = l[i + 1]  # 匹配信息

        if info == 'M':

            read_tran += read[read_start - 1: read_start + num - 1]
            score_tran += score[read_start - 1: read_start + num - 1]
            read_start = read_start + num  # 新起始位置

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
            # ref_start不变

        elif info == 'D':
            read_tran += ' '
            score_tran += ' '
            # read_start不变

            ref_tran += ref[ref_start - 1: ref_start + num - 1]
            ref_start = ref_start + num

        elif info == 'H':
            # read已截掉H部分
            continue
        else:
            continue

    return [ref_tran, read_tran, score_tran, l]


'''
 功能： 求出read落在ref的HVR的长度
 参数： Ref的HVR位点列表position, read比对到ref的起始位点start, 结束位点end
 返回： read落在Ref的HVR的片段长度比例 （read落在可变区长度/可变区总长度）  hvrLen
'''


def HVRlen(position, start, end):
    sumHVR = 1  # 防止出现除0错误
    for (x, y) in position:
        sumHVR += y - x + 1

    hvrLen = 0
    space = []  # 使用position位点分割read比对位置
    space.append(start)

    for (x0, x1) in position:
        if x0 >= start and x0 <= end:
            space.append(x0)

        if x1 >= start and x1 <= end:
            space.append(x1)

    space.append(end)

    if len(space) > 2:  # 说明有read落在HVR区域

        # 判断start是否落在HVR区域
        flag = False  # 默认认为start不落在HVR区域
        for (x0, x1) in position:
            if x0 <= start and start >= x1:
                flag = True
                break

        if flag == True:  # start落在变异区
            for i in range(0, len(space), 2):
                if i + 1 < len(space):
                    hvrLen += space[i + 1] - space[i] + 1
        else:
            for i in range(1, len(space), 2):
                if i + 1 < len(space):
                    hvrLen += space[i + 1] - space[i] + 1

    return float(hvrLen + 1) / sumHVR


'''
计算sam下read到ref的归属度
'''


def SamLikelihood():
    # 删除example.reads.fq.sam中未比对上的read, 生成example.sam
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

        # 1、取出ref的HVR列表
        position = []  # HVR列表
        file_HVR = open('Data/Result/HVR.txt', 'r')
        for refnox in file_HVR:
            refno = refnox.strip().split()[0]

            if refno == samx[2]:
                for HVR in refnox.strip().split()[1:]:
                    [x0, x1] = HVR.split(',')
                    position.append([eval(x0), eval(x1)])
                break
        file_HVR.close()

        # 2、取出ref比对序列
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

        # 3、将read与ref的碱基碱基序列转化
        score_phred = samx[10]  # 第10位是比对质量分数
        list = ReadRefTrans(samx[9], ref, samx[5], score_phred, eval(samx[3]))  # 第10位是比对质量分数 第3位是ref比对的起始位置

        # 4、计算read在ref的起始与终止位置，计算HVR权重w
        w = 0.0  # HVR权重值
        align_start = eval(samx[3])  # 比对到ref的起始位置
        info = CigerTrans(samx[5])  # 比对信息转化
        align_len = 0
        for j in range(0, len(info), 2):
            align_len += info[j]
        align_end = align_start + align_len - 1  # 比对到ref的终止位置

        w = HVRlen(position, align_start, align_end)  # 求出read落在ref的HVR的长度

        # 5、计算似然, 使得转换后的ref、read、score每个位置均不为空
        # list[0]是ref_tran  list[1]是read_tran  list[2]是score_tran  list[3]是比对信息
        ref_tran = list[0]
        read_tran = list[1]
        score_tran = list[2]
        lsigar = list[3]
        multi = 1  # 错误率乘积

        # 取出三者的最小值
        l = [len(ref_tran), len(read_tran), len(score_tran)]

        for j in range(min(l)):  # x是每个碱基的质量分数
            if ref_tran[j] != ' ' and read_tran[j] != ' ' and score_tran[j] != ' ':
                q = (math.pow(10, -((ord(score_tran[j]) - 33)) / 10))  # 错误率

                if ref_tran[j] == read_tran[j]:
                    multi = multi * (1.0 - q)

                else:
                    multi = multi * q / 3.0

        # 融合gap
        gapsum = 1.0  # 防止分母为0
        for j in range(0, len(lsigar), 2):
            if lsigar[j + 1] in ['S', 'I', 'D', 'H']:
                # gapsum +=  math.log(math.exp(lsigar[j]+1))/math.log(2)
                gapsum += math.log(lsigar[j] + 1) / math.log(2)

        if gapsum == 0.0:
            like = multi + w
        else:
            like = (multi + w) / gapsum  # like = (multi +  w) / gapsum

        file_readlike.write(samx[0] + ' ' + str(like) + ' ')
        for i in range(1,len(samx),1):
            file_readlike.write(samx[i]+' ')
        file_readlike.write('\n')

    file_readlike.close()
    file_sam.close()


'''
功能：从samLike.txt统计计算mu，sigma, 统计出低似然的read号，并存入lowRead.txt
参数：过滤threold
返回：低似然read号列表
'''
def LowReadLikelihoodNo(T):

    Recovery_Read = set()#待回收read
    lowLikeRead = set() #low likelihood read

    file_in = open('Data/Result/samLike.txt', 'r')
    for like in file_in:
        likex = like.strip().split()
        if eval(likex[1]) <= T:
            lowLikeRead.add(likex[0])
            #M-Map, No D I S H
            #回收全比对read
            C={'D','I','S','H'}
            if len(set(likex[2]) & C)==0:
                Recovery_Read.add(likex[0])
            #回收对比对read
            for i in range(12, len(likex), 1):
                if 'XA:Z:' in likex[i]:
                    Recovery_Read.add(likex[0])

    file_in.close()


    #去除example.sam中的低归属度read
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




    #将待回收read存入Recovery.fastq 格式：read号 序列 质量值
    file_lowRead = open('Data/Result/Recovery.fastq', 'w')
    for x in Recovery_Read:
        fq1_flag = False

        #判断是否在fq1中
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

        #判断是否在fq2中
        if fq1_flag==False:  #不在fq1中
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




init()  # read总数量
SamLikelihood()  # 计算归属度
T = 0.44
LowReadLikelihoodNo(T)  # 过滤低归属read，得example.reads.sam
