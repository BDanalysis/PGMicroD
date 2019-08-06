# coding:utf-8
import pysam

'''
功能： 计算物种浓度
参数： read总量
返回： 
'''


Dict={}   #存储预测结果，格式: 物种名称:浓度


def RefFre():
    N = 0
    # 统计read总量
    file = open('Data/example.reads.fq', 'r')
    for line in file:
        N = N + 1
    file.close()
    N = N / 2

    # 记录浓度
    file_out = open('Data/Result/RefFre.txt', 'w')
    file_refno = open('Data/Ref_Align/exist_Species.fasta', 'r')

    ref = ''  # 记录物种序列

    line = file_refno.readline()
    while line != '':
        if line[0] == '>':
            i = line.strip().split()[0][1:]  # 物种号

            name = (line.strip().split()[1]).split(';')[-1]  # 物种

            line = file_refno.readline()
            ref = line.strip()

            sum = 0  # 物种下read的总和

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



RefFre()  # 计算ref的频率
