# coding:utf-8
from sklearn import svm
from sklearn.preprocessing import StandardScaler

'''
@功能：通过训练数据得出SVM模型，读入测试数据测试
@日期：2018-12-18
'''

'''
训练SVM模型
'''


def SVMPredict():
    # 读数据
    X_train = []
    Y_train = []
    file_in = open('Data/TrainSet/TrainDataSet.txt', 'r')
    for line in file_in:
        x = []
        linex = line.strip().split()
        Y_train.append(eval(linex[0]))

        for i in range(4, 7, 1):
            x.append(eval(linex[i]))
        X_train.append(x)
    file_in.close()

    # 预测
    X_test = []  # 预测数据
    y_test = []  # 存储预测的结果
    y_refNo = []  # ref编号
    file_in = open('Data/PredictData/DataSet.txt', 'r')
    for line in file_in:
        x = []
        linex = line.strip().split()
        y_refNo.append(linex[0])

        for i in range(1, 4, 1):
            x.append(eval(linex[i]))
        X_test.append(x)
    file_in.close()

    # 训练分类器
    sc = StandardScaler()
    sc.fit(X_train)
    X_train = sc.transform(X_train)
    X_test = sc.transform(X_test)
    clf = svm.SVC(kernel='rbf', C=1.0, random_state=0, gamma=0.2)
    clf.fit(X_train, Y_train)

    y_test = clf.predict(X_test)

    #将存在的物种写入Data/Result/exist_Species_No.txt
    file_in=open('Data/Result/exist_Species_No.txt','w')
    for i in range(0, len(y_test), 1):
        if y_test[i] != 0.0:
            file_in.write(y_refNo[i]+'\n')
    file_in.close()


#取出预测出的物种序列
def Exist_Species_Fasta():
    file_out=open('Data/Result/exist_Species.fasta','w')
    file_in=open('Data/Result/exist_Species_No.txt','r')
    for No in file_in:

        file_ref=open('Data/Ref_Align/TotalRef1.fasta','r')
        line=file_ref.readline()
        while line!='':
            if line[0]=='>':
                linex=line.strip().split()
                if '>'+No.strip() == linex[0]:
                    file_out.write(line)
                    line=file_ref.readline()
                    file_out.write(line)
                    break
            line=file_ref.readline()

        file_ref.close()


    file_in.close()
    file_out.close()



SVMPredict()  #预测存在的物种编号
Exist_Species_Fasta()      #取出预测出的物种序列

