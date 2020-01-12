# coding:utf-8
from sklearn import svm
from sklearn.preprocessing import StandardScaler
from sklearn.externals import joblib  


def SVMPredict():
    dirs = 'Model'
    Grid = joblib.load(dirs + '/svm.pkl')


    
 
    X_test = []  
    y_refNo = [] 
    file_in = open('Data/PredictData/DataSet.txt', 'r')
    for line in file_in:
        x = []
        linex = line.strip().split()
        y_refNo.append(linex[0])

        for i in range(1, 4, 1):
            x.append(eval(linex[i]))
        X_test.append(x)
    file_in.close()

    y_test = Grid.predict(X_test)  


    

    file_in=open('Data/Result/exist_Species_No.txt','w')
    for i in range(0, len(y_test), 1):
        if y_test[i] == 1:
            file_in.write(y_refNo[i]+'\n')
    file_in.close()
    


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



SVMPredict() 
Exist_Species_Fasta()     
