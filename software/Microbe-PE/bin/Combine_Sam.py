# coding:utf-8


def Comebine_Sam():
    file_out=open('Data/Ref_Align/pre_Density.sam','w')

    file_in1=open('Data/Ref_Align/example.reads.sam','r') 
    for line1 in file_in1:
        file_out.write(line1)
    file_in1.close()

    file_in2=open('Data/Ref_Align/Recovery.sam','r')
    for line2 in file_in2:
        if line2[0]!='@':
            file_out.write(line2)
    file_in2.close()

    file_out.close()

Comebine_Sam()
