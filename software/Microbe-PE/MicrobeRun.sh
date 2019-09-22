date +'start time=%Y-%m-%d %k:%M:%S.%N'
echo "================================================================================"


if [ $# -ne 4 ];then
echo "command formate: sh MicrobeRun.sh FastaDir Read1Dir Read2Dir ResultDir"
exit
fi

cp $1 ./input/TotalRef.fasta
cp $2 ./input/example.reads1.fq
cp $3 ./input/example.reads2.fq




#将input的fasta和fastq移动到工作目录
cp ./input/TotalRef.fasta ./bin/Data/Ref_Align/TotalRef1.fasta
cp ./input/example.reads1.fq ./bin/Data/example.reads1.fq
cp ./input/example.reads2.fq ./bin/Data/example.reads2.fq



#获取HVR
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'AGYGGCGNACGGGTGAGTAA' -outfile ./bin/Data/HVRFile/V2.fuzznuc
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'CCTACGGGAGGCAGCAG' -outfile ./bin/Data/HVRFile/V3.fuzznuc
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'AYTGGGYDTAAAGNG' -outfile ./bin/Data/HVRFile/V4.fuzznuc
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'AGGATTAGATACCCT' -outfile ./bin/Data/HVRFile/V5.fuzznuc
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'TCGAtGCAACGCGAAGAA' -outfile ./bin/Data/HVRFile/V6.fuzznuc
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'GYAACGAGCGCAACCC' -outfile ./bin/Data/HVRFile/V7.fuzznuc
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'ATGGCTGTCGTCAGCT' -outfile ./bin/Data/HVRFile/V8.fuzznuc



#运行程序
cd bin
sh ProgramRun.sh

#将结果移动到output
cd ../
cp ./bin/Data/Result/RefFre.txt ./output/RefResult.txt
cp ./output/RefResult.txt $4




#删除中间结果文件
rm -rf ./input/TotalRef.fasta
rm -rf ./input/example.reads1.fq
rm -rf ./input/example.reads2.fq
rm -rf ./output/RefResult.txt

rm -rf ./bin/Data/HVRFile/*
rm -rf ./bin/Data/PredictData/*
rm -rf ./bin/Data/Ref_Align/*
rm -rf ./bin/Data/Result/*
rm -rf ./bin/Data/example.reads1.fq
rm -rf ./bin/Data/example.reads2.fq




echo "================================================================================"
date +'end time=%Y-%m-%d %k:%M:%S.%N'
