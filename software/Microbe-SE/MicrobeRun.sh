date +'start time=%Y-%m-%d %k:%M:%S.%N'
echo "================================================================================"


if [ $# -ne 2 ];then
echo "command formate: sh MicrobeRun.sh ReadDir ResultDir"
exit
fi



cp $1 ./input/example.reads.fq


#将input的fasta和fastq移动到工作目录
cp ./input/TotalRef.fasta ./bin/Data/Ref_Align/TotalRef1.fasta
cp ./input/example.reads.fq ./bin/Data/example.reads.fq


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
cp ./output/RefResult.txt $2


#删除中间结果文件
rm -rf ./input/example.reads.fq
rm -rf ./output/RefResult.txt

rm -rf ./bin/Data/HVRFile/*
rm -rf ./bin/Data/PredictData/*
rm -rf ./bin/Data/Ref_Align/*
rm -rf ./bin/Data/Result/*
rm -rf ./bin/Data/example.reads.fq






echo "================================================================================"
date +'end time=%Y-%m-%d %k:%M:%S.%N'