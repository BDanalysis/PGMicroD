date +'start time=%Y-%m-%d %k:%M:%S.%N'
echo "================================================================================"


if [ $# -ne 3 ];then
echo "command formate: sh MicrobeRun.sh FastaDir ReadDir ResultDir"
exit
fi


mkdir ./input
mkdir ./output
mkdir ./bin/Data
mkdir ./bin/Data/HVRFile
mkdir ./bin/Data/PredictData
mkdir ./bin/Data/Ref_Align
mkdir ./bin/Data/Result



cp $1 ./input/TotalRef.fasta
cp $2 ./input/example.reads.fq


# move fasta and fastq file into workshop space
cp ./input/TotalRef.fasta ./bin/Data/Ref_Align/TotalRef1.fasta
cp ./input/example.reads.fq ./bin/Data/example.reads.fq



fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'AGYGGCGNACGGGTGAGTAA' -outfile ./bin/Data/HVRFile/V2.fuzznuc
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'CCTACGGGAGGCAGCAG' -outfile ./bin/Data/HVRFile/V3.fuzznuc
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'AYTGGGYDTAAAGNG' -outfile ./bin/Data/HVRFile/V4.fuzznuc
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'AGGATTAGATACCCT' -outfile ./bin/Data/HVRFile/V5.fuzznuc
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'TCGAtGCAACGCGAAGAA' -outfile ./bin/Data/HVRFile/V6.fuzznuc
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'GYAACGAGCGCAACCC' -outfile ./bin/Data/HVRFile/V7.fuzznuc
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'ATGGCTGTCGTCAGCT' -outfile ./bin/Data/HVRFile/V8.fuzznuc




fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'TGCTGCCTCCCGTAGGAGT' -outfile ./bin/Data/HVRFile/V2_1.fuzznuc
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'ATTACCGCGGCTGCTGG' -outfile ./bin/Data/HVRFile/V3_1.fuzznuc
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'TACNVGGGTATCTAATCC' -outfile ./bin/Data/HVRFile/V4_1.fuzznuc
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'CCGTCAATTCCTTTGAGTTT' -outfile ./bin/Data/HVRFile/V5_1.fuzznuc
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'ACATTTCACAACACGAGCTGACGA' -outfile ./bin/Data/HVRFile/V6_1.fuzznuc
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'GTAGCRCGTGTGTMGCCC' -outfile ./bin/Data/HVRFile/V7_1.fuzznuc
fuzznuc -sequence ./bin/Data/Ref_Align/TotalRef1.fasta -pattern 'ACGGGCGGTGTGTAC' -outfile ./bin/Data/HVRFile/V8_1.fuzznuc





cd bin
sh ProgramRun.sh


cd ../
cp ./bin/Data/Result/RefFre.txt ./output/RefResult.txt
cp ./output/RefResult.txt $3


rm -rf ./input
rm -rf ./output
rm -rf ./bin/Data/HVRFile
rm -rf ./bin/Data/PredictData
rm -rf ./bin/Data/Ref_Align
rm -rf ./bin/Data/Result
rm -rf ./bin/Data




echo "================================================================================"
date +'end time=%Y-%m-%d %k:%M:%S.%N'
