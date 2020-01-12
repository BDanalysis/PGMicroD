bwa index Data/Ref_Align/TotalRef1.fasta
bwa mem Data/Ref_Align/TotalRef1.fasta Data/example.reads1.fq Data/example.reads2.fq > Data/Ref_Align/example.reads.fq.sam


python3 ExtractRef.py

python3 ReadLikelihood.py

samtools view -Sb Data/Ref_Align/example.reads.sam > Data/Ref_Align/example.reads1.bam
samtools sort -O bam -o Data/Ref_Align/example.reads.bam Data/Ref_Align/example.reads1.bam
samtools index Data/Ref_Align/example.reads.bam

python3 RefCharacter.py
python3 ModelPredict.py


cp Data/Result/exist_Species.fasta Data/Ref_Align/exist_Species.fasta
cp Data/Result/Recovery.fastq Data/Ref_Align/Recovery.fastq

bwa index Data/Ref_Align/exist_Species.fasta
bwa mem Data/Ref_Align/exist_Species.fasta Data/Ref_Align/Recovery.fastq > Data/Ref_Align/Recovery.sam

python3 Combine_Sam.py

samtools view -Sb Data/Ref_Align/pre_Density.sam > Data/Ref_Align/pre_Density1.bam
samtools sort -O bam -o Data/Ref_Align/pre_Density.bam Data/Ref_Align/pre_Density1.bam
samtools index Data/Ref_Align/pre_Density.bam

python3 Density.py
