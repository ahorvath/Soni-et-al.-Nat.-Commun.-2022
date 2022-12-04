PREFIX_INPUT=.
PREFIX_OUTPUT=.
PICARD=/home/admin/Programs/picard.jar
JAVA=java
JAR="${JAVA} -jar -Xmx16g"
TRIMMO=/home/admin/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar
GENOME=/data/index/pombase/ASM294v2.fa
INDEX=/data/index/pombase/ASM294v2

#conda activate khmer
#unique-kmers.py ${GENOME} -k 75 #=> 11693438
SIZE=12291456
SAMPLES=(<sample1> <sample2> <sample3>)
mkdir -p {bigwigs,hisat2}
for I in 0 1 2 3 ; do
       echo ${INDEX}
       SAMPLE=${SAMPLES[${I}]}
       echo ${SAMPLE}

       hisat2 --rg-id ${SAMPLE}  --no-spliced-alignment  --phred33 -p 30 --new-summary -x ${INDEX} -1 FQ/${SAMPLE}_1.fastq.gz -2 FQ/${SAMPLE}_2.fastq.gz > hisat2/${SAMPLE}_unsorted.sam 2> hisat2/${SAMPLE}.stat
       samtools sort -@ 30 hisat2/${SAMPLE}_unsorted.sam > hisat2/${SAMPLE}.bam
       rm -rfv hisat2/${SAMPLE}_unsorted.sam
      #cp /media/admin/Data/hisat2/${SAMPLE}.bam hisat2
       samtools index -@ 30 hisat2/${SAMPLE}.bam
	
       for NORM in RPKM; do
               bamCoverage -p 40 --filterRNAstrand forward -b hisat2/${SAMPLE}.bam -o bigwigs/${SAMPLE}.r_${NORM}.bw --binSize 1 --normalizeUsing ${NORM} --effectiveGenomeSize ${SIZE} 
               bamCoverage -p 40 --filterRNAstrand reverse -b hisat2/${SAMPLE}.bam -o bigwigs/${SAMPLE}.f_${NORM}.bw --binSize 1 --normalizeUsing ${NORM} --effectiveGenomeSize ${SIZE}  
       done
done

