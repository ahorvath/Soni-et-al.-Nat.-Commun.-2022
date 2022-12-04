### head for the job submission
#!/bin/bash
#PBS -q gpuvolta
#PBS -l ngpus=4 
#PBS -l ncpus=48
#PBS -P xc17
#PBS -l walltime=8:00:00
#PBS -l mem=128GB
#PBS -l jobfs=100GB
#PBS -l storage=gdata/xc17

module load samtools/1.10
module load bedtools/2.28.0
### map to the pombe genome
export PATH="/tools/minimap2/:${PATH}"
chromosome_fasta="/index/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa"
project_dir="Nanopore/"
export PATH="/bin/:${PATH}"
mkdir -p ${project_dir}/qc
for lib in <sample-names>; do 

	echo ${lib}

	data_folder="raw_data/${lib}"
	output_folder="${project_dir}/${lib}/output"
	mkdir -p ${output_folder}
	analysis_directory="${project_dir}/${lib}/k5"
	mkdir -p ${analysis_directory}
#	zcat ${output_folder}/*.fastq.gz > ${output_folder}/all.fastq
	minimap2 -ax splice -k7 ${chromosome_fasta} ${output_folder}/all.fastq > ${analysis_directory}/${lib}_alignment.sam 2> ${analysis_directory}/${lib}_alignment.log

	samtools view -@ 40 -bS ${analysis_directory}/${lib}_alignment.sam |\
	samtools sort -@ 40 > ${analysis_directory}/${lib}_alignment.bam
	samtools index ${analysis_directory}/${lib}_alignment.bam
	cat ${output_folder}/all.fastq | grep -c "read" > ${project_dir}/qc/${lib}_k5_all_reads.tsv

	bedtools bamtobed -i ${analysis_directory}/${lib}_alignment.bam | awk '!seen[$4]++' | wc -l > ${project_dir}/qc/${lib}_k5_uniq_mapped.tsv
done

