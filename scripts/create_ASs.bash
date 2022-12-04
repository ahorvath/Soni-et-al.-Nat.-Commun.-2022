module load bedtools 
module load R/4.0.0

SAMPLE=sp_st3237_RED1KO_polyA_RNAseq
for FC in 1.2 1.5 2; do
	echo ${FC}
	cat ${SAMPLE}.?_FC${FC}_min100.bed > ${SAMPLE}_FC${FC}_min100.bed

	awk -F"\t" '{OFS="\t"; print $1,($6=="+" ? $2:($3-1)),($6=="+"?($2+1):$3),$4,$5,$6}' ${SAMPLE}_FC${FC}_min100.bed > ${SAMPLE}_FC${FC}_min100_starpos.bed
	awk -F"\t" '{OFS="\t"; print $1,($6=="-" ? $2:($3-1)),($6=="+"?($2+1):$3),$4,$5,$6}' ${SAMPLE}_FC${FC}_min100.bed > ${SAMPLE}_FC${FC}_min100_endpos.bed

	intersectBed -S -a CUTs_minus_PROMPTs_FC${FC}.bed -b ../ONT_final/pombe_mRNAs_TSS.bed  > ASs_FC${FC}.bed

	for bam in `ls /g/data1a/xc17/Fischer/Nanopore/renamed_bam/ | grep Red | grep .bam$`; do
		name=`basename ${bam} .bam`
		echo ${name}
		intersectBed -wa -wb -bed -b /g/data1a/xc17/Fischer/Nanopore/renamed_bam/${bam} -a ASs_FC${FC}.bed > ${name}_ASs_FC${FC}.bed 
	done

#	Rscript infer_poly_tails_based_on_regions_group_by_regions.R CUTs_minus_PROMPTs_FC${FC}.bed > `basename ${bed} .bed`_numbers.log
#	Rscript infer_poly_tails_based_on_regions_group_by_regions.R CUTs_min100_FC${FC}.bed > `basename ${bed} .bed`_numbers.log
#	Rscript infer_poly_tails_based_on_regions_group_by_regions.R PROMPTs_closest_pm250togene_max1500_FC${FC}.bed > `basename ${bed} .bed`_numbers.log
	Rscript infer_poly_tails_based_on_regions_group_by_regions_v2.R ASs_FC${FC}.bed > `basename ${bed} .bed`_numbers.log
done
