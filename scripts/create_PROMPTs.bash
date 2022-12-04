module load bedtools 
module load R/4.0.0

Rscript calc_ratios_covnorm_chrI_II.R

SAMPLE=sp_st3237_RED1KO_polyA_RNAseq
#for FC in 1.2 1.5 2; do
	echo ${FC}
	cat ${SAMPLE}.?_FC${FC}_min100.bed > ${SAMPLE}_FC${FC}_min100.bed

	awk -F"\t" '{OFS="\t"; print $1,($6=="+" ? $2:($3-1)),($6=="+"?($2+1):$3),$4,$5,$6}' ${SAMPLE}_FC${FC}_min100.bed > ${SAMPLE}_FC${FC}_min100_starpos.bed

	closestBed -D ref -S -a <(sortBed -i ${SAMPLE}_FC${FC}_min100_starpos.bed) -b <(sortBed -i ../ONT_final/pombe_mRNAs_TSS.bed) | awk -F"\t" '$13 >= -250 && $13 <= 250 && $5 <= 1500' | awk -F"\t" '{OFS="\t"; print $4,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | sed 's/:/\t/;s/-/\t/' | cut -f1-6 > PROMPTs_closest_pm250togene_max1500_FC${FC}.bed

	ln -s ${SAMPLE}_FC${FC}_min100.bed CUTs_min100_FC${FC}.bed
	subtractBed -s -A -a ${SAMPLE}_FC${FC}_min100.bed -b PROMPTs_closest_pm250togene_max1500_FC${FC}.bed > CUTs_minus_PROMPTs_FC${FC}.bed

	for bam in `ls /g/data1a/xc17/Fischer/Nanopore/renamed_bam/ | grep Red | grep .bam$`; do
		name=`basename ${bam} .bam`
		echo ${name}
		intersectBed -wa -wb -bed -b /g/data1a/xc17/Fischer/Nanopore/renamed_bam/${bam} -a PROMPTs_closest_pm250togene_max1500_FC${FC}.bed > ${name}_PROMPTs_closest_pm250togene_max1500_FC${FC}.bed 
		intersectBed -wa -wb -bed -b /g/data1a/xc17/Fischer/Nanopore/renamed_bam/${bam} -a CUTs_min100_FC${FC}.bed > ${name}_CUTs_min100_FC${FC}.bed 
		intersectBed -wa -wb -bed -b /g/data1a/xc17/Fischer/Nanopore/renamed_bam/${bam} -a CUTs_minus_PROMPTs_FC${FC}.bed > ${name}_CUTs_minus_PROMPTs_FC${FC}.bed 
	done

	Rscript infer_poly_tails_based_on_regions_group_by_regions_v2.R CUTs_minus_PROMPTs_FC${FC}.bed > `basename ${bed} .bed`_numbers.log
	Rscript infer_poly_tails_based_on_regions_group_by_regions_v2.R CUTs_min100_FC${FC}.bed > `basename ${bed} .bed`_numbers.log
	Rscript infer_poly_tails_based_on_regions_group_by_regions_v2.R PROMPTs_closest_pm250togene_max1500_FC${FC}.bed > `basename ${bed} .bed`_numbers.log
done
