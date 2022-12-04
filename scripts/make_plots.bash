

for bed in CUTs_min100_FC1.2.bed  CUTs_min100_FC1.5.bed  CUTs_min100_FC2.bed  PROMPTs_closest_pm250togene_max1500_FC1.2.bed  PROMPTs_closest_pm250togene_max1500_FC1.5.bed  PROMPTs_closest_pm250togene_max1500_FC2.bed; do
	echo ${bed}

	Rscript infer_poly_tails_based_on_regions_group_by_regions_v2.R > `basename ${bed} .bed`_numbers.bed
done
