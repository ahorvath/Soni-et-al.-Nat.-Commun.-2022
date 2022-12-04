library(tailfindr)
ff <- c("<path-to-fast5>")

setwd("<workdir>")

#for (i in 1:length(ff)) {
	filepath <- ff[i]
	print(filepath)
	sample.id <- unlist(as.data.frame(strsplit(ff, split = "/"))[7,i])
	print(sample.id)
	df <- find_tails(fast5_dir = filepath,
                 save_dir = './tailfindr',
                 csv_filename = paste0(sample.id, "_rna_tails.csv"),
                 num_cores = 96,
		 basecall_group = 'Basecall_1D_000',
                 save_plots = F,
                 plot_debug_traces = TRUE,
                 plotting_library = 'rbokeh')

}
