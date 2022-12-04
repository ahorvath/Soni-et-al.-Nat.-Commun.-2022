library(rtracklayer)

bw.lst <- list()

comp <- read.table("compare5.tsv")

for (sample.id in unique(c(comp[,1], comp[,2]))) {
        print(sample.id)

        bw.f <- paste0("covnorm/", sample.id, ".covnorm.f.bw")
        bw.r <- paste0("covnorm/", sample.id, ".covnorm.r.bw")

        bw.lst[paste0(sample.id,".f")] <- import(BigWigFile(bw.f), as = 'RleList')
        bw.lst[paste0(sample.id,".r")] <- import(BigWigFile(bw.r), as = 'RleList')
}



for (i in 1:nrow(comp)) {
	print(comp[i,])

	for (direction in c("f", "r")) {
		print(direction)
		IP.sample <- paste0(comp[i,1], ".", direction)
		noIP.sample <- paste0(comp[i,2], ".", direction)

		IP.bw <- bw.lst[[IP.sample]]
		noIP.bw <- bw.lst[[noIP.sample]]

		ratio.list <- list()
		for (seq.name in names(IP.bw)) {
	                values <- log2(as.vector(IP.bw[[seq.name]])+1) - log2(as.vector(noIP.bw[[seq.name]])+1)
			ratio.list[[seq.name]] <- Rle(values)
		}
		ratio.bw <- RleList(ratio.list)
		export.bw(ratio.bw, paste0("covnorm/", IP.sample, ".covnorm.ratio.bw"))
	}
}

