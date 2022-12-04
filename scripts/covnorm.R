library(rtracklayer)

tt <- read.table("samples_YNK3.lst")

for (sample.id in tt[,1]) {
#        sample.id <- tt[1,2]
        print(sample.id)

        bw.f <- list.files("../renamed_bw", pattern = paste0(sample.id, ".f.bw"), full.names = TRUE)
        bw.r <- list.files("../renamed_bw", pattern = paste0(sample.id, ".r.bw"), full.names = TRUE)

        bw.total.f <- import.bw(BigWigFile(bw.f),  as = 'RleList')
        bw.total.r <- import.bw(BigWigFile(bw.r),  as = 'RleList')

	sum.cov.f <- sum(as.vector(bw.total.f$I)) + sum(as.vector(bw.total.f$II))
	sum.cov.r <- sum(as.vector(bw.total.r$I)) + sum(as.vector(bw.total.r$II))

        norm.bw.total.f <- bw.total.f
        norm.bw.total.r <- bw.total.r

	for(seq.name in names(norm.bw.total.f)) {
		print(seq.name)
		norm.bw.total.f[[seq.name]]@values <- norm.bw.total.f[[seq.name]]@values/(sum.cov.f/1e8)
		norm.bw.total.r[[seq.name]]@values <- norm.bw.total.r[[seq.name]]@values/(sum.cov.r/1e8)
	}

        export.bw(norm.bw.total.f, paste0("covnorm/", sample.id, ".covnorm.f.bw"))
        export.bw(norm.bw.total.r, paste0("covnorm/", sample.id, ".covnorm.r.bw"))
}

