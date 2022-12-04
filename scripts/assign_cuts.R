library(rtracklayer)
library(GenomicRanges)
sample.id <- "sp_st3237_RED1KO_polyA_RNAseq"

bw.name.f <- paste0("covnorm/", sample.id, ".f.covnorm.ratio.bw")
bw.name.r <- paste0("covnorm/", sample.id, ".r.covnorm.ratio.bw")

bw.f <- import.bw(BigWigFile(bw.name.f),  as = 'RleList')
bw.r <- import.bw(BigWigFile(bw.name.r),  as = 'RleList')
FC <- 1.2
filtered.bw.f <- bw.f
filtered.bw.r <- bw.r
flank.size <- 25
min.length <- 100
#max.length <- 1000
for(seq.name in grep("^ERCC", names(bw.f), invert = T, value = T)) {
	print(seq.name)
	filtered.bw.f[[seq.name]]@values <- ifelse(bw.f[[seq.name]]@values < log2(FC), 0, bw.f[[seq.name]]@values)
	filtered.bw.r[[seq.name]]@values <- ifelse(bw.r[[seq.name]]@values < log2(FC), 0, bw.r[[seq.name]]@values)
}

export.bw(filtered.bw.f, paste0("covnorm/", sample.id, ".f.covnorm.ratio.FC", FC, ".bw"))
filtered.gr.f <- as(filtered.bw.f, "GRanges")
filtered.gr.f <- filtered.gr.f[filtered.gr.f$score > 0]
filtered.gr.f <- filtered.gr.f[grep("^ERCC", seqnames(filtered.gr.f), invert = T)]
bed.gr.f <- filtered.gr.f
start(bed.gr.f) <- start(bed.gr.f) - flank.size
end(bed.gr.f) <- end(bed.gr.f) + flank.size
bed.gr.f <- reduce(trim(bed.gr.f))
strand(bed.gr.f) <- "+"
start(bed.gr.f) <- start(bed.gr.f) + flank.size
end(bed.gr.f) <- end(bed.gr.f) - flank.size

bed.gr.f.name <- paste0(seqnames(bed.gr.f), ":", start(bed.gr.f), "-", end(bed.gr.f))
out.f.bed <- cbind(as.data.frame(bed.gr.f)[1:3], bed.gr.f.name, as.data.frame(bed.gr.f)[4:5])

write.table(out.f.bed, file = paste0(sample.id, ".f_FC", FC, ".bed"), quote = F, sep = "\t", row.names = F, col.names = F)
len.out.f.bed <- subset(out.f.bed, width >= min.length)

write.table(len.out.f.bed, file = paste0(sample.id, ".f_FC", FC, "_min", min.length, ".bed"), quote = F, sep = "\t", row.names = F, col.names = F)

export.bw(filtered.bw.r, paste0("covnorm/", sample.id, ".r.covnorm.ratio.FC", FC, ".bw"))
filtered.gr.r <- as(filtered.bw.r, "GRanges")
filtered.gr.r <- filtered.gr.r[filtered.gr.r$score > 0]
filtered.gr.r <- filtered.gr.r[grep("^ERCC", seqnames(filtered.gr.r), invert = T)]
bed.gr.r <- filtered.gr.r
start(bed.gr.r) <- start(bed.gr.r) - flank.size
end(bed.gr.r) <- end(bed.gr.r) + flank.size
bed.gr.r <- reduce(trim(bed.gr.r))
strand(bed.gr.r) <- "-"
start(bed.gr.r) <- start(bed.gr.r) + flank.size
end(bed.gr.r) <- end(bed.gr.r) - flank.size
bed.gr.r.name <- paste0(seqnames(bed.gr.r), ":", start(bed.gr.r), "-", end(bed.gr.r))
out.r.bed <- cbind(as.data.frame(bed.gr.r)[1:3], bed.gr.r.name, as.data.frame(bed.gr.r)[4:5])
write.table(out.r.bed, file = paste0(sample.id, ".r_FC", FC, ".bed"), quote = F, sep = "\t", row.names = F, col.names = F)
len.out.r.bed <- subset(out.r.bed, width >= min.length)
write.table(len.out.r.bed, file = paste0(sample.id, ".r_FC", FC, "_min", min.length, ".bed"), quote = F, sep = "\t", row.names = F, col.names = F)

