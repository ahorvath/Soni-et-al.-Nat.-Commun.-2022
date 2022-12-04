library(Rsamtools)
library(ggplot2)
library(dplyr)

sample.ids <- c("Sp_Red1DEL-FTP_polyARNA_NP_DirectRNA_111120", "Sp_Red1-FTP_polyARNA_NP_DirectRNA_111120", "Sp_Red1WFtoAA-FTP_polyARNA_NP_DirectRNA_111120", "Sp_Red1DEL-FTP_RNAIP_NP_DirectRNA_051120", "Sp_Red1-FTP_RNAIP_NP_DirectRNA_051120", "Sp_Red1WFtoAA-FTP_RNAIP_NP_DirectRNA_051120")
args <- commandArgs(T)
# predicted regions to be quantified
bed.file <- args[1]
bed.name <- gsub(".bed", "", bed.file) 
print(bed.name)

all.tail.table <- NULL
all.tail.table.na0 <- NULL
for (orig.sample.id in sample.ids) {
	print(orig.sample.id)
	#orig.sample.id <- "Sp_Red1DEL-FTP_polyARNA_NP_DirectRNA_111120"
	sample.id <- gsub("Sp_|_NP_DirectRNA_111120|_NP_DirectRNA_051120", "", orig.sample.id)
	tail.table.file <- paste0("../", orig.sample.id, "_tails.csv")
	region.table.file <- paste0(orig.sample.id, "_", bed.file)
	tail.table <- read.table(tail.table.file, sep = ",", header = T)[,c(1,5)]
	print(region.table.file)
	region.table <- read.table(region.table.file, sep = "\t", header = F)[,c(4,10)]
	colnames(region.table) <- c("region", "read_id")
	merged.table <- merge(tail.table, region.table, by = c("read_id"))
	sample.tail.table <- merged.table %>% group_by(region) %>% summarise(mean_tail_length = mean(tail_length, na.rm = T), median_tail_length = median(tail_length, na.rm = T))
	sample.tail.table <- sample.tail.table %>%  filter(!is.na(mean_tail_length))
	cat("sample.tail.table:", dim(sample.tail.table), "\n")
	cat("sample.tail.table:", summary(sample.tail.table), "\n")
	sample.tail.table.na0 <- merged.table %>% group_by(region) %>% mutate(tail_length = tidyr::replace_na(tail_length, 0)) %>% summarise(mean_tail_length = mean(tail_length), median_tail_length = median(tail_length))
	all.tail.table <- bind_rows(all.tail.table, bind_cols(Condition = sample.id, Selection = ifelse(length(grep("IP", sample.id)) == 1, "IP", "noIP"), sample.tail.table))
	all.tail.table.na0 <- bind_rows(all.tail.table.na0, bind_cols(Condition = sample.id, sample.tail.table.na0))
}

#gg <- ggplot(all.tail.table, aes(Condition, mean_tail_length, fill = Selection)) + geom_boxplot() + ylim(0, 200) +  geom_point(position = position_dodge(width = 0.75),aes(group = Selection))  + theme_bw()
#gg.violin <- ggplot(all.tail.table, aes(Condition, mean_tail_length, fill = Selection)) + geom_violin() + ylim(0, 200) +  geom_point(position = position_dodge(width = 0.75),aes(group = Selection)) + theme_bw()

#ggsave(gg, filename = paste0("gw_RED1KO_illumina_polyA_mean_boxplot_nat_",bed.name,"_y200.pdf"), width = 12)
#ggsave(gg.violin, filename = paste0("gw_RED1KO_illumina_mean_polyA_violin_nat_",bed.name,"_y200.pdf"), width = 12)

gg <- ggplot(subset(all.tail.table, Selection == "IP"), aes(Condition, mean_tail_length, fill = Selection)) + geom_boxplot() + ylim(0, 200) + geom_point(position = position_dodge(width = 0.75),aes(group = Selection))  + theme_bw()
gg.violin <- ggplot(subset(all.tail.table, Selection == "IP"), aes(Condition, mean_tail_length, fill = Selection)) + geom_violin() + ylim(0, 200) + geom_point(position = position_dodge(width = 0.75),aes(group = Selection)) + theme_bw()
write.table(subset(all.tail.table, Selection == "IP"), paste0("gw_RED1KO_illumina_polyA_mean_boxplot_nat_",bed.name,"_IP_y200.tsv"),sep = "\t", quote = F, row.names = F)
ggsave(gg, filename = paste0("gw_RED1KO_illumina_polyA_mean_boxplot_nat_",bed.name,"_IP_y200.pdf"), width = 7)
ggsave(gg.violin, filename = paste0("gw_RED1KO_illumina_mean_polyA_violin_nat_",bed.name,"_IP_y200.pdf"), width = 7)

