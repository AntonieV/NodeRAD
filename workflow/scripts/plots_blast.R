library("tidyverse")
library("stringr")

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

blast <- snakemake@input[[1]]
blast_data <- read_tsv(blast, col_names = TRUE, trim_ws = TRUE)
names(blast_data) <- c('res_id','sim_id','identity', 'len_alignment', 'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score')

blast_data <- blast_data %>%
  mutate(res_id=as.character(res_id)) %>%
  separate(res_id, c("sample", "individual", "locus"), sep = "\\|", extra="merge", convert=TRUE) #%>%

for(i in blast_data$locus) {
  if(nchar(strsplit(i, "-")[[1]][1])<5){
    blast_data$locus[blast_data$locus==i] <- str_replace(blast_data$locus[blast_data$locus==i], 'LOC', 'LOC0')
  }
}

plot <- ggplot(data = blast_data, aes(y=locus, x=identity)) +
  geom_bar(width = 1.0, position = "dodge", stat="identity", aes(fill = identity), colour="Black") +
  scale_fill_gradient2(low="mediumpurple3", mid="steelblue3", high="green3", midpoint = 99, name = "Identity") +
  ggtitle("Indentity [%] of loci identified by NodeRAD vs. simulated loci") +
  xlab("Identity [%]") +
  ylab("Locus") +
  theme_minimal() +
  theme(aspect.ratio = 2.5/1.5, plot.title = element_text(hjust = 0.5), legend.position = "right", legend.key.size = unit(0.8, "cm"), axis.text.y = element_text(hjust = 0))
plot
ggsave(snakemake@output[["ident"]], width = 7, height = 7)

plot<-qplot(x=blast_data$identity,
         fill=..count..,
         geom="histogram",
         binwidth = 1,
         main = "Histogram of correctly identified loci by NodeRAD",
         col=I("black"),
         xlab = "Identity [%]",
         ylab = "Count",
         alpha=I(.8)) +
  scale_fill_gradient(low="purple3", high="green3", name = "Count") +
  theme(aspect.ratio = 2.5/1.5, plot.title = element_text(hjust = 0.5), legend.position = "right")
plot
ggsave(snakemake@output[["ident_hist"]], width = 7, height = 7)

plot <- ggplot(data = blast_data, aes(x=locus, y=bit_score, group = individual)) +
  geom_line(color = "gray70") +
  geom_point(aes(color = bit_score), size =3) +
  geom_point(shape = 1,size = 3, colour = "black") +
  scale_color_gradient(low="red", high="green", name = "Bit score") +
  ggtitle("Bitscores of loci identified by NodeRAD vs. simulated loci") +
  xlab("Locus") +
  ylab("Bit score") +
  theme_minimal() +
  theme(axis.text.x = element_text(color = "black", size = 7, angle = 90, hjust = 0, face = "plain"), plot.title = element_text(hjust = 0.5), legend.position = "right", legend.key.size = unit(0.4, "cm"), axis.text.y = element_text(hjust = 0))
plot
ggsave(snakemake@output[["bit_scores"]], width = 7, height = 7)

plot <- ggplot(data = blast_data, aes(x=locus, y=evalue, group = individual)) +
  geom_line(color = "gray70") +
  geom_point(aes(color = evalue), size =3) +
  geom_point(shape = 1,size = 3, colour = "black") +
  scale_color_gradient(low="green", high="red", name = "E-value") +
  ggtitle("E-Values of loci identified by NodeRAD vs. simulated loci") +
  xlab("Locus") +
  ylab("E-Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(color = "black", size = 7, angle = 90, hjust = 0, face = "plain"), plot.title = element_text(hjust = 0.5), legend.position = "right", legend.key.size = unit(0.4, "cm"), axis.text.y = element_text(hjust = 0))
plot
ggsave(snakemake@output[["evalues"]], width = 7, height = 7)
