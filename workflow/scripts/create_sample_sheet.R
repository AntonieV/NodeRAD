library("tidyverse")
setwd("../../../NodeRAD")
### little script to create samples.tsv

sample <- data.frame(
  sample = c("A", "B", "C", "D"),
  individual = c("Individual_16", "Individual_20", "Individual_23", "Individual_24"),
  barcode_1 = c("CCGTCC", "GTTTCG",  "ACTGAT", "ATTCCT"),  # p5 barcodes if PE reads
  barcode_2 = c(rep("ATCACG", 4)), # p7 barcodes if PE reads
  spacer_1 = c("AC","", rep("C", 2)),  # p5 spacer/insert if PE reads
  spacer_2 = c(rep("", 4))  # p7 spacer/insert if PE reads
)

# write for config in .test
write_tsv(sample,"./.test/config/samples.tsv")

# write for global config
write_tsv(sample, "./config/samples.tsv")
