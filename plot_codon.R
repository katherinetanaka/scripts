#### Plot codon usage ####

# Notes: Input is a Kazusa-style codon usage table in which tabs are replaced 
# by newlines. pAsa4(b|c) tables are from
# all non pseudo CDS and created in Artemis. A449 table is from 
# http://www.kazusa.or.jp/codon/

prepCodon <- function(file, tag){
  table <- read.table(file, stringsAsFactors = FALSE)
  names(table) <- c("codon", "per1000", "count")
  table$freq <- as.numeric(sapply(sapply(table[[2]], FUN = strsplit, split = "(", fixed = TRUE),
                       FUN = "[[", 1))
  table$tag <- factor(tag)
  table <- table[, c(1,4,5)]
  return(table)
  
}

A449 <- prepCodon("Desktop/pAsa4_codon/A449_usage2", "A449")
pAsa4 <- prepCodon("Desktop/pAsa4_codon/pAsa4_usage2", "pAsa4")
pAsa4b <- prepCodon("Desktop/pAsa4_codon/pAsa4b_usage2", "pAsa4b")
pAsa4c <- prepCodon("Desktop/pAsa4_codon/pAsa4c_usage2", "pAsa4c")
pSN254b <- prepCodon("Desktop/pAsa4_codon/pSN254b_usage2", "pSN254b")
pAsa5 <- prepCodon("Desktop/pAsa4_codon/pAsa5b_usage2", "pAsa5")

A449$codon <- pAsa4$codon

all <- rbind(A449, pAsa4, pAsa4b, pAsa4c, pSN254b, pAsa5)


#### Drawing ####
library(ggplot2)
a <- ggplot(all, aes(x = codon, y = freq, fill = tag))
a + geom_bar(stat = "identity", position = "dodge")
