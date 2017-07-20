# November 21st, 2016
# Partial vs complete ISs in 01-B526
# Used command lastal -f TAB db FILE > output

setwd("~/Desktop/IS_count_01-B526/")
tabRaw <- read.table("all_vs_IS.tab", stringsAsFactors=FALSE)
names(tabRaw) <- c("score", "db", "startDB", "alignedDB", "strandDB", "lengthDB",
                   "query", "startQ", "alignedQ", "strandQ", "lengthQ",
                   "string", "EG", "Evalue")

ISAS5 <- tabRaw[tabRaw$db == 'ISAS5-CP000646_10', ]

sortedISAS5 <- ISAS5[order(ISAS5$alignedDB, decreasing = TRUE), ]

# On peut le faire systématiquement pour toutes les IS en même temps
percent <- tabRaw$alignedDB/tabRaw$lengthDB
tabRaw$percent <- percent

factdb <- factor(tabRaw$db)
tabRaw$factdb <- factdb

# Comment on fait pour compter les totales et les partielles?
tabRaw <- tabRaw[order(tabRaw$factdb), ]

TotalPartial <- data.frame(IS = levels(factdb), complete = NA, partial = NA)
TPvector <- ifelse(tabRaw$percent > 0.9, "complete", "partial")
tabRaw$TP <- factor(TPvector)

TotalPartial <- (tabRaw$factdb, tabRaw$TP)
# Ngggg

Distri <- tapply(tabRaw$percent, INDEX = tabRaw$factdb, FUN = summary)

# Ouin, mon premier graphique était pas mal, mais il faut peut-être enlever le
# bruit

tabNotRaw <- tabRaw[tabRaw$db %in% c('ISAS10-CP000646_7', 'ISAS11-pAsal1',
                                     'ISAs19', 'ISAs21', 'ISAS3-ISAs5',
                                     'ISAs33-ISAS6', 'ISAs34', 'ISAS35-ISBst12-like',
                                     'ISAS4-ISAs32', 'ISAS5-CP000646_10', 
                                     'ISAS6', 'ISAS8'), ]
sorted526 <- tabRaw[order(tabRaw[,2], tabRaw[,7], tabRaw[,8]), ]
write.csv(sorted526, file = "sorted526.csv")
# Visualisation
library(ggplot2)
library(RColorBrewer)

a <- ggplot(tabNotRaw, aes(x = percent))
a + geom_density() + facet_grid(factdb ~ .)

### 13 juillet 2017: plus de visualisation
b <- ggplot(tabNotRaw, aes(x = percent, fill = factdb))
b + geom_histogram(bins = 20) +
  scale_fill_brewer(palette = "Paired", labels = c("ISAS10", "ISAS11", "ISAs19", "ISAs21", "ISAS3", "ISAS6", "ISAs34",
                                                   "ISAS35", "ISAS4", "ISAS5", "ISAS6", "ISAS8"),
                    guide =  guide_legend(title = "IS")) +
  labs(x = "pourcentage de l'IS aligné", y = "Nombre d'IS")


c <- ggplot(tabRaw, aes(x = percent, fill = factdb))
c + geom_histogram(bins = 20, na.rm = TRUE)

#### A449! ####
setwd("~/Desktop/IS_count_01-B526/")
tabRaw449 <- read.table("A449_vs_IS.tab", header=FALSE, stringsAsFactors=FALSE)
names(tabRaw449) <- c("score", "db", "startDB", "alignedDB", "strandDB", "lengthDB",
                   "query", "startQ", "alignedQ", "strandQ", "lengthQ",
                   "string", "EG", "Evalue")

# On peut le faire systématiquement pour toutes les IS en même temps
percent <- tabRaw449$alignedDB/tabRaw449$lengthDB
tabRaw449$percent <- percent

factdb <- factor(tabRaw449$db)
tabRaw449$factdb <- factdb

tabRaw449 <- tabRaw449[order(tabRaw449$factdb), ]

TotalPartial <- data.frame(IS = levels(factdb), complete = NA, partial = NA)
TPvector <- ifelse(tabRaw449$percent > 0.9, "complete", "partial")
tabRaw449$TP <- factor(TPvector)

# Ouin, mon premier graphique était pas mal, mais il faut peut-être enlever le
# bruit (en ne gardant que les IS pour lesquelles 1 complète)

tabNotRaw449 <- tabRaw449[tabRaw449$db %in% c("ISAS1", "ISAS10-CP000646_7",
                                              "ISAS11-pAsal1", "ISAs19",
                                              "ISAS2", "ISAS3-ISAs5",
                                              "ISAs33-ISAS6", "ISAs34",
                                              "ISAS4-ISAs32", "ISAS5-CP000646_10",
                                              "ISAS6", "ISAS8", "ISAS9"), ]
sorted449 <- tabRaw449[order(tabRaw449[,2], tabRaw449[,7], tabRaw449[,8]), ]
write.csv(sorted449, file = "sortedA449.csv")
# Visualisation
library(ggplot2)
library(RColorBrewer)

a <- ggplot(tabNotRaw, aes(x = percent))
a + geom_density() + facet_grid(factdb ~ .)

### 13 juillet 2017: plus de visualisation
b <- ggplot(tabNotRaw, aes(x = percent, fill = factdb))
b + geom_histogram(bins = 20) +
  scale_fill_brewer(palette = "Paired", labels = c("ISAS10", "ISAS11", "ISAs19", "ISAs21", "ISAS3", "ISAS6", "ISAs34",
                                                   "ISAS35", "ISAS4", "ISAS5", "ISAS6", "ISAS8"),
                    guide =  guide_legend(title = "IS")) +
  labs(x = "pourcentage de l'IS aligné", y = "Nombre d'IS")


c <- ggplot(tabRaw, aes(x = percent, fill = factdb))
c + geom_histogram(bins = 20, na.rm = TRUE)
