# Lets play with my alignment: what ISs are in pAsa5, 01-B522??

setwd("Desktop/projet_special/")


#### FUNCTIONS ####
loadpAsa5 <- function(filename){
  myrawtable <- read.table(filename, stringsAsFactors = FALSE)
  colnames(myrawtable) <- c("tag", "pAsa5", "refstart", "alRlength", "refstrand",
                         "contig", "cstart", "cRlength", "cstring", "clength")
  # I prefer start and stop
  refstop <- myrawtable$refstart + myrawtable$alRlength
  cstop <- myrawtable$cstart + myrawtable$cRlength
  
  myrawtable[, 4] <- refstop
  names(myrawtable[, 4]) <- "refstop"
  myrawtable[, 8] <- cstop
  names(myrawtable[, 8]) <- "cstop"
  
  return(myrawtable)
}
rawB522 <- read.table("01-B522/pAsa5.tab", stringsAsFactors = FALSE)

colnames(rawB522) <- c("tag", "pAsa5", "refstart", "alRlength", "refstrand",
                       "contig", "cstart", "cRlength", "cstring", "clength")
# I prefer start and stop
refstop <- rawB522$refstart + rawB522$alRlength
cstop <- rawB522$cstart + rawB522$cRlength

rawB522[, 4] <- refstop
names(rawB522[, 4]) <- c("refstop")
rawB522[, 8] <- cstop
names(rawB522[, 8]) <- c("cstop")

#### Formating done? ####
# Now I have to remove some alignments that I know don't fit
# Contig 135 is legit (close circle), contig 139 is second rep,
# contig 152 is prob pAsa9
goodB522 <- rawB522[rawB522$contig != "NZ_MIIM01000152.1", ]
goodB522 <- goodB522[-2, ]
# Ahha! B522 is a 454 sequencing! Reads are longer! IS are more covered!!!! So much
# it makes contigs overlap ?!?!?!?
# Contig 166 is ISAS5, but chromosomal
goodB522 <- goodB522[goodB522$contig != "NZ_MIIM01000166.1", ]
goodB522 <- goodB522[-6, ]
# Contig 153 is pAsa9
goodB522 <- goodB522[goodB522$contig != "NZ_MIIM01000153.1", ]
goodB522 <- goodB522[c(-9, -10), ]
goodB522 <- goodB522[c(-10, -11), ]
goodB522 <- goodB522[c(-16, -17, -18), ]
# Contig 155 is pAsa9, with low complexity gap of 200pb
goodB522 <- goodB522[goodB522$contig != "NZ_MIIM01000155.1", ]
goodB522 <- goodB522[goodB522$contig != "NZ_MIIM01000154.1", ]



allpAsa5 <- goodB522

library(ggplot2)
a <- ggplot(goodB522, aes(y = (mean(c(refstart, refstop))), x = tag, fill = tag))
a + geom_crossbar(aes(ymin = refstart, ymax = refstop, width = 0.1))



#### Good, another? 2009-157-K5####
raw157K5 <- loadpAsa5("2009-157K5/pAsa5.tab")
# There is potentially an insertion at begin contig96(in repA???)
# Potential deletion between 5b and 3b
# Remove: NZ_MIIQ01000097.1 (solo ISAS5), 
# Remove: row 7, 9
good157K5 <- raw157K5[c(-7, -9), ]
good157K5 <- good157K5[good157K5$contig != "NZ_MIIQ01000097.1", ]
a <- ggplot(good157K5, aes(y = (mean(c(refstart, alRlength))), x = tag, fill = tag))
a + geom_crossbar(aes(ymin = refstart, ymax = alRlength, width = 0.1))
# I miss some seq between traI and ISAS11B. Gotta look for them