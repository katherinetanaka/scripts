*###############################################################################
#### Multiple alignment to variant base call -- pAsa4 ##########################
# Katherine Tanaka, May 12th 2015 ##

#### Attribution ####
# This script is part of the following publication:
# Tanaka KH, Vincent AT, Trudel MV, Paquet VE, Frenette M, Charette SJ:
# The mosaic architecture of Aeromonas salmonicida subsp. salmonicida
# pAsa4 plasmid and its consequences on antibiotic resistance.
# PeerJ. 2016;4:e2595.

# and is reproduced according to its Creative Commons Licence
# https://creativecommons.org/licenses/by/4.0/

library(zoo)
library(reshape2)
library(ggplot2)

## Use with markx0 files, it's good!

stretcher <- readLines("stretcher/stretcher_output")

# Separate stretcher in head (the information) and Seq (the sequences)
endHead <- grep("#={1,250}", stretcher)[2] # To get end of last boxe
stretcherH <- stretcher[1:endHead]
stretcherSeq <- stretcher[(endHead+1):length(stretcher)]

# Extract sequences from stretcherSeq
find <- regexpr("\\s[ATGCatgc-]{1,50}$", stretcherSeq)
see <- regmatches(stretcherSeq, find)

extract <- strsplit(see, " ", fixed = TRUE)

# All odds: sequence 1. All evens, sequence 2, + get rid of the space
seq1 <- sapply(extract[seq(from = 1, to = length(extract)-1, by = 2)],
               FUN = "[[", 2)
seq2 <- sapply(extract[seq(from = 2, to = length(extract), by = 2)],
               FUN = "[[", 2)

seq1 <- paste0(seq1, collapse = "")
seq2 <- paste0(seq2, collapse = "")

# as seq1-track: 
mutDF <- data.frame(strsplit(seq1, "", fixed = TRUE),
                    strsplit(seq2, "", fixed = TRUE),
                    Type = rep(NA, length(seq1)), stringsAsFactors = FALSE,
                    length = rep(NA, length(seq1)))
names(mutDF) <- c("Name1", "Name2", "Type", "length")

##### Insertions are the trickiest ones ####
ins <- grep("-", mutDF$Name1, fixed = TRUE)
names(ins) <- rep(NA, length(ins))

# How to get insertion by groups
i <- 1
names(ins)[1] <- i

for(j in 2:length(ins)){
  if(ins[j] - ins[j-1] == 1){
    names(ins)[j] <- i
  } else {
    i <- i + 1
    names(ins)[j] <- i
  }
}

# Do notation for each group. // i is max group number now
for(k in 1:i){
  mutDF[(ins[names(ins) == k][1])-1, 3] <- "ins"
  mutDF[(ins[names(ins) == k][1])-1, 4] <- length(ins[names(ins) == k])
}

# We can now shrink Seq1 to its original length, no insertion
mutDFclean <- mutDF[-ins, ]

##### Now the deletions ####
del <- grep("-", mutDFclean$Name2, fixed = TRUE)
names(del) <- rep(NA, length(del))

l <- 1
names(del)[1] <- l

for(m in 2:length(del)){
  if(del[m] - del[m-1] == 1){
    names(del)[m] <- l
  } else {
    l <- l + 1
    names(del)[m] <- l
  }
}

for(n in 1:l){
  mutDFclean[(del[names(del) == n]), 3] <- "del"
  mutDFclean[(del[names(del) == n][1]), 4] <- length(del[names(del) == n])
  
}

#### Substitutions ####  # By the way ins + del = total - sub
sub <- mutDFclean$Name1 == mutDFclean$Name2 | mutDFclean$Name2 == "-"
whichSub <- grep("FALSE", sub, fixed = TRUE)

# Build bank for transition, transversion is everything else, since no "-"
transition <- c("AG", "GA", "CT", "TC")
# Write transition and transversion
for(p in 1:length(whichSub)){
  if(is.na(mutDFclean[whichSub[p], 3])){
    if(paste0(mutDFclean[whichSub[p], 1], mutDFclean[whichSub[p], 2]) %in% 
         transition){
      mutDFclean[whichSub[p], 3] <- "transi"
    } else {
      mutDFclean[whichSub[p], 3] <- "transver"
    }
  }
}

#### Preparing for ggplot ####
mutDFclean$Type <- factor(mutDFclean$Type)
mutDFclean <- cbind(mutDFclean, 1:length(mutDFclean$Name1))
names(mutDFclean) <- c("Name1", "Name2", "Type", "length", "position")

# Lets use sliding windows

test <- rollapply(mutDFclean$Type, width = 1000, FUN = function(x) length(grep("transi", x, fixed = TRUE)))
test2 <- rollapply(mutDFclean$Type, width = 1000, FUN = function(x) length(grep("transver", x, fixed = TRUE)))
test3 <- rollapply(mutDFclean$Type, width = 1000, FUN = function(x) length(grep("del", x, fixed = TRUE)))
test4 <- rollapply(mutDFclean$Type, width = 1000, FUN = function(x) length(grep("ins", x, fixed = TRUE)))

sliding <- cbind(501:(length(test)+500), test, test2, test3, test4)
sliding <- data.frame(sliding)
names(sliding) <- c("window", "transitions", "transversions",
                    "deletion", "insertion")
slidingGG <- melt(sliding, id.vars = "window",
                  measure.vars = c("transitions", "transversions", "insertion"))

### Sliding gg2

names(test3) <- 501:(length(test3)+500)
deletionGG <- data.frame(window = names(test3), values = test3)
slidingGG2 <- merge(slidingGG, deletionGG, by = "window")

# Better: draw big deletions
bigdel <- test3[test3 > 500]

f <- ggplot(data = slidingGG2, aes(x = window, y = value, 
                                   color = factor(variable), 
                                   alpha = 5/(values+1)))

f + 
  geom_line() +
  scale_x_continuous(breaks=seq(0, 182000, by = 10000)) +
  ylim(c(0, 100)) 
theme(panel.background = element_rect(fill = "white", color = "grey20"),
      panel.grid.major = element_line(color = "grey70"),
      panel.grid.minor = element_line(color = "grey70")) 

