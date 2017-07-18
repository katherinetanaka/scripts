*########################################################
#### Making heatmaps from tabular tBlastn results    ####
*########################################################
library(ggplot2)

#### Attribution ####
# This script is part of the following publication:
# Tanaka KH, Vincent AT, Trudel MV, Paquet VE, Frenette M, Charette SJ:
# The mosaic architecture of Aeromonas salmonicida subsp. salmonicida
# pAsa4 plasmid and its consequences on antibiotic resistance.
# PeerJ. 2016;4:e2595.

# and is reproduced according to its Creative Commons Licence
# https://creativecommons.org/licenses/by/4.0/

# Note from creator:
#  CDS were batch-submitted to NCBI website for tBlastn analysis using 
#  nr/nt database, resulting in multiple csv files. Files were merged by hand
#  and various analyses/experimentations were done to achieve a consistent and
#  readable heatmap (ie: kmeans groups, criteria for filtration). Only 
#  relevant code lines were extracted to write this document.

#### Steps that were done manually (ordering x axis) ####
xaxispAsa4 <- c("ASA_P4G049", "ASA_P4G048", "ASA_P4G047", "ASA_P4G046", 
           "ASA_P4G045",   
        "ASA_P4G044", "ASA_P4G043", "ASA_P4G042", "ASA_P4G041", "ASA_P4G040",   
        "ASA_P4G038", "ASA_P4G037", "ASA_P4G036", "ASA_P4G035", "ASA_P4G034",   
        "ASA_P4G033", "ASA_P4G032", "ASA_P4G031", "ASA_P4G030", "ASA_P4G029",   
        "ASA_P4G028", "ASA_P4G027", "ASA_P4G026", "ASA_P4G025", "ASA_P4G024",   
        "ASA_P4G023", "ASA_P4G022", "ASA_P4G021", "ASA_P4G020", "ASA_P4G019",   
       "ASA_P4G017", "ASA_P4G016", "ASA_P4G015", "ASA_P4G014", "ASA_P4G013",   
       "ASA_P4G012", "ASA_P4G011-10", "ASA_P4G009", "ASA_P4G008", "ASA_P4G007",            "ASA_P4G006", "ASA_P4G005", "ASA_P4G004", "ASA_P4G003", "ASA_P4G178",   
        "ASA_P4G177", "ASA_P4G176", "ASA_P4G175", "ASA_P4G174", "ASA_P4G173",   
       "ASA_P4G172", "ASA_P4G171", "ASA_P4G170", "ASA_P4G169", "ASA_P4G168",   
       "ASA_P4G166", "ASA_P4G165", "ASA_P4G152", "ASA_P4G153", "ASA_P4G154",   
        "ASA_P4G155", "ASA_P4G156", "ASA_P4G157", "ASA_P4G158", "traK",   
        "ASA_P4G160", "ASA_P4G161", "ASA_P4G162", "ASA_P4G163", "ASA_P4G151",   
       "ASA_P4JF001",   "ASA_P4JF002",   "ASA_P4JF003",   "ASA_P4JF004",   "ASA_P4JF005",  
       "ASA_P4JF006",   "ASA_P4JF007",   "ASA_P4JF008",   "ASA_P4JF009",   "ASA_P4B522.1", 
        "ASA_P4B522.3",  "ASA_P4B522.4",  "ASA_P4JF012",   "hypothetical",  "phage",  
       "putative",   "ASA_P4G150", "ASA_P4G149", "ASA_P4G148", "ASA_P4G147",   
       "ASA_P4G146", "ASA_P4G145", "ASA_P4G144", "ASA_P4G143", "ASA_P4G142",   
        "ASA_P4G141", "ASA_P4G139", "ASA_P4G138", "ASA_P4G137", "ASA_P4G136",   
       "ASA_P4G135", "ASA_P4G134", "ASA_P4G133", "ASA_P4G132", "ASA_P4G131",   
       "ASA_P4G130", "ASA_P4G129", "ASA_P4G128", "ASA_P4G127", "ASA_P4G126",   
        "ASA_P4G124", "ASA_P4G123", "ASA_P4G122", "ASA_P4G121", "ASA_P4G120",   
       "ASA_P4G119", "ASA_P4G118", "ASA_P4G117", "ASA_P4G116", "ASA_P4G115",            "ASA_P4G114", "ASA_P4G113", "ASA_P4G112", "ASA_P4G111", "ASA_P4G110",   
         "ASA_P4G109", "ASA_P4G108", "ASA_P4G107", "ASA_P4G106", "ASA_P4G105",   
        "ASA_P4G104", "ASA_P4G103", "ASA_P4G100", "ASA_P4G099", "ASA_P4G098",   
       "ASA_P4G097", "ASA_P4G096", "ASA_P4G095", "ASA_P4G094", "ASA_P4G093",   
       "ASA_P4G092", "ASA_P4G091", "ASA_P4G090", "ASA_P4G089", "ASA_P4G088",   
        "ASA_P4G087", "ASA_P4G002", "ASA_P4G001", "ASA_P4G081", "ASA_P4G080",   
        "ASA_P4G079", "ASA_P4G078", "ASA_P4G077", "ASA_P4G076", "ASA_P4G075",   
       "ASA_P4G074", "ASA_P4G073", "ASA_P4G072", "ASA_P4G071", "ASA_P4G070",   
       "ASA_P4G069", "tus",     "ASA_P4G067", "ASA_P4G066", "ASA_P4G065",   
       "ASA_P4G064", "ASA_P4G063", "ASA_P4G062", "ASA_P4G061", "ASA_P4G060",   
       "ASA_P4G059", "ASA_P4G058", "ASA_P4G057", "ASA_P4G056", "ASA_P4G055",   
       "ASA_P4G054", "ASA_P4G053", "ASA_P4G052", "ASA_P4G051", "ASA_P4G050")

#### tBlastnAll: function for dataframe creation from multiple csv file ####
tBlastnAll <- function(directory, ...){
  
  dir = paste0("./", directory)
  
  #### Exception handling ####
    # directory must be a character string, belongs to directory list
  if (!dir %in% list.dirs()) {
    stop("Argument directory must be a sub-directory of working one.")
  }
  #### End exception handling ####
  
  #### Main ####
  ## Get all files collided into one dataframe
  files <- list.files(path = dir, pattern = "*.csv")
  setwd(dir)
  temp.table <- do.call("rbind",
                    lapply(files,
                          function(x) read.csv(x, stringsAsFactors = FALSE,
                                               header = FALSE)))
  setwd("..")
  
  ## Data purification
  # Use strsplit to separate GI from Accession number
  subject <- strsplit(temp.table[, 2], split = "|", fixed = TRUE)
  
  # Use Accession numbers instead
  accession.true <- sapply(subject, FUN="[", 4)
  
  # Simplification 
  temp.table.s <- cbind(temp.table[ , c(1, 3, 4, 12, 13)],
                        accession.true, stringsAsFactors = FALSE)
  colnames(temp.table.s) <- c("query.id", "identity", "positive", "qstart",
                              "qend", "Accession")
  
  ## Get best entry by sequence and GI combined
  names <- unique(temp.table.s$query.id)
  final.table <- as.data.frame(NULL)

  for(i in 1:length(names)) {
    subset <- temp.table.s[temp.table.s$query.id == names[i], ]
    final.table <- rbind(final.table, subset[!duplicated(subset$Accession), ])
}
  out <- list(table = final.table)
  class(out) <- "tBlastnRAW"
  return(out)
}
#### End of make the all table method with best query.id/Accession combo ####

#### Filter: make remove poor match ####
filter <- function(data, minMatch) {
    #### Exception handling ####
    # must use a tBlastnRAW class
    if (class(data) != "tBlastnRAW") {
      stop("Use tBlastnAll to merge tables together.")
    }
    if (typeof(minMatch) != "integer") {
      stop("Enter an integer at minMatch")
    }
    #### End exception handling ####
    mytable <- data[["table"]]
    # Filtering
    entries <- rep(NA, length(unique(mytable$Accession)))
    names(entries) <- unique(mytable$Accession)
    
    for (i in 1:length(entries)) {
      nb.match <- length(mytable$Accession[
        mytable$Accession %in% names(entries)[i]])
      
      if (nb.match >= minMatch) {
        entries[i] <- nb.match
      }
    }
    entries.less <- entries[!is.na(entries)]
    table.cut <- mytable[mytable$Accession %in% names(entries.less), ]
    out <- list(table = table.cut, minMatch = minMatch)
    class(out) <- "tBlastnFilter"
    return(out)
}
#### End of Filter ####

#### Cluster: use kmeans to order y axis ####
clusterK <- function(data, Ngroups) {
  #### Exception handling ####
  # must use a tBlastnFilter class
  if (class(data) != "tBlastnFilter") {
    stop("Use filter to remove poorly represented Acc. no.")
  }
  if (typeof(Ngroups) != "integer") {
    stop("Enter an integer at Ngroups")
  }
  #### End exception handling ####
  dfLong <- data$table
  # Make table to kmeans
  tableAdist <- xtabs(identity ~ Accession + query.id, data = dfLong)
  trable <- matrix(tableAdist, nrow = attr(tableAdist, "dim"),
                   ncol = attr(tableAdist, "dim"))
  rownames(trable) <- attr(tableAdist, "dimnames")$Accession
  
  # kmeans
  k <- kmeans(trable, Ngroups, nstart = 40)
  tres.df <- data.frame(Accession = names(k$cluster), group = k$cluster)
  tres.df <- tres.df[order(tres.df$group), ]
  yaxis <- factor(tres.df$Accession, levels = tres.df$Accession)
  
  dfLong.group <- merge(dfLong, tres.df, by = "Accession", sort = FALSE)
  dfLong.group$Accession <- factor(dfLong.group$Accession,
                                   levels = levels(yaxis))
  dfLong.group$group <- factor(dfLong.group$group, levels = 1:Ngroups)
  out <- list(table = dfLong.group, yaxis = yaxis, Ngroups = Ngroups)
  class(out) <- "TBprocessed"
  return(out)
}
#### End of functions ####

# For pAsa4
mypAsa4 <- tBlastnAll("tblastn_pAsa4")
FpAsa4 <- filter(mypAsa4, 4L)
FpAsa4$table$query.id <- factor(FpAsa4$table$query.id, levels = xaxispAsa4)

PpAsa4 <- clusterK(FpAsa4, 5L) # Due to kmeans alloting random order to
# groups, results may vary

a <- ggplot(PpAsa4$table, aes(x = query.id,
                             y = Accession, fill = group, alpha = identity))
a + geom_tile() +
  scale_x_discrete("ORFs", labels = NULL) +
  scale_y_discrete("Hits", labels = NULL) +
  theme(panel.background = element_rect(fill = "white"))

