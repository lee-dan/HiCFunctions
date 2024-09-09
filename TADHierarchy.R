# This program generates a hierarchy list for each TAD, detailing each of its 'children' TADs.

# Format:
# Rscript TADHierarchy.R --inputFile ... --outputDirectory ... --threshold ...

# Args:
# --inputFile: The .bedpe file containing the list of TADs.
# --outputDirectory**: The directory in which the TAD hierarchy will be placed in.
# -threshold**: The minimum overlap percentage [0 to 1] for a TAD to be considered a 'child' TAD of another TAD.
# The recommended value is 1.0.

# Example:
# Rscript TADHierarchy.R --inputFile sampleTADs.bedpe --outputDirectory outDir --threshold 1.0

suppressWarnings({
  suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts = FALSE))
  library(magrittr, warn.conflicts = FALSE)
  library(dplyr, warn.conflicts = FALSE)
  library(progress, warn.conflicts = FALSE)
  library(data.table, warn.conflicts = FALSE)
})

# Read user-inputted parameters:
args <- commandArgs(trailingOnly = TRUE)
argsList <- paste(unlist(args), collapse = ' ')
listoptions <- unlist(strsplit(argsList, '--'))[-1]
options.args <- sapply(listoptions, function(x){
  unlist(strsplit(x, ' '))[-1]
}, simplify = FALSE)
options.names <- sapply(listoptions, function(x){
  option <- unlist(strsplit(x, ' '))[1]
})
names(options.args) <- unlist(options.names)
inFile <- options.args$inputFile
TADthreshold <- as.double(options.args$threshold)
outDir <- options.args$outputDirectory

# Create TAD class
setClass("TAD", slots = list(seqname = "character", 
                             start = "numeric", 
                             end = "numeric", 
                             index = "character",
                             children = "character"))

# Extract TAD information
TADList <- read.table(file = inFile, sep="\t")
TADList <- TADList %>%
  select(V1, V2, V3)
TADList$index <- as.character(1:nrow(TADList))
colnames(TADList) <- c("seqnames", "start", "end", "index")

# Create TAD hierarchy array containing each TAD with its respective
# "children" TADs
TADhierarchy <- list()

# Iterate through each chromosome to generate its respective TAD hierarchy
# List of unique chromosomes in TAD file
uniqueChromosomes <- unique(TADList$seqnames)

# Progress Bar:
pb <- progress_bar$new(format = paste("Processing TADs: [:bar] :percent eta: :eta"),
total = length(uniqueChromosomes),
clear = FALSE)

for (i in 1:length(uniqueChromosomes)) {
  pb$tick()

  # Create TAD hierarchy array for current chromosome to be appended to 
  # TAD hierarchy
  TADhierachyLevel <- list()
  
  # Data frame of all TADs of current chromosome
  TADsOfSameChromosome = TADList[TADList$seqnames == uniqueChromosomes[i],]
  
  # Sort by widths (descending)
  TADsOfSameChromosome <- TADsOfSameChromosome[order(TADsOfSameChromosome$end - TADsOfSameChromosome$start, decreasing=TRUE),]
  
  # Iterate through each TAD of current chromosome
  for (j in 1:nrow(TADsOfSameChromosome)) {
    
    # Extract information from TAD
    TADseqname <- TADsOfSameChromosome[j, "seqnames"]
    TADstart <- TADsOfSameChromosome[j, "start"]
    TADend <- TADsOfSameChromosome[j, "end"]
    TADindex <- TADsOfSameChromosome[j, "index"]
    
    # Form TAD object
    TADobj <- new("TAD", 
                  seqname = TADseqname, 
                  start = TADstart, 
                  end = TADend,
                  index = TADindex,
                  children=".")
    
    # If TAD hierarchy is not empty, iterate through list of TADs and 
    # determine whether overlap ratio is above the threshold value
    if (length(TADhierachyLevel) > 0) {
      for (k in 1:length(TADhierachyLevel)) {
        otherTADstart <- TADhierachyLevel[[k]]@start
        otherTADend <- TADhierachyLevel[[k]]@end
        
        overlapStart = max(TADstart, otherTADstart)
        overlapEnd = min(TADend, otherTADend)
        
        overlapDist <- overlapEnd - overlapStart
        
        # Determine overlap proportion
        overlapProportion <- overlapDist / (TADend - TADstart)
        
        # Check whether overlap proportion is greater than threshold
        if (overlapProportion >= TADthreshold) {
          
          # Append current TAD to children of parent TAD
          if (TADhierachyLevel[[k]]@children == ".") {
            TADhierachyLevel[[k]]@children <- paste("TAD_", TADindex, sep="")
          } else {
            TADhierachyLevel[[k]]@children <- paste(TADhierachyLevel[[k]]@children, 
                                                    "; ", "TAD_", TADindex, sep="")
          }
        }
      }
    }
    
    # Append current TAD to TAD hierarchy of current chromosome
    TADhierachyLevel <- TADhierachyLevel %>%
      append(TADobj)
  }
  # Append TAD hierarchy of current chromosome to list of all TAD hierarchies
  TADhierarchy <- TADhierarchy %>%
    append(TADhierachyLevel)
}

# Convert TAD hierarchy to dataframe
TADhierarchyAsDF <- data.frame(matrix(ncol = 5, nrow = 0)) 
for (tad in TADhierarchy) {
  TADhierarchyAsDF[nrow(TADhierarchyAsDF) + 1,] = c(as.character(tad@seqname),
                                                    tad@start,
                                                    tad@end,
                                                    tad@index,
                                                    tad@children)
}
colnames(TADhierarchyAsDF) <- c("seqnames", "start", "end", "index", "children")

# Convert TAD hierarchy dataframe to G Ranges object for sorting
TADhierarchyGR <- makeGRangesFromDataFrame(TADhierarchyAsDF,
                                           keep.extra.columns = TRUE,
                                           seqnames.field = "seqnames",
                                           start.field = "start",
                                           end.field = "end")
TADhierarchyGR <- sortSeqlevels(TADhierarchyGR)
TADhierarchyGR <- sort(TADhierarchyGR)

# Convert to .bedpe file format
TADhierarchyDF <- as.data.frame(TADhierarchyGR)
finalTADDF <- data.frame(TADhierarchyDF$seqnames,
                         TADhierarchyDF$start,
                         TADhierarchyDF$end,
                         TADhierarchyDF$seqnames,
                         TADhierarchyDF$start,
                         TADhierarchyDF$end,
                         TADhierarchyDF$index,
                         TADhierarchyDF$children)
finalTADDF$TADhierarchyDF.index <- paste("TAD_", finalTADDF$TADhierarchyDF.index, sep="")

# Saving TAD hierarchy
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive=TRUE)
}
fileName <- paste(outDir, "/TADHierarchy.bedpe", sep="")
writeLines(paste("#chr1", "\t",
                 "x1", "\t",
                 "x2", "\t",
                 "chr2", "\t",
                 "y1", "\t",
                 "y2", "\t",
                 "identifier",
                 "children"),
           paste(outDir, "/TADHierarchy.bedpe", sep="")
)
write.table(finalTADDF, 
            fileName, 
            row.names=FALSE, 
            append = TRUE,
            sep = "\t", 
            col.names = FALSE, 
            quote = FALSE)
