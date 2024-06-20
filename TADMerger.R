# Daniel Lee
# Program to merge TADs at different resolutions

library(GenomicRanges)
library(magrittr)
library(dplyr)

#__________Merge TADs__________

# Receive input on TAD directory and resolutions 
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
inDir <- options.args$inputDirectory
resolutions <- options.args$resolutions
TADthreshold <- as.double(options.args$threshold)
outDir <- options.args$outputDirectory

for (resolution in resolutions) {
  resolution <- paste(resolution, "_blocks.bedpe")
}

# Create an empty "merged list" of unique TADs across all resolutions
uniqueTADs <- GRanges()

# Iterate through each resolution 
for (resolution in resolutions) {
  print(paste(resolution, "kb", sep=""))
  # Create a label for each resolution (i.e. "10kb_resolution", 
  # "25kb_resolution", etc.)
  shortFileName <- paste(sub("000_blocks.bedpe.*", "", resolution), "kb_resolution", sep="")
  
  # Read the data from the .bedpe file to be stored as a dataframe
  TADdf <- read.table(file=normalizePath(paste(inDir, "/", resolution, "_blocks.bedpe", sep="")), sep="\t")
  colnames(TADdf) <- c("chr1",	"x1", "x2", 
                       "chr2", "y1", "y2", 
                       "name",	"score", "strand1",	"strand2", "color", "score", "uVarScore", "lVarScore", "upSign", "loSign")
  # Correct for formatting
  TADdf$chr1[substr(TADdf$chr1, 1, 1) != "c"] <- paste("chr", TADdf$chr1, sep="") 
  
  # Create unique identifiers for each TAD, consisting of its resolution and its line number in original file
  TADidentifiers <- data.frame(matrix(ncol = 1, nrow = 0))
  colnames(TADidentifiers) <- "identifiers"
  for (i in 1:nrow(TADdf)) {
    TADidentifiers[nrow(TADidentifiers) + 1,] = paste(shortFileName, "_", i, "-", shortFileName, "_", i, sep="")
  }
  
  # Append metadata to TAD dataframe
  TADdf <- cbind(TADdf, TADidentifiers)
  
  # Iterate through each TAD
  for (i in 1:nrow(TADdf)) {
    
    # Extract information about the TAD
    TADseqname = TADdf[i, "chr1"]
    TADstart = TADdf[i, "x1"]
    TADend = TADdf[i, "x2"]
    TADidentifier = TADdf[i, "identifiers"]
    
    # Convert into Genomic Ranges object
    TADinfo <- data.frame(TADseqname, TADstart, TADend, TADidentifier)
    colnames(TADinfo) <- c("seqname", "start", "end", "identifier")
    TADasGR <- makeGRangesFromDataFrame(TADinfo,
                                        keep.extra.columns = TRUE,
                                        seqnames.field = "seqname",
                                        start.field = "start",
                                        end.field = "end")
    
    #________Determine whether or not the TAD is unique________
    
    # First, find all TADs in the list of merged TADs thus far that overlap with 
    # the given TAD
    overlaps <- findOverlaps(TADasGR, uniqueTADs, select = "all")
    
    # If there are no overlapping TADs, then this TAD is unique
    if (length(overlaps) == 0) {
      # Add TAD to list of merged TADs
      uniqueTADs <- uniqueTADs %>%
        append(TADasGR)
    } else {
      # Otherwise, iterate through the overlapping TADs:
      overlappingTADs <- subjectHits(overlaps)
      isTADUnique <- TRUE
      
      for (index in overlappingTADs) {
        # Check if overlap proportion of both TADs is less than the threshold
        otherTADstart <- start(uniqueTADs)[index]
        otherTADend <- end(uniqueTADs)[index]
        
        overlapStart <- max(TADstart, otherTADstart)
        overlapEnd = min(TADend, otherTADend)
        overlapDist <- overlapEnd - overlapStart
        
        overlapProportion1 <- overlapDist / (TADend - TADstart)
        overlapProportion2 <- overlapDist / (otherTADend - otherTADstart)
        
        if (overlapProportion1 > TADthreshold && overlapProportion2 > TADthreshold) {
          isTADUnique = FALSE
          
          # Merge TADs by taking the lower of the two start positions, 
          # and the higher of the two end positions
          start(uniqueTADs)[index] = min(TADstart, otherTADstart)
          if (TADstart < otherTADstart) uniqueTADs$identifier[index] = gsub(".*-", paste(sub("-.*", "", TADidentifier), "-", sep=""), uniqueTADs$identifier[index]) 
          
          end(uniqueTADs)[index] = max(TADend, otherTADend)
          if (TADend > otherTADend) uniqueTADs$identifier[index] = gsub("-.*", paste("-", sub(".*-", "", TADidentifier), sep=""), uniqueTADs$identifier[index])
          
          break
        }
      }
      
      # If no overlaps greater than [TADthreshold] were found, then the TAD is unique
      # and will be added to the list of unique TADs.
      if (isTADUnique) {
        uniqueTADs <- uniqueTADs %>%
          append(TADasGR)
      }
    }
  }
}
# Sort uniqueTADs
uniqueTADs <- sortSeqlevels(uniqueTADs)
uniqueTADs <- sort(uniqueTADs)

# Format as .bedpe file
uniqueTADsDF <- as.data.frame(uniqueTADs)
finalDF <- data.frame(uniqueTADsDF$seqnames,
                      uniqueTADsDF$start,
                      uniqueTADsDF$end, 
                      uniqueTADsDF$seqnames,
                      uniqueTADsDF$start,
                      uniqueTADsDF$end,
                      uniqueTADsDF$identifier)
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive=TRUE)
}
writeLines(paste("#chr1", "\t",
                 "x1", "\t",
                 "x2", "\t",
                 "chr2", "\t",
                 "y1", "\t",
                 "y2", "\t",
                 "name"),
           paste(outDir, "/mergedTADs.bedpe", sep="")
)

# Save list of uniqueTADs
write.table(finalDF, normalizePath(paste(outDir, "/mergedTADs.bedpe", sep="")), row.names=FALSE, col.names = FALSE, append = TRUE, quote = FALSE, sep = "\t")
