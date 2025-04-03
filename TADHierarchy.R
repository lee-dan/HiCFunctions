#!/usr/bin/env Rscript
# TADHierarchy.R
# This program generates a hierarchy list for each TAD, detailing its child TADs based on a specified overlap threshold.
# Author: Daniel Lee
# Date: April 2, 2025
# Version: 1.0.0

# Load required packages with proper error handling
required_packages <- c("GenomicRanges", "magrittr", "dplyr", "progress", "data.table")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required but not installed. Please install it first."))
  }
}

# Suppress warnings and package startup messages
suppressWarnings({
  suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts = FALSE))
  library(magrittr, warn.conflicts = FALSE)
  library(dplyr, warn.conflicts = FALSE)
  library(progress, warn.conflicts = FALSE)
  library(data.table, warn.conflicts = FALSE)
})

#' Parse command line arguments
#' @return A list of parsed arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    stop("No arguments provided. Use --help for usage information.")
  }
  
  # Check for help flag
  if ("--help" %in% args || "-h" %in% args) {
    cat("Usage: Rscript TADHierarchy.R --inputFile <file> --outputDirectory <dir> --threshold <value>\n")
    cat("\nArguments:\n")
    cat("  --inputFile: The .bedpe file containing the list of TADs\n")
    cat("  --outputDirectory: Directory where the TAD hierarchy will be saved\n")
    cat("  --threshold: Minimum overlap percentage [0 to 1] for a TAD to be considered a 'child' TAD (recommended: 1.0)\n")
    cat("  --help, -h: Show this help message\n")
    quit(status = 0)
  }
  
  argsList <- paste(unlist(args), collapse = ' ')
  listoptions <- unlist(strsplit(argsList, '--'))[-1]
  options.args <- sapply(listoptions, function(x) {
    unlist(strsplit(x, ' '))[-1]
  }, simplify = FALSE)
  options.names <- sapply(listoptions, function(x) {
    option <- unlist(strsplit(x, ' '))[1]
  })
  names(options.args) <- unlist(options.names)
  
  # Validate required arguments
  required_args <- c("inputFile", "outputDirectory", "threshold")
  missing_args <- required_args[!required_args %in% names(options.args)]
  if (length(missing_args) > 0) {
    stop(paste("Missing required arguments:", paste(missing_args, collapse = ", ")))
  }
  
  return(options.args)
}

#' Validate input parameters
#' @param args List of parsed arguments
#' @return Validated arguments
validate_args <- function(args) {
  # Validate input file
  if (!file.exists(args$inputFile)) {
    stop(paste("Input file does not exist:", args$inputFile))
  }
  
  # Validate threshold
  threshold <- as.double(args$threshold)
  if (is.na(threshold) || threshold < 0 || threshold > 1) {
    stop("Threshold must be a number between 0 and 1")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(args$outputDirectory)) {
    dir.create(args$outputDirectory, recursive = TRUE)
    message(paste("Created output directory:", args$outputDirectory))
  }
  
  return(list(
    inFile = args$inputFile,
    outDir = args$outputDirectory,
    threshold = threshold
  ))
}

#' Create TAD class
create_tad_class <- function() {
  setClass("TAD", slots = list(
    seqname = "character", 
    start = "numeric", 
    end = "numeric", 
    index = "character",
    children = "character"
  ))
}

#' Extract TAD information from input file
#' @param input_file Path to the input file
#' @return Data frame with TAD information
extract_tad_info <- function(input_file) {
  tryCatch({
    TADList <- read.table(file = input_file, sep = "\t")
    TADList <- TADList %>%
      select(V1, V2, V3)
    TADList$index <- as.character(1:nrow(TADList))
    colnames(TADList) <- c("seqnames", "start", "end", "index")
    return(TADList)
  }, error = function(e) {
    stop(paste("Error reading input file:", e$message))
  })
}

#' Process TADs for a specific chromosome
#' @param TADsOfSameChromosome Data frame of TADs for a specific chromosome
#' @param threshold Overlap threshold
#' @return List of TAD objects with hierarchy information
process_chromosome_tads <- function(TADsOfSameChromosome, threshold) {
  # Sort by widths (descending)
  TADsOfSameChromosome <- TADsOfSameChromosome[order(TADsOfSameChromosome$end - TADsOfSameChromosome$start, decreasing = TRUE),]
  
  # Create TAD hierarchy array for current chromosome
  TADhierachyLevel <- list()
  
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
                  children = ".")
    
    # If TAD hierarchy is not empty, iterate through list of TADs and 
    # determine whether overlap ratio is above the threshold value
    if (length(TADhierachyLevel) > 0) {
      for (k in 1:length(TADhierachyLevel)) {
        otherTADstart <- TADhierachyLevel[[k]]@start
        otherTADend <- TADhierachyLevel[[k]]@end
        
        overlapStart <- max(TADstart, otherTADstart)
        overlapEnd <- min(TADend, otherTADend)
        
        overlapDist <- overlapEnd - overlapStart
        
        # Determine overlap proportion
        overlapProportion <- overlapDist / (TADend - TADstart)
        
        # Check whether overlap proportion is greater than threshold
        if (overlapProportion >= threshold) {
          # Append current TAD to children of parent TAD
          if (TADhierachyLevel[[k]]@children == ".") {
            TADhierachyLevel[[k]]@children <- paste("TAD_", TADindex, sep = "")
          } else {
            TADhierachyLevel[[k]]@children <- paste(TADhierachyLevel[[k]]@children, 
                                                    "; ", "TAD_", TADindex, sep = "")
          }
        }
      }
    }
    
    # Append current TAD to TAD hierarchy of current chromosome
    TADhierachyLevel <- TADhierachyLevel %>% append(TADobj)
  }
  
  return(TADhierachyLevel)
}

#' Convert TAD hierarchy to dataframe
#' @param TADhierarchy List of TAD objects
#' @return Data frame with TAD hierarchy information
convert_hierarchy_to_df <- function(TADhierarchy) {
  TADhierarchyAsDF <- data.frame(matrix(ncol = 5, nrow = 0)) 
  for (tad in TADhierarchy) {
    TADhierarchyAsDF[nrow(TADhierarchyAsDF) + 1,] <- c(
      as.character(tad@seqname),
      tad@start,
      tad@end,
      tad@index,
      tad@children
    )
  }
  colnames(TADhierarchyAsDF) <- c("seqnames", "start", "end", "index", "children")
  return(TADhierarchyAsDF)
}

#' Save TAD hierarchy to BEDPE file
#' @param TADhierarchyDF Data frame with TAD hierarchy information
#' @param output_file Path to the output file
save_hierarchy_to_bedpe <- function(TADhierarchyDF, output_file) {
  # Convert to .bedpe file format
  finalTADDF <- data.frame(
    TADhierarchyDF$seqnames,
    TADhierarchyDF$start,
    TADhierarchyDF$end,
    TADhierarchyDF$seqnames,
    TADhierarchyDF$start,
    TADhierarchyDF$end,
    TADhierarchyDF$index,
    TADhierarchyDF$children
  )
  finalTADDF$TADhierarchyDF.index <- paste("TAD_", finalTADDF$TADhierarchyDF.index, sep = "")
  
  # Write header
  writeLines(
    paste("#chr1", "\t", "x1", "\t", "x2", "\t", "chr2", "\t", "y1", "\t", "y2", "\t", "identifier", "\t", "children"),
    output_file
  )
  
  # Write data
  write.table(
    finalTADDF, 
    output_file, 
    row.names = FALSE, 
    append = TRUE,
    sep = "\t", 
    col.names = FALSE, 
    quote = FALSE
  )
  
  message(paste("Saved TAD hierarchy to:", output_file))
}

#' Main function
main <- function() {
  # Parse and validate arguments
  args <- parse_args()
  params <- validate_args(args)
  
  # Create TAD class
  create_tad_class()
  
  # Extract TAD information
  message("Reading TAD information from input file...")
  TADList <- extract_tad_info(params$inFile)
  
  # Create TAD hierarchy array containing each TAD with its respective "children" TADs
  TADhierarchy <- list()
  
  # List of unique chromosomes in TAD file
  uniqueChromosomes <- unique(TADList$seqnames)
  
  # Progress Bar
  pb <- progress_bar$new(
    format = paste("Processing TADs: [:bar] :percent eta: :eta"),
    total = length(uniqueChromosomes),
    clear = FALSE
  )
  
  # Process each chromosome
  for (i in 1:length(uniqueChromosomes)) {
    pb$tick()
    
    # Data frame of all TADs of current chromosome
    TADsOfSameChromosome <- TADList[TADList$seqnames == uniqueChromosomes[i],]
    
    # Process TADs for this chromosome
    TADhierachyLevel <- process_chromosome_tads(TADsOfSameChromosome, params$threshold)
    
    # Append TAD hierarchy of current chromosome to list of all TAD hierarchies
    TADhierarchy <- TADhierarchy %>% append(TADhierachyLevel)
  }
  
  # Convert TAD hierarchy to dataframe
  message("Converting TAD hierarchy to dataframe...")
  TADhierarchyAsDF <- convert_hierarchy_to_df(TADhierarchy)
  
  # Convert TAD hierarchy dataframe to G Ranges object for sorting
  TADhierarchyGR <- makeGRangesFromDataFrame(
    TADhierarchyAsDF,
    keep.extra.columns = TRUE,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end"
  )
  TADhierarchyGR <- sortSeqlevels(TADhierarchyGR)
  TADhierarchyGR <- sort(TADhierarchyGR)
  
  # Save TAD hierarchy
  output_file <- file.path(params$outDir, "TADHierarchy.bedpe")
  save_hierarchy_to_bedpe(as.data.frame(TADhierarchyGR), output_file)
  
  message("TAD hierarchy generation completed successfully")
}

# Run the main function with error handling
tryCatch({
  main()
}, error = function(e) {
  message("Error: ", e$message)
  quit(status = 1)
})
