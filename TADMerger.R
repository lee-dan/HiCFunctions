#!/usr/bin/env Rscript
# TADMerger.R
# This program merges TADs at different resolutions based on a maximum allowed threshold of shared area.
# Author: Daniel Lee
# Date: April 2, 2025
# Version: 1.0.0

# Format:
# Rscript TADMerger.R --inputDirectory ... --outputDirectory ... --resolutions ... --threshold ...

# Args:
# --inputDirectory: The directory which contains the .bedpe files of TADs
# --outputDirectory: The directory in which the merged TADs will be placed in
# --resolutions: The resolutions of TADs to be merged [e.g. 10000, 25000, 50000, etc.]
# --threshold: If two TADs share more than [threshold] of their 'linear area' (distance from
# start index to end index), then the TADs will be 'merged'. The recommended value is 0.7 (i.e. 70%).
# E.g. The boundaries of TAD A are [0, 100,000], and the boundaries of TAD B are [25,000, 110,000].
# If [threshold] = 0.7, then TAD A and TAD B will be merged, forming a single TAD with boundaries [0, 110,000].
# Their shared area is [100,000] - [25,000] = [75,000]. [75,000] > 0.7 x [100,000] (the area of TAD A), and
# [75,000] > 0.7 x [85,000] (the area of TAD B).

# Example:
# Rscript TADMerger.R --inputDirectory inDir --resolutions 10000 25000 --threshold 0.7 --outputDirectory outDir

# Load required packages with proper error handling
required_packages <- c("GenomicRanges", "magrittr", "dplyr", "progress")
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
    cat("Usage: Rscript TADMerger.R --inputDirectory <dir> --outputDirectory <dir> --resolutions <values> --threshold <value>\n")
    cat("\nArguments:\n")
    cat("  --inputDirectory: The directory which contains the .bedpe files of TADs\n")
    cat("  --outputDirectory: The directory in which the merged TADs will be placed in\n")
    cat("  --resolutions: The resolutions of TADs to be merged [e.g. 10000, 25000, 50000, etc.]\n")
    cat("  --threshold: If two TADs share more than [threshold] of their 'linear area', then the TADs will be 'merged'\n")
    cat("  --help, -h: Show this help message\n")
    quit(status = 0)
  }

  argsList <- paste(unlist(args), collapse = " ")
  listoptions <- unlist(strsplit(argsList, "--"))[-1]
  options.args <- sapply(listoptions, function(x) {
    unlist(strsplit(x, " "))[-1]
  }, simplify = FALSE)
  options.names <- sapply(listoptions, function(x) {
    option <- unlist(strsplit(x, " "))[1]
  })
  names(options.args) <- unlist(options.names)

  # Validate required arguments
  required_args <- c("inputDirectory", "outputDirectory", "resolutions", "threshold")
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
  # Validate input directory
  if (!dir.exists(args$inputDirectory)) {
    stop(paste("Input directory does not exist:", args$inputDirectory))
  }

  # Validate threshold
  threshold <- as.double(args$threshold)
  if (is.na(threshold) || threshold < 0 || threshold > 1) {
    stop("Threshold must be a number between 0 and 1")
  }

  # Validate resolutions
  resolutions <- args$resolutions
  if (length(resolutions) == 0) {
    stop("At least one resolution must be specified")
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(args$outputDirectory)) {
    dir.create(args$outputDirectory, recursive = TRUE)
    message(paste("Created output directory:", args$outputDirectory))
  }

  return(list(
    inDir = args$inputDirectory,
    outDir = args$outputDirectory,
    resolutions = resolutions,
    threshold = threshold
  ))
}

#' Process a TAD file
#' @param file_path Path to the TAD file
#' @param resolution Resolution of the TAD file
#' @param uniqueTADs List of unique TADs
#' @param threshold Threshold for merging TADs
#' @return Updated list of unique TADs
process_tad_file <- function(file_path, resolution, uniqueTADs, threshold) {
  # Create a label for each resolution (i.e. "10kb_resolution",
  # "25kb_resolution", etc.)
  shortFileName <- paste(sub("000_blocks.bedpe.*", "", resolution), "kb_resolution", sep = "")

  # Read the data from the .bedpe file to be stored as a dataframe
  TADdf <- read.table(file = normalizePath(file_path), sep = "\t")
  colnames(TADdf) <- c(
    "chr1", "x1", "x2",
    "chr2", "y1", "y2",
    "name", "score", "strand1", "strand2", "color", "score", "uVarScore", "lVarScore", "upSign", "loSign"
  )

  # Correct for formatting
  TADdf$chr1[substr(TADdf$chr1, 1, 1) != "c"] <- paste("chr", TADdf$chr1, sep = "")

  # Create unique identifiers for each TAD, consisting of its resolution and its line number in original file
  TADidentifiers <- data.frame(matrix(ncol = 1, nrow = 0))
  colnames(TADidentifiers) <- "identifiers"
  for (i in 1:nrow(TADdf)) {
    TADidentifiers[nrow(TADidentifiers) + 1, ] <- paste(shortFileName, "_", i, "-", shortFileName, "_", i, sep = "")
  }

  # Append metadata to TAD dataframe
  TADdf <- cbind(TADdf, TADidentifiers)

  # Progress Bar:
  pb <- progress_bar$new(
    format = paste("Processing TADs at ", shortFileName, " [:bar] :percent eta: :eta", sep = ""),
    total = nrow(TADdf), clear = FALSE
  )
  pb$tick(0)

  # Iterate through each TAD
  for (i in 1:nrow(TADdf)) {
    pb$tick()

    # Extract information about the TAD
    TADseqname <- TADdf[i, "chr1"]
    TADstart <- TADdf[i, "x1"]
    TADend <- TADdf[i, "x2"]
    TADidentifier <- TADdf[i, "identifiers"]

    # Convert into Genomic Ranges object
    TADinfo <- data.frame(TADseqname, TADstart, TADend, TADidentifier)
    colnames(TADinfo) <- c("seqname", "start", "end", "identifier")
    TADasGR <- makeGRangesFromDataFrame(TADinfo,
      keep.extra.columns = TRUE,
      seqnames.field = "seqname",
      start.field = "start",
      end.field = "end"
    )

    # Determine whether or not the TAD is unique
    # First, find all TADs in the list of merged TADs thus far that overlap with
    # the given TAD
    suppressWarnings({
      overlaps <- findOverlaps(TADasGR, uniqueTADs, select = "all")
    })

    # If there are no overlapping TADs, then this TAD is unique
    if (length(overlaps) == 0) {
      # Add TAD to list of merged TADs
      suppressWarnings({
        uniqueTADs <- uniqueTADs %>%
          append(TADasGR)
      })
    } else {
      # Otherwise, iterate through the overlapping TADs:
      overlappingTADs <- subjectHits(overlaps)
      isTADUnique <- TRUE

      for (index in overlappingTADs) {
        # Check if overlap proportion of both TADs is less than the threshold
        otherTADstart <- start(uniqueTADs)[index]
        otherTADend <- end(uniqueTADs)[index]

        overlapStart <- max(TADstart, otherTADstart)
        overlapEnd <- min(TADend, otherTADend)
        overlapDist <- overlapEnd - overlapStart

        overlapProportion1 <- overlapDist / (TADend - TADstart)
        overlapProportion2 <- overlapDist / (otherTADend - otherTADstart)

        if (overlapProportion1 > threshold && overlapProportion2 > threshold) {
          isTADUnique <- FALSE

          # Merge TADs by taking the lower of the two start positions,
          # and the higher of the two end positions
          start(uniqueTADs)[index] <- min(TADstart, otherTADstart)
          if (TADstart < otherTADstart) {
            uniqueTADs$identifier[index] <- gsub(".*-", paste(sub("-.*", "", TADidentifier), "-", sep = ""), uniqueTADs$identifier[index])
          }

          end(uniqueTADs)[index] <- max(TADend, otherTADend)
          if (TADend > otherTADend) {
            uniqueTADs$identifier[index] <- gsub("-.*", paste("-", sub(".*-", "", TADidentifier), sep = ""), uniqueTADs$identifier[index])
          }

          break
        }
      }

      # If no overlaps greater than threshold were found, then the TAD is unique
      if (isTADUnique) {
        suppressWarnings({
          uniqueTADs <- uniqueTADs %>%
            append(TADasGR)
        })
      }
    }
  }

  return(uniqueTADs)
}

#' Save TADs to a BEDPE file
#' @param uniqueTADs List of unique TADs
#' @param output_file Path to the output file
save_tads_to_bedpe <- function(uniqueTADs, output_file) {
  # Sort uniqueTADs
  uniqueTADs <- sortSeqlevels(uniqueTADs)
  uniqueTADs <- sort(uniqueTADs)

  # Format as .bedpe file
  uniqueTADsDF <- as.data.frame(uniqueTADs)
  finalDF <- data.frame(
    uniqueTADsDF$seqnames,
    uniqueTADsDF$start,
    uniqueTADsDF$end,
    uniqueTADsDF$seqnames,
    uniqueTADsDF$start,
    uniqueTADsDF$end,
    uniqueTADsDF$identifier
  )

  # Write header
  writeLines(
    paste("#chr1", "\t", "x1", "\t", "x2", "\t", "chr2", "\t", "y1", "\t", "y2", "\t", "name"),
    output_file
  )

  # Write data
  write.table(
    finalDF,
    output_file,
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE,
    quote = FALSE,
    sep = "\t"
  )

  message(paste("Saved merged TADs to:", output_file))
}

#' Main function
main <- function() {
  # Parse and validate arguments
  args <- parse_args()
  params <- validate_args(args)

  # Convert resolutions into .bedpe file names
  resolution_files <- paste(params$resolutions, "_blocks.bedpe", sep = "")

  # Create an empty "merged list" of unique TADs across all resolutions
  uniqueTADs <- GRanges()

  # Process each resolution file
  for (resolution in resolution_files) {
    file_path <- file.path(params$inDir, resolution)

    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }

    message(paste("Processing file:", resolution))
    uniqueTADs <- process_tad_file(file_path, resolution, uniqueTADs, params$threshold)
  }

  # Save results
  output_file <- file.path(params$outDir, "mergedTADs.bedpe")
  save_tads_to_bedpe(uniqueTADs, output_file)

  message("TAD merging completed successfully")
}

# Run the main function with error handling
tryCatch(
  {
    main()
  },
  error = function(e) {
    message("Error: ", e$message)
    quit(status = 1)
  }
)
