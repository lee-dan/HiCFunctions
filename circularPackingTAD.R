#!/usr/bin/env Rscript
# circularPackingTAD.R
# This program creates an interactive circular packing visualization for a given set of TADs.
# Author: Daniel Lee
# Date: April 2, 2025
# Version: 1.0.0

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
    cat("Usage: Rscript circularPackingTAD.R --inputFile <file> --chromosome <chr> --start <pos> --end <pos> --merge <method> --title <text> --outputFile <name>\n")
    cat("\nArguments:\n")
    cat("  --inputFile: The .bedpe file containing the hierarchical list of TADs\n")
    cat("  --chromosome: Specifies which chromosome will be visualized\n")
    cat("  --start: The starting index (lower bound) of the genomic range to be visualized\n")
    cat("  --end: The ending index (upper bound) of the genomic range to be visualized\n")
    cat("  --merge: Specify the merging method: 'Bigger' or 'Smaller' (default: 'Bigger')\n")
    cat("  --title: The title of the chart to be created (NO SPACES)\n")
    cat("  --outputFile: The name of the output file to be created (omit file ending: .html)\n")
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
  required_args <- c("inputFile", "chromosome", "start", "end", "title", "outputFile")
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

  # Validate chromosome
  if (is.null(args$chromosome) || args$chromosome == "") {
    stop("Chromosome must be specified")
  }

  # Validate start and end positions
  start <- suppressWarnings(as.double(args$start))
  end <- suppressWarnings(as.double(args$end))

  if (is.na(start) || is.na(end)) {
    stop("Start and end positions must be numeric values")
  }

  if (start >= end) {
    stop("Start position must be less than end position")
  }

  # Validate merge method
  merge <- args$merge
  if (is.null(merge) || merge == "") {
    merge <- "Bigger"
  } else if (merge != "Bigger" && merge != "Smaller") {
    warning("Invalid merge method. Using default: 'Bigger'")
    merge <- "Bigger"
  }

  # Validate title
  if (is.null(args$title) || args$title == "") {
    stop("Title must be specified")
  }

  # Validate output file
  if (is.null(args$outputFile) || args$outputFile == "") {
    stop("Output file name must be specified")
  }

  return(list(
    rawData = args$inputFile,
    chromosome = args$chromosome,
    start = start,
    end = end,
    merge = merge,
    title = args$title,
    outFile = args$outputFile
  ))
}

#' Create TAD class
create_tad_class <- function() {
  setClass("TAD", slots = list(
    name = "character",
    id = "numeric",
    chr = "character",
    parent = "numeric",
    value = "numeric",
    start = "numeric",
    end = "numeric"
  ))
}

#' Read and filter TADs from input file
#' @param rawData Path to the input file
#' @param chromosome Chromosome to filter
#' @param start Start position
#' @param end End position
#' @return Data frame with filtered TADs
read_and_filter_tads <- function(rawData, chromosome, start, end) {
  tryCatch(
    {
      tads <- read.table(file = rawData, sep = "\t")
      colnames(tads) <- c(
        "chr1", "x1", "x2",
        "chr2", "y1", "y2",
        "identifier", "children"
      )

      # Filter TADs by chromosome and position
      filtered_tads <- tads[tads$chr1 == chromosome & tads$x1 >= start & tads$x2 <= end, ]

      if (nrow(filtered_tads) == 0) {
        stop(paste("No TADs found in the specified range:", chromosome, ":", start, "-", end))
      }

      return(filtered_tads)
    },
    error = function(e) {
      stop(paste("Error reading or filtering TADs:", e$message))
    }
  )
}

#' Get children TADs of a given TAD
#' @param tad TAD identifier
#' @param tads Data frame with TAD information
#' @return Vector of child TAD identifiers
get_children <- function(tad, tads) {
  index <- which(tads$identifier == tad)
  if (tads$children[index] == ".") {
    return(c())
  }
  return(strsplit(tads$children[index], "; ")[[1]])
}

#' Create edge data frame for TAD hierarchy
#' @param tads Data frame with TAD information
#' @return Data frame with edge information
create_edge_data <- function(tads) {
  # As a first pass, select all TADs that are independent (i.e. they have no relatives)
  master_data_edge <- data.frame(matrix(ncol = 2, nrow = 0))
  tadsThatAreChildren <- strsplit(gsub("[.,;]", "", toString(tads$children)), " ")[[1]][nchar(strsplit(gsub("[.,;]", "", toString(tads$children)), " ")[[1]]) > 1]

  for (tad in tads$identifier) {
    if (!(tad %in% tadsThatAreChildren)) {
      master_data_edge <- rbind(master_data_edge, c(" ", tad))
    }
  }
  colnames(master_data_edge) <- c("from", "to")

  # Iterate through the remaining TADs
  data_edge <- data.frame(matrix(ncol = 2, nrow = 0))

  # For each TAD, iterate through its children, and add edges that correspond
  # only to a TAD's DIRECT descendants
  for (i in 1:nrow(tads)) {
    children <- get_children(tads$identifier[i], tads)

    # To do so, iterate through all of its children TADS, and add only that tads
    # that are not children of any of the other children TADs
    subchildren <- c()
    if (length(children) != 0) {
      for (child in children) {
        subchildren <- append(subchildren, get_children(child, tads))
      }
    }

    for (child in children) {
      if (!(child %in% subchildren)) {
        data_edge <- rbind(data_edge, c(tads$identifier[i], child))
      }
    }
  }
  colnames(data_edge) <- c("from", "to")

  # Combine list of all edges thus far
  master_data_edge <- rbind(master_data_edge, data_edge)
  data_edge <- master_data_edge

  return(data_edge)
}

#' Process edges to handle duplicate parents
#' @param data_edge Data frame with edge information
#' @param tads Data frame with TAD information
#' @param merge Merge method ("Bigger" or "Smaller")
#' @return Processed edge data frame
process_edges <- function(data_edge, tads, merge) {
  # Iterate through list of edges, removing those with duplicate (or more) parents
  # that are not hierarchical, choosing either the largest or the smallest parent TAD
  uniqueOnes <- unique(data_edge$to)
  if (length(uniqueOnes) != length(data_edge$to)) {
    notUnique <- unique(data_edge$to[duplicated(data_edge$to) == TRUE])

    for (tad in notUnique) {
      rownums <- which(data_edge$to %in% tad)

      minSize <- .Machine$double.xmax
      maxSize <- -.Machine$double.xmax
      chosenParent <- 0

      for (row in rownums) {
        currentParent <- data_edge$from[row]
        if (merge == "Bigger") {
          if (tads[tads$identifier == currentParent, "x2"] - tads[tads$identifier == currentParent, "x1"] > maxSize) {
            maxSize <- tads[tads$identifier == currentParent, "x2"] - tads[tads$identifier == currentParent, "x1"]
            chosenParent <- row
          }
        } else if (merge == "Smaller") {
          if (tads[tads$identifier == currentParent, "x2"] - tads[tads$identifier == currentParent, "x1"] < minSize) {
            minSize <- tads[tads$identifier == currentParent, "x2"] - tads[tads$identifier == currentParent, "x1"]
            chosenParent <- row
          }
        }
      }

      rownums <- rownums[!rownums %in% c(chosenParent)]
      data_edge <- data_edge[-rownums, ]
    }
  }

  # Adjust row names after removing edges
  row.names(data_edge) <- 1:nrow(data_edge)

  # Add master TAD to edges list
  master <- data.frame(c(NA), c(" "))
  colnames(master) <- c("from", "to")
  master <- rbind(master, data_edge)
  data_edge <- master

  return(data_edge)
}

#' Create TAD objects from edge data
#' @param data_edge Data frame with edge information
#' @param tads Data frame with TAD information
#' @param chromosome Chromosome
#' @return List of TAD objects
create_tad_objects <- function(data_edge, tads, chromosome) {
  TADList <- c()

  for (i in 1:nrow(data_edge)) {
    if (data_edge$to[i] == " ") {
      value1 <- tads[nrow(tads), "x2"] - tads[1, "x1"]
      start1 <- tads[1, "x1"]
      end1 <- tads[nrow(tads), "x2"]
    } else {
      value1 <- tads[tads$identifier == data_edge$to[i], "x2"] - tads[tads$identifier == data_edge$to[i], "x1"]
      start1 <- tads[tads$identifier == data_edge$to[i], "x1"]
      end1 <- tads[tads$identifier == data_edge$to[i], "x2"]
    }

    TADobj <- new("TAD",
      name = data_edge$to[i],
      id = i,
      chr = chromosome,
      parent = match(data_edge$from[i], data_edge$to),
      value = value1,
      start = start1,
      end = end1
    )
    TADList <- append(TADList, TADobj)
  }

  return(TADList)
}

#' Convert TAD objects to JSON
#' @param TADList List of TAD objects
#' @return JSON string
convert_to_json <- function(TADList) {
  json <- "["

  index <- 1
  for (tad in TADList) {
    jsonConvert <- "\n\t{"
    jsonConvert <- paste(jsonConvert, "\n\t\t\"name\": \"", tad@name, "\",", sep = "")
    addParent <- FALSE
    if (!is.na(tad@parent)) {
      addParent <- TRUE
    }

    if (!addParent) { # This is the last
      jsonConvert <- paste(jsonConvert, "\n\t\t\"id\": ", tad@id, ",", sep = "")
      jsonConvert <- paste(jsonConvert, "\n\t\t\"chromosome\": \"", tad@chr, "\",", sep = "")
      jsonConvert <- paste(jsonConvert, "\n\t\t\"value\": ", tad@value, ",", sep = "")
      jsonConvert <- paste(jsonConvert, "\n\t\t\"start\": ", tad@start, ",", sep = "")
      jsonConvert <- paste(jsonConvert, "\n\t\t\"end\": ", tad@end, sep = "")
      jsonConvert <- paste(jsonConvert, "\n\t}", sep = "")
    } else { # Add parent
      jsonConvert <- paste(jsonConvert, "\n\t\t\"id\": ", tad@id, ",", sep = "")
      jsonConvert <- paste(jsonConvert, "\n\t\t\"chromosome\": \"", tad@chr, "\",", sep = "")
      jsonConvert <- paste(jsonConvert, "\n\t\t\"parent\": ", tad@parent, ",", sep = "")
      jsonConvert <- paste(jsonConvert, "\n\t\t\"value\": ", tad@value, ",", sep = "")
      jsonConvert <- paste(jsonConvert, "\n\t\t\"start\": ", tad@start, ",", sep = "")
      jsonConvert <- paste(jsonConvert, "\n\t\t\"end\": ", tad@end, sep = "")
      jsonConvert <- paste(jsonConvert, "\n\t}", sep = "")
    }

    json <- paste(json, jsonConvert, sep = "")
    index <- index + 1
    if (index <= length(TADList)) {
      json <- paste(json, ",", sep = "")
    }
  }

  json <- paste(json, "\n]", sep = "")
  return(json)
}

#' Create HTML file with circular packing visualization
#' @param json JSON data
#' @param title Chart title
#' @return HTML content
create_html_content <- function(json, title) {
  html_content <- paste(
    "
<!DOCTYPE html>
<html lang=\"en\">
\t<head>
\t\t<meta charset=\"utf-8\">
\t\t<title>Circle Packing Chart in JavaScript</title>
\t\t<script src=\"https://cdn.anychart.com/releases/8.11.0/js/anychart-core.min.js\"></script>
\t\t<script src=\"https://cdn.anychart.com/releases/8.11.0/js/anychart-circle-packing.min.js\"></script>
\t\t<script src=\"https://cdn.anychart.com/releases/8.11.0/js/anychart-data-adapter.min.js\"></script>
\t\t<style type=\"text/css\">
\t\t\thtml,
\t\t\tbody,
\t\t\t#container {
\t\t\t\twidth: 100%;
\t\t\t\theight: 100%;
\t\t\t\tmargin: 0;
\t\t\t\tpadding: 0;
\t\t\t}
\t\t</style>
\t</head>
\t<body>
\t\t<div id=\"container\"></div>
\t\t\t<script>

\t\t\t\tanychart.onDocumentReady(function () {

\t\t\t\t\t// load a json data file
\t\t\t\t\tvar data =
",
    json,
    ";
\t\t\t\t\t// add the data
\t\t\t\t\tvar treeData = anychart.data.tree(data, 'as-table');

\t\t\t\t\t// create a circle packing chart instance
\t\t\t\t\tvar chart = anychart.circlePacking(treeData);

\t\t\t\t\t//customize the tooltip
\t\t\t\t\t//.toExponential()
\t\t\t\t\tchart
\t\t\t\t\t\t.tooltip()
\t\t\t\t\t\t.useHtml(true)
\t\t\t\t\t\t.format(function () {
\t\t\t\t\t\t\treturn '<div>'
\t\t\t\t\t\t\t\t\t+ '<span>Chromosome: ' + this.item.get('chromosome') + '</span><br/>'
\t\t\t\t\t\t\t\t\t+ '<span>Size: ' + this.value + '</span><br/>'
\t\t\t\t\t\t\t\t\t+ '<span>Start: ' + this.item.get('start') + '</span><br/>'
\t\t\t\t\t\t\t\t\t+ '<span>End: ' + this.item.get('end') + '</span>'
\t\t\t\t\t\t\t\t\t+ '</div>'
\t\t\t\t\t\t});

\t\t\t\t\t// add a chart title
\t\t\t\t\tchart
\t\t\t\t\t\t.title()
\t\t\t\t\t\t.enabled(true)
\t\t\t\t\t\t.useHtml(true)
\t\t\t\t\t\t.text(
\t\t\t\t\t\t\t'<span style = \"color: #112B3C;font-weight:600;font-size:18px;\">",
    title,
    "</span>\'
\t\t\t\t\t\t);

\t\t\t\t\t// customize the appearance
\t\t\t\t\t// Adjust thickness to 1 or 0.25

\t\t\t\t\tchart.background('#ffffff');
\t\t\t\t\tchart
\t\t\t\t\t\t.hovered()
\t\t\t\t\t\t.stroke(function () {
\t\t\t\t\t\t\treturn {
\t\t\t\t\t\t\t\tthickness: 0.7,
\t\t\t\t\t\t\t};
\t\t\t\t\t\t});
\t\t\t\t\tchart
\t\t\t\t\t\t.stroke(function () {
\t\t\t\t\t\t\treturn {
\t\t\t\t\t\t\t\tthickness: 0.3,
\t\t\t\t\t\t\t};
\t\t\t\t\t\t});

\t\t\t\t\t// customize the labels
\t\t\t\t\tchart
\t\t\t\t\t\t.labels()
\t\t\t\t\t\t.fontSize('14')
\t\t\t\t\t\t.fontColor('#696969')
\t\t\t\t\t\t.textShadow('none')
\t\t\t\t\t\t.anchor('center-top').offsetY('-3%');

\t\t\t\t\tchart.labels()
\t\t\t\t\t\t.background()
\t\t\t\t\t\t.enabled(true)
\t\t\t\t\t\t.fill(\"#f6f6f6 0.8\")
\t\t\t\t\t\t.stroke(\"#888888\")
\t\t\t\t\t\t.corners(5);

\t\t\t\t\t// specify the container element id
\t\t\t\t\tchart.container('container');

\t\t\t\t\t// initiate the drawing of the chart
\t\t\t\t\tchart.draw();

\t\t\t\t\t}
\t\t\t\t);
\t</script>
</html>
",
    sep = ""
  )

  return(html_content)
}

#' Save HTML file
#' @param html_content HTML content
#' @param output_file Path to the output file
save_html_file <- function(html_content, output_file) {
  tryCatch(
    {
      fileConn <- file(paste(output_file, ".html", sep = ""))
      writeLines(html_content, fileConn)
      close(fileConn)
      message(paste("Saved HTML file to:", paste(output_file, ".html", sep = "")))
    },
    error = function(e) {
      stop(paste("Error saving HTML file:", e$message))
    }
  )
}

#' Main function
main <- function() {
  # Parse and validate arguments
  args <- parse_args()
  params <- validate_args(args)

  # Create TAD class
  create_tad_class()

  # Read and filter TADs
  message("Reading and filtering TADs...")
  tads <- read_and_filter_tads(params$rawData, params$chromosome, params$start, params$end)

  # Create edge data
  message("Creating TAD hierarchy...")
  data_edge <- create_edge_data(tads)

  # Process edges
  data_edge <- process_edges(data_edge, tads, params$merge)

  # Create TAD objects
  message("Creating TAD objects...")
  TADList <- create_tad_objects(data_edge, tads, params$chromosome)

  # Convert to JSON
  message("Converting to JSON...")
  json <- convert_to_json(TADList)

  # Create HTML content
  message("Creating HTML content...")
  html_content <- create_html_content(json, params$title)

  # Save HTML file
  save_html_file(html_content, params$outFile)

  message("Circular packing visualization completed successfully")
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
