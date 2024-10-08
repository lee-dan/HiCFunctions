# This program creates an interactive circular packing visualization for a given set of TADs. 

# Format:
# Rscript circularPackingTAD.R --inputFile ... --chromosome ... --start ... --end ... --merge ... --title ... --outputFile ...

# Args:
# --inputFile**: The .bedpe file containing the hierarchal list of TADs.
# --chromosome**: Specifies which chromosome will be visualized.
# --start**: The starting index (lower bound) of the genomic range to be visualized
# --end**: The ending index (upper bound) of the genomic range to be visualized
# --merge**: Specify the *merging* method:
# I.e. For TADs that are a subTAD of two unique, nonhierarchal TADs (that is, neither one is a child of the other), the program will choose either the bigger parent TAD or the smaller parent TAD. 
# Enter "Bigger" or "Smaller". The default is "Bigger".
# --title**: The title of the chart to be created (NO SPACES)
# --outputFile**: The name of the output file to be created (omit file ending: .html)

# Example:
# Rscript circularPackingTAD.R --inputFile TADHierarchy.bedpe --chromosome chr19 --start 43000000 --end 47000000 --merge Bigger --title TAD_Circular_Packing_Chart --outputFile TAD_Circular_Packing_Chart

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

rawData <- options.args$inputFile
chromosome <- options.args$chromosome
start <- as.double(options.args$start)
end <- as.double(options.args$end)
merge <- options.args$merge
title <- options.args$title
outFile <- options.args$outputFile

# Catch for formatting issues
if (merge != "Bigger" && merge != "Smaller") {
  merge = "Bigger"
}

# Save TADs as a data frame within selected range
tads <- read.table(file = rawData, sep="\t")
colnames(tads) <- c("chr1", "x1", "x2",
                    "chr2", "y1", "y2",
                    "identifier", "children")
tads <- tads[tads$chr1 == chromosome & tads$x1 >= start & tads$x2 <= end,]  

# As a first pass, select all TADs that are independent (i.e. they have no relatives)
master_data_edge <- data.frame(matrix(ncol = 2, nrow = 0))
tadsThatAreChildren <- strsplit(gsub('[.,;]', '', toString(tads$children)), " ")[[1]][nchar(strsplit(gsub('[.,;]', '', toString(tads$children)), " ")[[1]]) > 1]
for (tad in tads$identifier) {
  if (!(tad %in% tadsThatAreChildren)) {
    master_data_edge <- rbind(master_data_edge, c(" ", tad))
  }
}
colnames(master_data_edge) <- c("from", "to")

# Iterate through the remaining TADs 
data_edge <- data.frame(matrix(ncol = 2, nrow = 0))

#---Function to retrieve all children TADs of a given TAD---
getChildren <- function(tad) {
  index <- which(tads$identifier == tad)
  if (tads$children[index] == ".") {
    return (c())
  }
  return(strsplit(tads$children[index], "; ")[[1]])
}

# For each TAD, iterate through its children, and add edges that correspond
# only to a TAD's DIRECT descendants
for (i in 1:nrow(tads)) {
  children <- getChildren(tads$identifier[i])
  
  # To do so, iterate through all of its children TADS, and add only that tads
  # that are not children of any of the other children TADs
  subchildren <- c()
  if (length(children) != 0) {
    for (child in children) {
      subchildren <- append(subchildren, getChildren(child))
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

# Iterate through list of edges, removing those with duplicate (or more) parents
# that are not hierarchical, choosing either the largest or the smallest parent TAD
uniqueOnes <- unique(data_edge$to)
if (length(uniqueOnes) != length(data_edge$to)) {
  notUnique <- unique(data_edge$to[duplicated(data_edge$to) == TRUE])
  
  for (tad in notUnique) {
    rownums <- which(data_edge$to %in% tad)
    
    minSize = .Machine$double.xmax
    maxSize = -.Machine$double.xmax
    chosenParent = 0
    for (row in rownums) {
      currentParent = data_edge$from[row]
      if (merge == "Bigger") {
        if (tads[tads$identifier == currentParent, "x2"] - tads[tads$identifier == currentParent, "x1"] > maxSize) {
          maxSize = tads[tads$identifier == currentParent, "x2"] - tads[tads$identifier == currentParent, "x1"]
          chosenParent = row
        }
      } else if (merge == "Smaller") {
        if (tads[tads$identifier == currentParent, "x2"] - tads[tads$identifier == currentParent, "x1"] < minSize) {
          minSize = tads[tads$identifier == currentParent, "x2"] - tads[tads$identifier == currentParent, "x1"]
          chosenParent = row
        }
      }
    }
    rownums <- rownums[! rownums %in% c(chosenParent)]
    data_edge <- data_edge[-rownums,]
  }
}

# Adjust row names after removing edges
row.names(data_edge) <- 1:nrow(data_edge)

# Add master TAD to edges list
master <- data.frame(c(NA), c(" "))
colnames(master) <- c("from", "to")
master <- rbind(master, data_edge)
data_edge <- master

# Create TAD class
setClass("TAD", slots = list(name = "character",
                             id = "numeric",
                             chr = "character",
                             parent = "numeric",
                             value = "numeric",
                             start = "numeric",
                             end = "numeric"))
  
# Create List of Tads
TADList <- c()
for (i in 1:nrow(data_edge)) {
  if (data_edge$to[i] == " ") {
    value1 = tads[nrow(tads), "x2"] - tads[1, "x1"]
    start1 = tads[1, "x1"]
    end1 = tads[nrow(tads), "x2"]
  } else {
    value1 = tads[tads$identifier == data_edge$to[i], "x2"] - tads[tads$identifier == data_edge$to[i], "x1"]
    start1 = tads[tads$identifier == data_edge$to[i], "x1"]
    end1 = tads[tads$identifier == data_edge$to[i], "x2"]
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

# Create JSON file representation
json <- "["

index <- 1
for (tad in TADList) {
  jsonConvert <- "\n\t{"
  jsonConvert <- paste(jsonConvert, "\n\t\t\"name\": \"", tad@name, "\",", sep="")
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
  index = index + 1
  if (index <= length(TADList)) {
    json <- paste(json, ",", sep = "")
  }
}
json <- paste(json, "\n]", sep = "")

# Create an HTML file using circular packing functions from cdn.anychart.com 
jsonFile <- paste(
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
\t\t\t\t\t\t\t\'<span style = \"color: #112B3C;font-weight:600;font-size:18px;\">",
title,
"</span>\'
\t\t\t\t\t\t);

\t\t\t\t\t// customize the appearance
\t\t\t\t\t// Adjust thickness to 1 or 0.25

\t\t\t\t\tchart.background(\'#ffffff\');
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
\t\t\t\t\t\t.fontSize(\'14\')
\t\t\t\t\t\t.fontColor(\'#696969\')
\t\t\t\t\t\t.textShadow(\'none\')
\t\t\t\t\t\t.anchor(\'center-top\').offsetY(\'-3%\');

\t\t\t\t\tchart.labels()
\t\t\t\t\t\t.background()
\t\t\t\t\t\t.enabled(true)
\t\t\t\t\t\t.fill(\"#f6f6f6 0.8\")
\t\t\t\t\t\t.stroke(\"#888888\")
\t\t\t\t\t\t.corners(5);	

\t\t\t\t\t// specify the container element id
\t\t\t\t\tchart.container(\'container\');

\t\t\t\t\t// initiate the drawing of the chart
\t\t\t\t\tchart.draw();

\t\t\t\t\t}
\t\t\t\t);
\t</script>
</html>
",
sep = "")

# Write html file
fileConn <- file(paste(outFile, ".html", sep=""))
writeLines(jsonFile, fileConn)
close(fileConn)