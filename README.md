# HiCFunctions

[![R](https://img.shields.io/badge/R-3.6.0+-blue.svg)](https://www.r-project.org/)
[![GenomicRanges](https://img.shields.io/badge/GenomicRanges-1.40.0+-green.svg)](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)

A comprehensive set of R scripts for analyzing and visualizing Topologically Associating Domains (TADs) from Hi-C data.

## Overview

HiCFunctions provides a suite of tools for processing, merging, and visualizing TADs at different resolutions. These tools are designed to work with standard .bedpe files and produce hierarchical representations of TADs.

## Installation

```bash
# Clone the repository
git clone https://github.com/daniellee/HiCFunctions.git

# Navigate to the repository
cd HiCFunctions

# Install required R packages
R -e "install.packages(c('GenomicRanges', 'magrittr', 'dplyr', 'progress', 'data.table'), repos='https://cran.rstudio.com/')"
```

## Usage

### 1. TAD Merger [TADMerger.R]

This program merges TADs at different resolutions based on a maximum allowed threshold of shared area.

#### Command Format:

```bash
Rscript TADMerger.R --inputDirectory <dir> --outputDirectory <dir> --resolutions <values> --threshold <value>
```

#### Parameters:

- **--inputDirectory**: Directory containing the .bedpe files of TADs
- **--outputDirectory**: Directory where the merged TADs will be saved
- **--resolutions**: Resolutions of TADs to be merged (e.g., 10000, 25000, 50000)
- **--threshold**: If two TADs share more than [threshold] of their 'linear area' (distance from start index to end index), then the TADs will be 'merged'
  - Example: The boundaries of TAD A are [0, 100,000], and the boundaries of TAD B are [25,000, 110,000]. If [threshold] = 0.7, then TAD A and TAD B will be merged, forming a single TAD with boundaries [0, 110,000]. Their shared area is [100,000] - [25,000] = [75,000]. [75,000] > 0.7 x [100,000] (the area of TAD A), and [75,000] > 0.7 x [85,000] (the area of TAD B).
  - **Recommended value: 0.7**

#### Example:

```bash
Rscript TADMerger.R --inputDirectory inDir --resolutions 10000 25000 --threshold 0.7 --outputDirectory outDir
```

#### Output:

The program generates a .bedpe file named 'mergedTADs.bedpe' in the specified output directory. The last column, 'name', is an identifier for each TAD, indicating from which original TAD the lower boundary was obtained from, and which original TAD the upper boundary was obtained from.

### 2. TAD Hierarchy [TADHierarchy.R]

This program generates a hierarchy list for each TAD, detailing each of its 'children' TADs.

#### Command Format:

```bash
Rscript TADHierarchy.R --inputFile <file> --outputDirectory <dir> --threshold <value>
```

#### Parameters:

- **--inputFile**: The .bedpe file containing the list of TADs
- **--outputDirectory**: Directory where the TAD hierarchy will be saved
- **--threshold**: The minimum overlap percentage [0 to 1] for a TAD to be considered a 'child' TAD of another TAD
  - **Recommended value: 1.0**

#### Example:

```bash
Rscript TADHierarchy.R --inputFile sampleTADs.bedpe --outputDirectory outDir --threshold 1.0
```

#### Output:

The program generates a .bedpe file named 'TADHierarchy.bedpe' in the specified output directory.

### 3. Circular Packing Visualization [circularPackingTAD.R]

This program creates an interactive circular packing visualization for a given set of TADs. It requires a TAD Hierarchy file (created from TAD Hierarchy).

#### Command Format:

```bash
Rscript circularPackingTAD.R --inputFile <file> --chromosome <chr> --start <pos> --end <pos> --merge <method> --title <text> --outputFile <name>
```

#### Parameters:

- **--inputFile**: The .bedpe file containing the hierarchical list of TADs
- **--chromosome**: Specifies which chromosome will be visualized
- **--start**: The starting index (lower bound) of the genomic range to be visualized
- **--end**: The ending index (upper bound) of the genomic range to be visualized
- **--merge**: Specify the merging method:
  - For TADs that are a subTAD of two unique, nonhierarchical TADs (that is, neither one is a child of the other), the program will choose either the bigger parent TAD or the smaller parent TAD
  - Options: "Bigger" or "Smaller" (default: "Bigger")
- **--title**: The title of the chart to be created (NO SPACES)
- **--outputFile**: The name of the output file to be created (omit file ending: .html)

#### Example:

```bash
Rscript circularPackingTAD.R --inputFile TADHierarchy.bedpe --chromosome chr19 --start 43000000 --end 47000000 --merge Bigger --title TAD_Circular_Packing_Chart --outputFile TAD_Circular_Packing_Chart
```

#### Output:

The program generates an HTML file with the specified name that can be opened in a web browser to view the interactive visualization.

## Sample Files

The repository includes sample files in the "Examples" directory:

- Sample input files from 4DN
- Merged TADs examples
- TAD Hierarchy examples
- Circular Packing visualization examples

## Dependencies

- R (>= 3.6.0)
- GenomicRanges
- magrittr
- dplyr
- progress
- data.table

## Author

**Daniel Lee** - [GitHub](https://github.com/lee-dan)

## Acknowledgments

This project was developed for analyzing and visualizing Topologically Associating Domains (TADs) from Hi-C data. Special thanks to the 4DN Consortium for providing sample data.
