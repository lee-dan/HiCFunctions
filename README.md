# HiCFunctions
## Set of programs to operate on TADs

## Availability
The included programs are available via **command line** on this GitHub repository.

## Usage

### 1) TAD Merger [TADMerger.R]
This program merges TADs at different resolutions based on a maximum allowed threshold of shared area.

#### Format:
```
Rscript TADMerger.R --inputDirectory ... --outputDirectory ... --resolutions ... --threshold ...
```
* **--inputDirectory**: The directory which contains the .bedpe files of TADs
* **--outputDirectory**: The directory in which the merged TADs will be placed in
* **--resolutions**: The resolutions of TADs to be merged [e.g. 10000, 25000, 50000, etc.]
* **--threshold**: If two TADs share more than [threshold] of their 'linear area' (distance from start index to end index), then the TADs will be 'merged'. 
    * Eg. The boundaries of TAD A are [0, 100,000], and the boundaries of TAD B are [25,000, 110,000]. If [threshold] = 0.7, then TAD A and TAD B will be merged, forming a single TAD with boundaries [0, 110,000]. Their shared area is [100,000] - [25,000] = [75,000]. [75,000] > 0.7 x [100,000] (the area of TAD A), and [75,000] > 0.7 x [85,000] (the area of TAD B).
    *  **The recommended value is 0.7.**

#### Example:
```
Rscript TADMerger.R --inputDirectory inDir --resolutions 10000 25000 --threshold 0.7 --outputDirectory outDir
```
#### Output:
The program will generate a .bedpe file named 'mergedTADs.bedpe' in the _outputDirectory_. The last column entitled, 'name', is an identifier for each TAD, indicating from which original TAD the lower boundary was obtained from, and which original TAD the upper boundary was obtained from.

### 2) TAD Hierarchy [TADHierarchy.R]
This program generates a hierarchy list for each TAD, detailing each of its 'children' TADs.

#### Format:
```
Rscript TADHierarchy.R --inputFile ... --outputDirectory ... --threshold ...
```
* **--inputFile**: The .bedpe file containing the list of TADs.
* **--outputDirectory**: The directory in which the TAD hierarchy will be placed in.
* **--threshold**: The minimum overlap percentage [0 to 1] for a TAD to be considered a 'child' TAD of another TAD.
    *  **The recommended value is 1.0**. 

#### Example:
```
Rscript TADHierarchy.R --inputFile sampleTADs.bedpe --outputDirectory outDir --threshold 1.0
```
#### Output:
The program will generate a .bedpe file named 'TADHierarchy.bedpe' in the _outputDirectory_. 

### 3) Circular Packing Visualization [circularPackingTAD.R]
This program creates an interactive circular packing visualization for a given set of TADs. 
It requires a TAD Hierarchy file (created from *2) TAD Hierarchy*)

#### Format:
```
Rscript circularPackingTAD.R --inputFile ... --chromosome ... --start ... --end ... --merge ... --title ... --outputFile ...
```

* **--inputFile**: The .bedpe file containing the hierarchal list of TADs.
* **--chromosome**: Specifies which chromosome will be visualized.
* **--start**: The starting index (lower bound) of the genomic range to be visualized
* **--end**: The ending index (upper bound) of the genomic range to be visualized
* **--merge**: Specify the *merging* method:
    * I.e. For TADs that are a subTAD of two unique, nonhierarchal TADs (that is, neither one is a child of the other), the program will choose either the bigger parent TAD or the smaller parent TAD. 
    * Enter "Bigger" or "Smaller". The default is "Bigger".
* **--title**: The title of the chart to be created (NO SPACES)
* **--outputFile**: The name of the output file to be created (omit file ending: .html)

#### Example:
```
Rscript circularPackingTAD.R --inputFile TADHierarchy.bedpe --chromosome chr19 --start 43000000 --end 47000000 --merge Bigger --title TAD_Circular_Packing_Chart --outputFile TAD_Circular_Packing_Chart
```
#### Output:
The program will generate a .html file named '*outputFile*'.

## Sample Files:
See the folder, "Examples", for a sample input file from 4DN, along with the "Merged TADs", "TAD Hierarchy", and "Circular Packing".
