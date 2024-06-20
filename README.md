# HiCFunctions
## Set of programs to operate on TADs and Loops

## Availability
The included programs are available via **command line** on this GitHub repository.

## Usage

### 1) TAD Merger [TADMerger.R]
This program merges TADs at different resolutions based on a maximum allowed threshold of shared area.

#### Format:
```
Rscript TADMerger.R --inputDirectory ... --outputDirectory ... --resolutions ... --threshold ...
```
*  **--inputDirectory**: The directory which contains the .bedpe files of TADs
* **--outputDirectory**: The directory in which the merged TADs will be placed in
*  **--resolutions**: The resolutions of TADs to be merged [e.g. 10000, 25000, 50000, etc.]
*  **--threshold**: The maximum percentage [0 to 1] of linear overlapped area between any two TADs in order for them to still be considered unique. That is, if the area shared by TAD A and TAD B is greater than or equal to [threshold] of the area of both TADs, then the two TADs will be merged. 
    * Eg. The boundaries of TAD A are [0, 100,000], and the boundaries of TAD B are [20,000, 90,000]. If the [threshold] = 0.7, then TAD A and TAD B will be merged.
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
*  **--inputFile**: The .bedpe file containing the list of TADs.
* **--outputDirectory**: The directory in which the TAD hierarchy will be placed in.
* **--threshold**: The minimum overlap percentage [0 to 1] for a TAD to be considered a 'child' TAD of another TAD.
    *  **The recommended value is 1.0**. 

#### Example:
```
Rscript TADHierarchy.R --inputFile sampleTADs.bedpe --outputDirectory outDir --threshold 1.0
```
#### Output:
The program will generate a .bedpe file named 'TADHierarchy.bedpe' in the _outputDirectory_. 

## Sample Files:
See the folder, "Examples", for a sample input file from 4DN, along with the "Merged TADs" and "TAD Hierarchy".