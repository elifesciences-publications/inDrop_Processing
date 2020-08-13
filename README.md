This code is associated with the paper from Rinaldi et al., "An atlas of cell types in the mouse epididymis and vas deferens". eLife, 2020. http://doi.org/10.7554/eLife.55474


# inDrop_Processing

Initial processing and mapping indrop data

## Getting Started

You will need the demultiplexed fastqs. An example bcl2fastq command is below for this library construction. Barcodes should be specified on the SampleSheet.csv file.

```
bcl2fastq --runfolder-dir ./Files --sample-sheet ./SampleSheet.csv --output-dir ./ --use-bases-mask y58n*,y*,I*,y16n* --mask-short-adapter-reads 0 --minimum-trimmed-read-length 0 --barcode-mismatches 1
```

### Prerequisites

Nextflow;
Python;
Samtools;
Tophat2

### Run the nextflow pipeline

Edit the paths on this script to your local paths

```
./run_nextflow.sh > ./run_log.txt &
```
