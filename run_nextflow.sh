#!/bin/sh

nextflow indrop_single_cell.nf \
--cleanLowEndUmis_path "./scripts/cleanLowEndUmis.py" \
--countUniqueAlignedBarcodes_fromFile_filePath "./scripts/countUniqueAlignedBarcodes_fromFile.py" \
--ESAT_path "./scripts/esat.v0.1_09.09.16_24.18.umihack.jar" \
--filter_lowCountBC_bam_print_py_filePath "./scripts/filter_lowCountBC_bam_print.py" \
--extractValidReadsPath "./scripts/extractValidReads_V3_ps.py" \
--genome "./mm10.fa" \
--gtfFilePath "./ucsc.gtf" \
--genomeIndexPath "./mm10" \
--gene_to_transcript_mapping_file "./mm10_trans2gene_NMonly.txt" \
--cutoff_reads_for_valid_cell "100" \
--mate_split "single" \
--determined_fastq "./*_{R1,R2,R3}_001.fastq.gz" \
--cellBarcodeFile "./scripts/gel_barcode1_list.txt"
