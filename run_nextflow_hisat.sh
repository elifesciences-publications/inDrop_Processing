#!/bin/sh

nextflow SingleCell-Indrop-Hisat.nf \
--cleanLowEndUmis_path "/project/umw_biocore/bin/singleCell/singleCellScripts/cleanLowEndUmis.py" \
--countUniqueAlignedBarcodes_fromFile_filePath "/project/umw_biocore/bin/singleCell/singleCellScripts/countUniqueAlignedBarcodes_fromFile.py" \
--ESAT_path "/project/umw_biocore/bin/singleCell/singleCellScripts/esat.v0.1_09.09.16_24.18.umihack.jar" \
--filter_lowCountBC_bam_print_py_filePath "/project/umw_biocore/bin/singleCell/singleCellScripts/filter_lowCountBC_bam_print.py" \
--extractValidReadsPath "/project/umw_biocore/bin/singleCell/singleCellScripts/extractValidReads_V3_ps.py" \
--cutoff_reads_for_valid_cell "3000" \
--mate_split "single" \
--cellBarcodeFile "/home/oy28w/ZSA_snRNA/Cellularbarcodes.txt" \
--determined_fastq "/home/oy28w/nextflowruns/determined_fastq_1sample/test/*_????_{R1,R2,R3}_001.fastq.gz" \
--hisat2_path "/project/umw_biocore/bin/hisat2/hisat2" \
--genomeIndexPath "/share/data/umw_biocore/genome_data/mousetest/mm10/mm10" \
--gene_to_transcript_mapping_file "/project/umw_biocore/bin/singleCell/singleCellFiles/mm10_trans2gene_NMonly.txt"
