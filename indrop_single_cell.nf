params.outdir = './indrop_merge'  

if (!params.cellBarcodeFile){params.cellBarcodeFile = ""} 
if (!params.mate_split){params.mate_split = ""} 
if (!params.genome){params.genome = ""} 
if (!params.determined_fastq){params.determined_fastq = ""} 
if (!params.genomeIndexPath){params.genomeIndexPath = ""} 
if (!params.gtfFilePath){params.gtfFilePath = ""} 
if (!params.countUniqueAlignedBarcodes_fromFile_filePath){params.countUniqueAlignedBarcodes_fromFile_filePath = ""} 
if (!params.filter_lowCountBC_bam_print_py_filePath){params.filter_lowCountBC_bam_print_py_filePath = ""} 
if (!params.cutoff_reads_for_valid_cell){params.cutoff_reads_for_valid_cell = ""} 
if (!params.ESAT_path){params.ESAT_path = ""} 
if (!params.gene_to_transcript_mapping_file){params.gene_to_transcript_mapping_file = ""} 
if (!params.cleanLowEndUmis_path){params.cleanLowEndUmis_path = ""} 
if (!params.extractValidReadsPath){params.extractValidReadsPath = ""} 

g_7_cellBarcodeFile = file(params.cellBarcodeFile) 
g_13_mate = params.mate_split
g_18_genome = file(params.genome) 
g_27_fastq_set_g_75 = Channel.fromPath(params.determined_fastq).toSortedList() 
g_34_genomeIndexPath = params.genomeIndexPath
g_36_gtfFilePath = file(params.gtfFilePath) 
g_41_script = file(params.countUniqueAlignedBarcodes_fromFile_filePath) 
g_43_script = file(params.filter_lowCountBC_bam_print_py_filePath) 
g_44_cutoff = params.cutoff_reads_for_valid_cell
g_50_script = file(params.ESAT_path) 
g_51_trans2gene_file = file(params.gene_to_transcript_mapping_file) 
g_52_script = file(params.cleanLowEndUmis_path) 
g_70_g_extractValid = file(params.extractValidReadsPath) 

// function that returns the sampleName by matching asterisk(*) of a glob pattern starting from the beginning of the name.
// eg. readPrefix('/some/data/file_alpha_S1_002_1.fa', 'file*_??_{001,002}_{1,2}.fa' ) Returns 'file_alpha'
def readPrefix( actual, template ) {
    final fileName = actual instanceof Path ? actual.name : actual.toString()
    def filePattern = template.toString()
    int p = filePattern.lastIndexOf('/')
    if( p != -1 ) filePattern = filePattern.substring(p+1)
    if( !filePattern.contains('*') && !filePattern.contains('?') ) 
        filePattern = '*' + filePattern 
    def regex = filePattern.replace('.','\\.').replace('*','(.*)').replace('?','(?:.?)').replace('{','(?:').replace('}',')').replace(',','|')
    def matcher = (fileName =~ /$regex/  )
    if( matcher.matches() ) { 
        def end = matcher.end(matcher.groupCount() )      
        def prefix = fileName.substring(0,end)
        while(prefix.endsWith('-') || prefix.endsWith('_') || prefix.endsWith('.') ) 
          prefix=prefix[0..-2]
        return prefix
    }
    return null
}

process extractValidReads {

publishDir params.outdir, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /validfastq\/.*.fastq$/) "valid_fastq/$filename"
}

input:
 file extractVcode from g_70_g_extractValid
 file cellBarcode from g_7_cellBarcodeFile
 set file(fastq1), file(fastq2), file(fastq3) from g_27_fastq_set_g_75.flatMap().buffer(size:3)

output:
 set val(sampleName), file("validfastq/*.fastq") into g_75_valid_fastq_g_73

script:
// readPrefix function returns the sampleName by matching asterisk(*) of a glob pattern starting from the beginning of the name.
// note: params.determined_fastq is pipeline input where reads are defined.
sampleName = readPrefix (fastq1, params.determined_fastq)
name = fastq1.toString() - ~/(\.fastq.gz)?(\.fq.gz)?(\.fastq)?(\.fq)?$/
"""
mkdir -p validfastq
python ${extractVcode} -i ${fastq1} -o ${name} -d validfastq -b ${cellBarcode} -u 8
"""
}


process Merge_Reads {

input:
 val mate from g_13_mate
 set val(name), file(reads) from g_75_valid_fastq_g_73.groupTuple()

output:
 set val(name), file("reads/*") into g_73_reads_g_68

script:
"""
mkdir reads
cat ${reads} > ${name}_mergeValid.fastq
mv ${name}_mergeValid.fastq reads/.
"""
}

readsPerFile  = 5000000 //* @input @description:"The number of reads per file"
//Since splitFastq operator requires flat file structure, first convert grouped structure to flat, execute splitFastq, and then return back to original grouped structure
//.map(flatPairsClosure).splitFastq(splitFastqParams).map(groupPairsClosure)

//Mapping grouped read structure to flat structure
def flatPairsClosure = {row -> (row[1] instanceof Collection) ? tuple(row[0], file(row[1][0]), file(row[1][1])) : tuple(row[0], file(row[1]))}
//Mapping flat read structure to grouped read structure
def groupPairsClosure = {row -> tuple(row[0], (row[2]) ? [file(row[1]), file(row[2])] : [file(row[1])])}

// if mate of split process different than rest of the pipeline, use "mate_split" as input parameter. Otherwise use default "mate" as input parameter
def mateParamName = (params.mate_split) ? "mate_split" : "mate"
def splitFastqParams;
if (params[mateParamName] != "pair"){
    splitFastqParams = [by: readsPerFile, file:true]
}else {
    splitFastqParams = [by: readsPerFile, pe:true, file:true]
}


process SplitFastq {

input:
 set val(name), file(reads) from g_73_reads_g_68.map(flatPairsClosure).splitFastq(splitFastqParams).map(groupPairsClosure)
 val mate from g_13_mate

output:
 set val(name), file("split/*") into g_68_split_fastq_g_72

"""
mkdir -p split
mv ${reads} split/.
"""
}

params_tophat = "" //textbox

process Map_Tophat2 {

publishDir params.outdir, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${newName}.bam$/) "mapped_reads/$filename"
}

input:
 set val(name), file(reads) from g_68_split_fastq_g_72
 file gtfPath from g_36_gtfFilePath
 file genome from g_18_genome
 val mate from g_13_mate
 val indexPath from g_34_genomeIndexPath

output:
 set val(name), file("${newName}.bam") into g_72_mapped_reads_g_56
 set val(name), file("${newName}_align_summary.txt") into g_72_summary_g_63
 set val(name), file("${newName}_unmapped.bam") into g_72_unmapped_reads

script:
nameAll = reads.toString()
nameArray = nameAll.split(' ')

if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file =  nameAll 
    runGzip = ''
}

"""
$runGzip
if [ "${mate}" == "pair" ]; then
    tophat2 -p 4 ${params_tophat}  --keep-tmp -G ${gtfPath} -o . ${indexPath} $file
else
    tophat2 -p 4 ${params_tophat}  --keep-tmp -G ${gtfPath} -o . ${indexPath} $file
fi

if [ -f unmapped.bam ]; then
    mv unmapped.bam ${newName}_unmapped.bam
else
    touch ${newName}_unmapped.bam
fi

mv accepted_hits.bam ${newName}.bam
mv align_summary.txt ${newName}_align_summary.txt
"""
}


process Merge_Bam {

publishDir params.outdir, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${oldname}.bam$/) "merged_bam/$filename"
}

input:
 set val(oldname), file(bamfiles) from g_72_mapped_reads_g_56.groupTuple()

output:
 set val(oldname), file("${oldname}.bam") into g_56_merged_bams_g_57

shell:
'''
num=$(echo "!{bamfiles.join(" ")}" | awk -F" " '{print NF-1}')
if [ "${num}" -gt 0 ]; then
samtools merge !{oldname}.bam !{bamfiles.join(" ")}
else
mv !{bamfiles.join(" ")} !{oldname}.bam
fi
'''
}


process sort_index {

input:
 set val(name), file(initial_alignment) from g_56_merged_bams_g_57

output:
 set val(name), file("*.bam") into g_57_bam_file_g_76
 set val(name), file("*.bai") into g_57_bam_index_g_76

"""
samtools sort ${initial_alignment} ${name}_sorted
samtools index ${name}_sorted.bam 
"""
}


process Count_cells {

publishDir params.outdir, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_count.txt$/) "cell_counts/$filename"
}

input:
 set val(oldname), file(sorted_bams) from g_57_bam_file_g_76
 set val(oldname), file(bams_index) from g_57_bam_index_g_76
 file script_path from g_41_script

output:
 set val(oldname), file("bam/*.bam") into g_76_sorted_bam_g_47
 set val(oldname), file("bam/*.bam.bai") into g_76_bam_index_g_47
 set val(oldname), file("*_count.txt") into g_76_count_file_g_47

"""
find  -name "*.bam" > filelist.txt
python ${script_path} -i filelist.txt -o ${oldname}_count.txt
mkdir bam
mv $sorted_bams bam/.
mv $bams_index bam/.
"""
}


process filter_lowCount {

input:
 file script_path from g_43_script
 set val(oldname), file(sorted_bams) from g_76_sorted_bam_g_47
 val cutoff from g_44_cutoff
 set val(name), file(count_file) from g_76_count_file_g_47
 set val(oldname), file(bam_index) from g_76_bam_index_g_47

output:
 set val(name), file("*_filtered.bam") into g_47_filtered_bam_g_48

"""
python ${script_path} -i ${sorted_bams} -b ${name}_count.txt -o ${name}_filtered.bam -n ${cutoff}
"""
}


process ESAT {

publishDir params.outdir, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*umi.distributions.txt$/) "UMI_distributions/$filename"
}

input:
 file esat_path from g_50_script
 file trans2gene_filepath from g_51_trans2gene_file
 set val(name), file(filtered_bam) from g_47_filtered_bam_g_48

output:
 set val(name), file("*umi.distributions.txt") into g_48_UMI_distributions_g_49

"""
find  -name "*.bam" | awk '{print "${name}\t"\$1 }' > filelist.txt
java \
-Xmx40g \
-jar ${esat_path} \
-alignments filelist.txt \
-out ${name}_esat.txt \
-geneMapping ${trans2gene_filepath} \
-task score3p \
-wLen 100 \
-wOlap 50 \
-wExt 1000 \
-sigTest .01 \
-multimap ignore \
-scPrep
"""
}


process UMI_Trim {

publishDir params.outdir, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_umiClean.txt$/) "UMI_count_final/$filename"
}

input:
 file cleanLowEndUmis_path from g_52_script
 set val(name), file(umi_dist) from g_48_UMI_distributions_g_49

output:
 set val(name), file("*_umiClean.txt") into g_49_UMI_clean

"""
python ${cleanLowEndUmis_path} \
-i ${umi_dist} \
-o ${name}_umiClean.txt \
-n 2
"""
}


process Tophat_Summary_Merge {

publishDir params.outdir, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}_tophat_sum.tsv$/) "tophat_summary/$filename"
}

input:
 set val(name), file(alignSum) from g_72_summary_g_63.groupTuple()
 val mate from g_13_mate

output:
 set val(name), file("${name}_tophat_sum.tsv") into g_63_report

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 
use Data::Dumper;

my %tsv;
my @headers = ();
my $name = "!{name}";

alteredAligned();

my @keys = keys %tsv;
my $summary = "$name"."_tophat_sum.tsv";
my $header_string = join("\\t", @headers);
`echo "$header_string" > $summary`;
foreach my $key (@keys){
	my $values = join("\\t", @{ $tsv{$key} });
	`echo "$values" >> $summary`;
}


sub alteredAligned
{
	my @files = qw(!{alignSum});
	my $multimappedSum;
	my $alignedSum;
	my $inputCountSum;
	push(@headers, "Sample");
    push(@headers, "Total Reads");
	push(@headers, "Multimapped Reads Aligned Tophat");
	push(@headers, "Unique Reads Aligned Tophat");
	foreach my $file (@files){
		my $multimapped;
		my $aligned;
		my $inputCount;
		chomp($aligned = `cat $file | grep 'Aligned pairs:' | awk '{sum=\\$3} END {print sum}'`);
		if ($aligned eq "") { # then it is single-end
		        chomp($inputCount = `cat $file | grep 'Input' | awk '{sum=\\$3} END {print sum}'`);
				chomp($aligned = `cat $file | grep 'Mapped' | awk '{sum=\\$3} END {print sum}'`);
				chomp($multimapped = `cat $file | grep 'multiple alignments' | awk '{sum+=\\$3} END {print sum}'`);
			}else{ # continue to pair end
			    chomp($inputCount = `cat $file | grep 'Input' | awk '{sum=\\$3} END {print sum}'`);
				chomp($multimapped = `cat $file | grep -A 1 'Aligned pairs:' | awk 'NR % 3 == 2 {sum+=\\$3} END {print sum}'`);
			}
        $multimappedSum += int($multimapped);
        $alignedSum += (int($aligned) - int($multimapped));
        $inputCountSum += int($inputCount);
        if ($alignedSum < 0){
            $alignedSum = 0;
        }
	}
	$tsv{$name} = [$name, $inputCountSum];
	push($tsv{$name}, $multimappedSum);
	push($tsv{$name}, $alignedSum);
}
'''

}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
