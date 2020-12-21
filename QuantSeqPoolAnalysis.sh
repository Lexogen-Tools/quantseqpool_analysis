#!/usr/bin/env bash

# this script assumes that the untrimmed, unprocessed fastq file are available as ./ ${sample} / R1.fastq.gz in a "fastq" sub-directory in the baseDir directory
# in detail, the result of the bash command "ls <baseDir>/fastq" should contain the name of the sample that is referenced with -s. 

fastq_raw=fastq
fastq_extracted=extracted
fastq_trimmed=trimmed
alignment_deduplicated=deduplicated

# UMI is fixed at N10

usage() { echo "Usage: $0 [-s|--sample <sampleName>] [-g|--gtfFile <gtfFile>] [-d|--starGenomeDir <dir>] (Optional: [-b|--baseDir <dir>] [-f|--feature <counting_feature>] [-a|--attribute <grouping attribute>] [-t|--threads] <# of threads>)" 1>&2; exit 1; }

my_needed_commands="umi_tools samtools cutadapt STAR gawk featureCounts"

missing_counter=0
for needed_command in $my_needed_commands; do
  if ! hash $needed_command >/dev/null 2>&1; then
    printf "Command not found in PATH or environment: %s\n" $needed_command >&2
    ((missing_counter++))
  fi
done

if ((missing_counter > 0)); then
  printf "Minimum %d commands are missing in PATH or environment, aborting\n" $missing_counter >&2
  exit 1
fi

baseDir=$(pwd)"/"
counting_feature="exon"
identifying_attribute="gene_id"
nrThreads=1

OPT=$(getopt -o s:g:d:u:b:f:a:t: --long sampleName:,gtfFile:,starGenomeDir:,baseDir:,feature:,attribute:,threads:, -- "$@")

eval set -- "$OPT"

while true; do
  case "$1" in
    -s | --sampleName ) sample="$2"; shift 2;;
    -g | --gtfFile ) gtfFile="$2"; shift 2;;
    -d | --starGenomeDir ) genomeDir="$2"; shift 2;;
    -b | --baseDir ) baseDir="$2"; shift 2;;
    -f | --feature ) counting_feature="$2"; shift 2;;
    -a | --attribute ) identifying_attribute="$2"; shift 2;;
    -t | --threads ) nrThreads="$2"; shift 2;;
    -- ) shift; break;;
    * ) usage;;
  esac
done

[ $# -gt 0 ] && { echo "no command line arguments outside options." 1>&2; usage; }

[ -z ${sample+x} ] && { echo "sample name has to be specified." 1>&2; usage; }

[ -z ${gtfFile+x} ] && { echo "gtfFile has to be specified." 1>&2; usage; }

[ -z ${genomeDir+x} ] && { echo "starGenomeDir has to be specified." 1>&2; usage; }

genomeDir=$(readlink -f ${genomeDir})
gtfFile=$(readlink -f ${gtfFile})
baseDir=$(readlink -f ${baseDir})

pushd $baseDir

# UMI extract
echo "extracting N10 UMIs from read 2"
mkdir -p ${fastq_extracted}/${sample}/
umi_tools extract --extract-method=string --bc-pattern X --bc-pattern2 NNNNNNNNNN -L ${fastq_extracted}/${sample}/extraction_log -S ${fastq_extracted}/${sample}/R1.fastq.gz -I ${fastq_raw}/${sample}/R1.fastq.gz --read2-in=${fastq_raw}/${sample}/R2.fastq.gz --read2-out=/dev/null

# trimming
echo "trimming"
mkdir -p ${fastq_trimmed}/${sample}/
cutadapt -m 20 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 ${fastq_extracted}/${sample}/R1.fastq.gz | cutadapt  -m 20 -O 3 --nextseq-trim=10  -a "r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - | cutadapt -m 20 -O 3 -a "r1polyA=A{18}" - | cutadapt -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed -o ${fastq_trimmed}/${sample}/R1.fastq.gz - 

# alignment
echo "alignment"
mkdir -p alignment/${sample}/
STAR --runThreadN ${nrThreads} --readFilesCommand zcat --genomeDir ${genomeDir} --readFilesIn ${fastq_trimmed}/${sample}/R1.fastq.gz --outFilterType BySJout --outFilterMultimapNmax 200 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --limitOutSJcollapsed 5000000 --limitIObufferSize 200000000 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate --outFileNamePrefix alignment/${sample}/ --limitBAMsortRAM 2000000000

# index
samtools index -@ ${nrThreads} alignment/${sample}/"Aligned.sortedByCoord.out.bam" 

# collapsing
echo "deduplication"
mkdir -p deduplicated/${sample}/
umi_tools dedup -I alignment/${sample}/"Aligned.sortedByCoord.out.bam" -S deduplicated/${sample}/"Aligned.sortedByCoord.out.bam" --multimapping-detection-method=NH --output-stats=deduplicated/${sample}/deduplicated.txt --log=deduplicated/${sample}/deduplication.log

tmpdir=$(mktemp -d ${sample}XXXX)

# counting unique alignments
echo "counting"
mkdir -p counting/${sample}/
featureCounts -s 1 -T ${nrThreads} -t ${counting_feature} -g ${identifying_attribute} -a ${gtfFile} -o ${tmpdir}/unique.count ${alignment_deduplicated}/${sample}/"Aligned.sortedByCoord.out.bam" &&  awk 'NR>=3{printf("%s\t%s\n",$1,$NF) > "counting/'${sample}'/unique.count"}' ${tmpdir}/unique.count && rm ${tmpdir}/unique.count*

# counting all alignments, including multimapping reads
featureCounts -s 1 -T ${nrThreads} -t ${counting_feature} -g ${identifying_attribute} -a ${gtfFile} -o ${tmpdir}/multimapper.count -M ${alignment_deduplicated}/${sample}/"Aligned.sortedByCoord.out.bam" &&  awk 'NR>=3{printf("%s\t%s\n",$1,$NF) > "counting/'${sample}'/multimapper.count"}' ${tmpdir}/multimapper.count && rm ${tmpdir}/multimapper.count*

# counting all alignments, including multimapping reads, but divide counts by the number of multimappers
featureCounts -s 1 -T ${nrThreads} -t ${counting_feature} -g ${identifying_attribute} -a ${gtfFile} -o ${tmpdir}/multimapper_avg.count -M --fraction ${alignment_deduplicated}/${sample}/"Aligned.sortedByCoord.out.bam" &&  awk 'NR>=3{printf("%s\t%s\n",$1,$NF) > "counting/'${sample}'/multimapper_avg.count"}' ${tmpdir}/multimapper_avg.count && rm ${tmpdir}/multimapper_avg.count*

rm -d ${tmpdir} # NOTE: this step fails if the directory is not empty, for instance if the previous commands did not remove the temporary files due to some kind of error

popd
