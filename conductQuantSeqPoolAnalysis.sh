#!/usr/bin/env bash

counting_feature=exon
identifying_attribute=gene_id
umidef=N10
nrThreads=1
baseDir=$(pwd)

usage() { echo "Usage: $0  [-g|--gtfFile <gtfFile>] [-d|--starGenomeDir <dir>] [-b|--baseDir <dir>] [-c|--completeRawR1 </path/to/R1.fastq.gz>] [-e|--completeRawR2 </path/to/R2.fastq.gz>] [-s|--samplesheet </path/to/SampleSheet.csv>] [-f|--feature <counting_feature>] [-a|--attribute <grouping attribute>] [-t|--threads] <# of threads>)" 1>&2; exit 1; }

OPT=$(getopt -o g:d:u:b:c:e:s:f:a:t: --long gtfFile:,starGenomeDir:,baseDir:,completeRawR1:,completeRawR2:,sampleSheet:,feature:,attribute:,threads:, -- "$@")

eval set -- "$OPT"

while true; do
  case "$1" in
    -g | --gtfFile ) gtfFile="$2"; shift 2;;
    -d | --starGenomeDir ) genomeDir="$2"; shift 2;;
    -b | --baseDir ) baseDir="$2"; shift 2;;
    -c | --completeRawR1 ) pathR1="$2"; shift 2;;
    -e | --completeRawR2 ) pathR2="$2"; shift 2;;
    -s | --sampleSheet ) samplesheet="$2"; shift 2;;
    -f | --feature ) counting_feature="$2"; shift 2;;
    -a | --attribute ) identifying_attribute="$2"; shift 2;;
    -t | --threads ) nrThreads="$2"; shift 2;;
    -- ) shift; break;;
    * ) echo "invalid state" ; usage;;
  esac
done

[ $# -gt 0 ] && { echo "no command line arguments outside options." 1>&2; usage; }

[ -z ${gtfFile+x} ] && { echo "gtfFile has to be specified." 1>&2; usage; }
[ -z ${genomeDir+x} ] && { echo "starGenomeDir has to be specified." 1>&2; usage; }
genomeDir=$(readlink -f ${genomeDir})
pathR1=$(readlink -f ${pathR1})
pathR2=$(readlink -f ${pathR2})
gtfFile=$(readlink -f ${gtfFile})
baseDir=$(readlink -f ${baseDir}) ; mkdir -p ${baseDir}
analysisScript=$(readlink -f ./QuantSeqPoolAnalysis.sh)

if [ -z ${samplesheet+x} ] ; then
    echo "no samplesheet provided; it is assumed that the samples were already properly demuliplexed per i7,i5,i1 index combination and written into fastq / \$\{sample\} / r1 (,fastq / \$\{sample\} / r2)."
    pushd $baseDir
    samples=($(basename -a $(find fastq/ -maxdepth 1 -mindepth 1 -type d | sort)| tr '\n' ' '))
else
    samplesheet=$(readlink -f ${samplesheet})
    samples=($(cut -d ',' -f $(awk -F ',' 'NR==1{ for(ii=1;ii<=NF;ii++){if($ii=="sample_name"){print ii}}}' ${samplesheet}) ${samplesheet} | tail -n +2 | tr '\n' ' '))

    pushd $baseDir
    mkdir -p fastq # try to create dir if not exists
    n_samples_in_fastqdir=$(comm -12 <( echo ${samples[@]} | tr ' ' '\n') <( basename -a $(find fastq/ -maxdepth 1 -mindepth 1 -type d | sort) 2>/dev/null ) | wc -l)
    if [ $n_samples_in_fastqdir != ${#samples[@]} ] 2>/dev/null ; then
        # not all samples present in fastq dir -> check if run was provided as argument, then process the run with standard settings for demultiplexing
        if [ -z ${pathR1+x} ] ; then
            echo "path to R1 and R2 has to be specified in order to be processed with idemux."
            echo "please use appropriate demultipexing-software available to you to obtain the reads."
            1>&2; usage;
        fi
        mkdir -p idemultiplexed/
        idemux --r1 ${pathR1} --r2 ${pathR2} --sample-sheet ${samplesheet} --out idemultiplexed/
        for sample in ${samples[@]}; do
            mkdir -p fastq/${sample} ;
            ln -sf ../../idemultiplexed/${sample}_R1.fastq.gz fastq/${sample}/R1.fastq.gz
            ln -sf ../../idemultiplexed/${sample}_R2.fastq.gz fastq/${sample}/R2.fastq.gz
        done ;
    fi
    # else, all samples present as directories in fastq dir, do not run idemux and attempt to process fastq files as if they were already demultiplexed with idemux.
fi

mkdir -p extracted trimmed alignment deduplicated counting
for sample in ${samples[@]}; do
    ${analysisScript} -s ${sample} -g ${gtfFile} -d ${genomeDir} -t $nrThreads
done
