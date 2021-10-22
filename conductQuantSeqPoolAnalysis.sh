#!/usr/bin/env bash

counting_feature=exon
identifying_attribute=gene_id
umidef=N10
nrThreads=1
outputDir=$(pwd)

usage() { echo "Usage: $0  [-g|--gtfFile <gtfFile>] [-d|--starGenomeDir <dir>] [-c|--rawR1 </path/to/R1.fastq.gz>] [-e|--rawR2 </path/to/R2.fastq.gz>] [-s|--samplesheet </path/to/SampleSheet.csv>] (Optional: [-o|--outputDir <dir>] [-f|--feature <counting_feature>] [-a|--attribute <grouping attribute>] [-t|--threads] <# of threads>)" 1>&2; exit 1; }

OPT=$(getopt -o g:d:o:c:e:s:f:a:t: --long gtfFile:,starGenomeDir:,output:,rawR1:,rawR2:,sampleSheet:,feature:,attribute:,threads:, -- "$@")

eval set -- "$OPT"

while true; do
  case "$1" in
    -g | --gtfFile ) gtfFile="$2"; shift 2;;
    -d | --starGenomeDir ) genomeDir="$2"; shift 2;;
    -o | --output ) outputDir="$2"; shift 2;;
    -c | --rawR1 ) pathR1="$2"; shift 2;;
    -e | --rawR2 ) pathR2="$2"; shift 2;;
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
[ -z ${samplesheet+x} ] && { echo "samplesheet has to be specified." 1>&2; usage; }
[ -z ${pathR1+x} ] && { echo "path to R1 has to be specified." 1>&2; usage; }
[ -z ${pathR2+x} ] && { echo "path to R2 has to be specified." 1>&2; usage; }
gtfFile=$(readlink -f ${gtfFile})
genomeDir=$(readlink -f ${genomeDir})
pathR1=$(readlink -f ${pathR1})
pathR2=$(readlink -f ${pathR2})
samplesheet=$(readlink -f ${samplesheet})
outputDir=$(readlink -f ${outputDir}) ; mkdir -p ${outputDir}
analysisScript=$(readlink -f $(dirname $(readlink -f $0))/QuantSeqPoolAnalysis.sh)

samples=($(cut -d ',' -f $(awk -F ',' 'NR==1{ for(ii=1;ii<=NF;ii++){if($ii=="sample_name"){print ii}}}' ${samplesheet}) ${samplesheet} | tail -n +2 | tr '\n' ' '))

pushd $outputDir

mkdir -p fastq # try to create dir if not exists
mkdir -p idemultiplexed/ # try to create dir if not exists
idemux --r1 ${pathR1} --r2 ${pathR2} --sample-sheet ${samplesheet} --out idemultiplexed/
for sample in ${samples[@]}; do
    mkdir -p fastq/${sample} ;
    ln -sf ../../idemultiplexed/${sample}_R1.fastq.gz fastq/${sample}/R1.fastq.gz
    ln -sf ../../idemultiplexed/${sample}_R2.fastq.gz fastq/${sample}/R2.fastq.gz
done ;

# analysis per sample
for sample in ${samples[@]}; do
    ${analysisScript} -s ${sample} -g ${gtfFile} -d ${genomeDir} -f ${counting_feature} -a ${identifying_attribute} -t $nrThreads
done

# summary of counts
for ctype in unique all all_avg-multimapper ; do
	(echo -e "id "${samples[@]} | sed 's/ /\t/g'; paste $(for s in ${samples[@]}; do ls counting/$s/${ctype}.tsv; done) | awk '{printf "%s",$1 ; for (ii=2;ii<=NF;ii=ii+2){printf "\t%s",$ii} ; printf "\n" }') > summary_${ctype}.tsv
done

popd
