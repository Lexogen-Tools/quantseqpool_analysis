#!/usr/bin/env bash

counting_feature=exon
identifying_attribute=gene_id
umidef=N10
nrThreads=1
outputDir=$(pwd)

#usage() { echo "Usage: $0  [-g|--gtfFile <gtfFile>] [-d|--starGenomeDir <dir>] [-c|--rawR1 </path/to/R1.fastq.gz>] [-e|--rawR2 </path/to/R2.fastq.gz>] [-s|--samplesheet </path/to/SampleSheet.csv>] (Optional: [-o|--outputDir <dir>] [-f|--feature <counting_feature>] [-a|--attribute <grouping attribute>] [-t|--threads] <# of threads>)" 1>&2; exit 1; }
usage() { echo "Usage: $0  [-g|--gtfFile <gtfFile>] [-d|--starGenomeDir <dir>] [-s|--samplesheet </path/to/SampleSheet.csv>] (Optional: [-o|--outputDir <dir>] [-f|--feature <counting_feature>] [-a|--attribute <grouping attribute>] [-t|--threads] <# of threads>)" 1>&2; exit 1; }

#OPT=$(getopt -o g:d:o:c:e:s:f:a:t: --long gtfFile:,starGenomeDir:,output:,rawR1:,rawR2:,sampleSheet:,feature:,attribute:,threads:, -- "$@")
OPT=$(getopt -o g:d:o:s:f:a:t: --long gtfFile:,starGenomeDir:,output:,sampleSheet:,feature:,attribute:,threads:, -- "$@")

eval set -- "$OPT"

while true; do
  case "$1" in
    -g | --gtfFile ) gtfFile="$2"; shift 2;;
    -d | --starGenomeDir ) genomeDir="$2"; shift 2;;
    -o | --output ) outputDir="$2"; shift 2;;
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
gtfFile=$(readlink -f ${gtfFile})
genomeDir=$(readlink -f ${genomeDir})
samplesheet=$(readlink -f ${samplesheet})
outputDir=$(readlink -f ${outputDir}) ; mkdir -p ${outputDir}
analysisScript=$(readlink -f $(dirname $(readlink -f $0))/QuantSeqPoolAnalysis.sh)

echo "Reading samplesheet ${samplesheet}"
while read -r sample fq1 fq2; do
  echo "  read: ${sample}, ${fq1}, ${fq2}"
  samples+=("$sample")
  fqs_R1+=($(realpath ${fq1}))
  fqs_R2+=($(realpath ${fq2}))
done < <(tail -n +2 ${samplesheet} | sed 's/,/\t/g')
if [ ${#samples} -lt 1 ]; then
  echo "Error: no samples in SDF file: ${samplesheet}!"
  exit 1
fi

pushd $outputDir

mkdir -p fastq # try to create dir if not exists

echo "${samples[@]}"
echo "R1 ${fqs_R1[@]}"
echo "R2 ${fqs_R2[@]}"

echo "Linking read files..."
for (( j=0; j<${#samples[@]}; j++ ));
do
    mkdir -p fastq/${samples[$j]} ;
    ln -sf ${fqs_R1[$j]} fastq/${samples[$j]}/R1.fastq.gz
    ln -sf ${fqs_R2[$j]} fastq/${samples[$j]}/R2.fastq.gz
done;

# analysis per sample
echo "Analysis per sample"
for sample in ${samples[@]}; do
    echo "  analyzing ${sample}"
    echo "${analysisScript} -s ${sample} -g ${gtfFile} -d ${genomeDir} -f ${counting_feature} -a ${identifying_attribute} -t $nrThreads"
    ${analysisScript} -s ${sample} -g ${gtfFile} -d ${genomeDir} -f ${counting_feature} -a ${identifying_attribute} -t $nrThreads
done

# summary of counts
echo "Summary of counts"
for ctype in unique all all_avg-multimapper ; do
    if [ -f counting/$s/${ctype}.tsv ]; then
        echo "Summarizing counts of ${ctype}"
        (echo -e "id "${samples[@]} | sed 's/ /\t/g'; paste $(for s in ${samples[@]}; do ls counting/$s/${ctype}.tsv; done) | awk '{printf "%s",$1 ; for (ii=2;ii<=NF;ii=ii+2){printf "\t%s",$ii} ; printf "\n" }') > summary_${ctype}.tsv
        echo "Summarized in ${summary_${ctype}.tsv}"
    fi
done

popd
