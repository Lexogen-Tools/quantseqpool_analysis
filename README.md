# QuantSeqPoolWF

This script is a basic analysis on QuantSeqPool data using demultiplexed fastq files OR demultiplexes this data using idemux from fastq files and then analyzes the samples.

The dependencies of the script are listed in environment.yml. You can use the file to install them using conda (anaconda, miniconda) into a QuantSeqPool environment e.g. with: 
```
conda env create -f environment.yml
conda activate QuantSeqPool
```
Alternatively, you can use all tools with your own installations, provided you have them in your PATH and the versions are compatible with the required versions. 

The user is required to have write access to the provided base directory or the current working directory, which acts as a temp directory for intermediate results and output directory.
The script assumes that the STAR genome was already built prior to the execution of the script. You can do this with a call like 
```
mkdir -p /path/to/star/genome/directory && STAR --runMode genomeGenerate --runThreadN 10 --genomeDir /path/to/star/genome/directory --genomeFastaFiles /path/to/reference/sequences/multifasta.fa --sjdbGTFfile /path/to/reference/annotation/in/gtf/format.gtf --sjdbOverhang 99
```
The STAR genome has to match the version of the STAR aligner binary that is used to align the samples. Please be aware that the STAR aligner requires a certain amount of RAM available for your processes depending on the size of your reference genome.
Also, please be aware that the script assumes the presence of exon features in the .gtf file. If this is not true e.g. in bacteria species, change the respective variable value to the correct feature type. If the reference annotation does not identify genes with gene_id, please change the respective variable value to the correct attribute name.

You can either provide a samplesheet to define the samples of the analysis or you place the raw, demultiplexed fastq files of each sample into a directory hierarchy as
```
 fastq / [sample1] / R1.fastq.gz
 fastq / [sample1] / R2.fastq.gz
 fastq / [sample2] / R1.fastq.gz
 fastq / [sample2] / R2.fastq.gz
 ...
```
into the [basedir]. In the former scenario the script will read the sample information from the samplesheet and create the fastq directory accordingly after demultiplexing the data with idemux.

Example call of the script:
```
./conductQuantSeqPoolAnalysis.sh -g /home/data/DAP_resources/hsa_GRCh38.102_ERCC_SIRV/annotation_organism_ercc_sirv_biotyped.gtf -d /home/data/DAP_resources/hsa_GRCh38.102_ERCC_SIRV/genomic_data/ -b workdir -s workdir/SampleSheet.csv -t 18
```
, where [workdir] either contains the fastq directory, in which the fastq files are stored as mentioned above.
Alternatively, you provide R1 and R2 via the respective arguments --completeRawR1 and --completeRawR2 and provide as sample sheet together with the working directory via -b. 

In both cases you need write permissions for the work directory
