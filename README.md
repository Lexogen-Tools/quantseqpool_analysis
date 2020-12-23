# QuantSeqPoolWF
**This analysis script is implemented as a bash routine and is intended to work in combination with Lexogen's QuantSeq-Pool Sample-Barcoded 3′ mRNA-Seq Library Prep Kit for Illumina**

QuantSeq-Pool is the optimal solution for gene expression profiling for large screening projects using sample barcoding, early pooling, and batch processing of up to 96 samples in one reaction providing a workflow that is easily scalable for multiplexing up to 36,864 samples. Each sample is identified via one of 96 Lexogen i1 sample indices (for further information please visit https://www.lexogen.com/quantseq-pool-sample-barcoded-3mrna-sequencing/). The analysis script extracts the samples according to their sample indicies from an pool of reads in a UDI barcode combination and is then further processed with bioinformatics tools to trim the raw data, align the reads and count the expressed genes.

## Requirements
The script [QuantSeqPoolAnalysis.sh](QuantSeqPoolAnalysis.sh) uses the following software tools and scripts in its routine:
```
    - idemux (https://github.com/Lexogen-Tools/idemux)
    - umi_tools (https://github.com/CGATOxford/UMI-tools)
    - cutadapt (https://cutadapt.readthedocs.io/en/stable/)
    - STAR aligner (https://github.com/alexdobin/STAR)
    - samtools (https://www.htslib.org/)
    - featureCounts (http://subread.sourceforge.net/)
    - GNU awk (https://www.gnu.org/software/gawk/)
```

The tools can either be installed and made available in the user's PATH or installed into a conda environment and used in this environment. If you use your own installations, please make sure that the versions of your installed tools correspond to the versions in [environment.yml](environment.yml). Also, please ensure that the tools are available in your path.

You can use conda (anaconda, miniconda) to install the tools into a new conda environment and use this environment, e.g. named 'QuantSeqPool', by calling 
```
$ conda env create -f /path/to/environment.yml && conda activate QuantSeqPool
```
## How to run this script
Provided you have the required tools installed, cloned the content of the git repository onto your machine and changed into this cloned directory, you can start the script with a call like
```
./conductQuantSeqPoolAnalysis.sh -g path/to/reference/annotation.gtf -d path/to/STAR/genome --rawR1 path/to/R1.fastq.gz --rawR2 path/to/R2.fastq.gz -s path/to/sample-sheet.csv -o path/to/output/directory -t [nr of parallel jobs]
```
**path/to/reference/annotation.gtf** is the path to an reference annotation in Ensembl-style gtf format (gff version 2). The file contains the annotations of features on a reference genome. The annotation has to match the reference that was used to generate the STAR genome.

**path/to/STAR/genome** is a pre-built STAR genome reference. This reference is supposed to match the organism that is to be analyzed with the the experiment.

**path/to/demultiplexed_data/R1.fastq.gz** is the path to a gzipped fastq file, which contains the first read of the read-pairs of the sample(s) of this experiment.

**path/to/demultiplexed_data/R2.fastq.gz** is the path to a gzipped fastq file, which contains the second read of the read-pairs of the sample(s) of this experiment.

**path/to/sample-sheet.csv** is the path to a sample sheet that lists the samples in the fastq file and their associated Lexogen UDI and sample barcodes (see section "Demultiplexing with idemux").

**path/to/output_directory** is the path to a directory to write intermediate and final analysis results.

**[nr of parallel jobs]** is the number of threads which are allowed for each tool call benefits from threading.
## The analysis routine
The script is a bash routine that processes a pair of i7/i5 demultiplexed, gzipped fastq files. The script allows to analyze up to 96 samples within a single file pair.

The following sections outline the analysis steps implemented in the script.
### Demultiplexing with idemux
```
$ idemux --r1 ${pathR1} --r2 ${pathR2} --sample-sheet ${samplesheet} --out idemultiplexed/
```
In addition to the use of i7/i5 indices, the library preparation introduces Lexogen's i1 sample barcodes (see https://www.lexogen.com/quantseq-pool-sample-barcoded-3mrna-sequencing/) to increase the number of simultaneously analyzed samples. The script starts after i7/i5 demultiplexing has been performed, e.g. with bcl2fastq. The first step in the script uses iDemux to demultiplex the gzipped fastq file according to a sample sheet which lists the names of the samples together with their associated Lexogen i1 sample barcode. 
**The sample sheet in this repository, [sample-sheet.csv](sample-sheet.csv), is filled with barcode combinations that match the supplied test case, please edit the file to enter the barcode combinations that are present in your experiment before you start the analysis**
Please note that sample names are used for the creation of temporary files and directories and must therefore conform with file name requirements in a *nix environment.

If you are familiar with the use of iDemux and know how to demultiplex your sequencing run into a single pair of fastq files, preserving the index information in your read headers, you can instead 

* provide this single demultiplexed fastq file pair to the analysis script and 
* list all the samples of your run in the sample-sheet.csv file together with their complete combination of i7, i5, and i1 barcode sequences

The analysis script will then perform demultiplexing of i7/i5 and i1 indices together and analyze all corresponding samples in one go. For further information on how to use iDemux and demultiplex into a single fastq file preserving the index information in read headers, please refer to the iDemuxCPP Readme (https://github.com/Lexogen-Tools/idemux).

### Extracting UMIs with umi_tools
```
$ umi_tools extract --extract-method=string --bc-pattern X --bc-pattern2 NNNNNNNNNN -L path/to/sample/extraction_log -S path/to/output/extracted/sample/R1.fastq.gz \
> -I path/to/output/fastq/raw/R1.fastq.gz --read2-in=path/to/output/fastq/sample/R2.fastq.gz --read2-out=/dev/null
```
umi_tools extracts the 10N long UMI sequence from the second read in the read pair and writes this information into the fastq read ID (sequence ID) of the first read. The second read has no more useable information and is not used past this analysis step.
### Trimming with cutadapt
```
$ cutadapt --quiet -m 20 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 path/to/output/extracted/sample/R1.fastq.gz | \
> cutadapt --quiet -m 20 -O 3 --nextseq-trim=10 -a "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - | \
> cutadapt --quiet -m 20 -O 3 -a "r1polyA=A{18}" - | \
> cutadapt --quiet -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed -o path/to/output/trimmed/sample/R1.fastq.gz -
```
Prior to alignment reads are trimmed with cutadapt. Reads which result in a length <20 or consist entirely of adapter sequences are removed.
### Alignment with the STAR aligner
```
$ STAR --runThreadN ${nrThreads} --readFilesCommand zcat --genomeDir /path/to/star/genome --readFilesIn path/to/output/trimmed/sample/R1.fastq.gz \
> --outFilterType BySJout --outFilterMultimapNmax 200 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
> --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --limitOutSJcollapsed 5000000 \
> --limitIObufferSize 200000000 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate --outFileNamePrefix path/to/output/alignment/sample/ \
> --limitBAMsortRAM 2000000000
```
The alignment with STAR requires a pre-built STAR aligner genome. If you do not have such genome for your reference, but you have the **genomic** reference sequences and the annotation of your target species you can use the STAR aligner to generate a genome with a call like
```
$ mkdir -p path/to/star/genome && STAR --runMode genomeGenerate --runThreadN 10 --genomeDir /path/to/star/genome --genomeFastaFiles /path/to/reference/sequences/multifasta.fa \
> --sjdbGTFfile /path/to/reference/annotation.gtf --sjdbOverhang 99
```
To ensure that the STAR genome version matches the version of the STAR binary of this script, we recommend to create the genome index files with the STAR binary that was listed in the conda environment. Please be aware that depending on the size of the reference you will require a significant amount of memory to properly run the alignment (about 30G of RAM on the human reference genome).
Further, if the target organism does not have ordinary splicing events you should omit the arguments --sjdbGTFfile and --sjdbOverang.
### Deduplication with umi_tools
```
$ umi_tools dedup -I path/to/output/alignment/sample/Aligned.sortedByCoord.out.bam -S path/to/output/deduplicated/sample/Aligned.sortedByCoord.out.bam \
> --multimapping-detection-method=NH --output-stats=path/to/output/deduplicated/sample/deduplicated.txt --log=path/to/output/deduplicated/sample/deduplication.log
```
umi_tools is used to remove duplicates of reads based on the location of the alignment and the UMI information in the read header that was extracted from Read 2 during the UMI extraction step.
### Counting with featureCounts
```
$ featureCounts -s 1 -T ${nrThreads} -t exon -g gene_id -a path/to/reference/annotation.gtf -o ${tmpdir}/unique.count path/to/output/deduplicated/sample/Aligned.sortedByCoord.out.bam && \
> awk 'NR>=3{printf("%s\t%s\n",$1,$NF) > "path/to/output/counting/'${sample}'/unique.count"}' ${tmpdir}/unique.count && rm ${tmpdir}/unique.count*
```
We use featureCounts to classify and count the alignments against the reference by using a matching annotation of the reference. The default setting of this quantification is to use the gene_id attribute of the exon features in the annotation. If your target species uses a different term for the expression of features please provide the arguments **--feature** and **--attribute** to tell the script which feature and which attribute it is supposed to use for this analysis step.
### Output
All results and intermediate results can be used for further downstream analysis. The end of this analysis script is the quantification with featureCounts, which is done in three versions 
- counts based on unique alignments to the reference (unique)
- counts based on all alignments against the reference (all)
- counts based on all alignments against the reference, but each alignment that maps multiple times is divided by its number of alignments (all_avg-multimapper)

The results of each counting variation are summarized into a file each, named summary_[type].tsv. The summary is a tabular-formatted text file with tab as field separator. The column headers contain the sample names, the first column contains the ids of the genes.
## Testing the pipeline

To verify that the pipeline has been properly installed, the gzipped fastq file pair in [Resources/fastq/](Resources/fastq/) can be used to test the analysis pipeline. Starting the pipeline with
```
$ ./conductQuantSeqPoolAnalysis.sh -g path/to/reference/annotation.gtf -d path/to/STAR/genome --rawR1 Resources/fastq/R1.fastq.gz --rawR2 Resources/fastq/R2.fastq.gz -s sample-sheet.csv -t 18 -o test
```
should result in summary files similar to the ones you find in [Resources/](Resources/). The example summary files were generated based on these fastq files and by using the reference genome hsa_GRCh38.102 from Ensembl:
- gtf (ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/)
- fasta (ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/)
