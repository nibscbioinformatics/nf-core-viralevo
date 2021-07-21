# nf-core-viralevo: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the 'results' directory after the pipeline has finished. All paths are relative to the top-level 'results' directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [Preprocessing](#preprocessing)
  * [FastQC](#fastqc) - Raw read QC
  * [Cutadapt](#cutadapt) - Adapter and quality trimming
  * [FastQC](#fastqc) - Trimmed read QC
* [Alignment](#alignment)
  * [BWA-mem](#bwa-mem) - Alignment to the reference genome
* [Alignment post-processing](#alignment-post-processing)
  * [SAMtools](#samtools) - Sort and index alignments, and stats 
  * [LoFreq indelqual](#lofreq_indelqual) - Indel quality insertion in the BAM files
  * [SAMtools](#samtools) - Index bams (containing indel quality) and reference genome fasta
* [Variant calling with LoFreq](#variant-calling-with-lofreq)
  * [LoFreq call-parallel](#lofreq_call-parallel) - Variant calling with LoFreq
* [Variant calling with iVAR](#variant-calling-with-ivar)
  * [iVAR trim](#ivar_trim) - Primer trimming, sort, index and stats on BAM files (containing indel quality)
  * [iVAR](#ivar) - Variant calling with iVAR
  * [TSV2VCF](#tsv2vcf) - Convert TSV files to VCF format
* [Variant annotation](#variant-annotation)
  * [SnpEff](#snpeff) - variant annotation
* [Variant table generation](#variant-table-generation)
* [Filtered VCF generation](#filtered-vcf-generation)   
* [Quality control](#quality-control)
  * [MultiQC](#multiqc) - Aggregate report describing results and QC from the pipeline
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Preprocessing

### FastQC

<details markdown="1">
<summary>Output files</summary>

* `QC/fastqc_raw/`
  * `*_fastqc.html`: FastQC report containing quality metrics for your untrimmed raw fastq files.
  * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

**NB:** The FastQC plots in this directory are generated relative to the raw, input reads. They may contain regions of low quality.
</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### Cutadapt

<details markdown="1">
<summary>Output files</summary>

* `QC/cutadapt/`
  * `*.fastq.gz`: The trimmed/modified fastq reads. These files are NOT published in the pipeline, therefore, you will not find them in the directory.
  * `*.cutadapt.log`: Cutadapt log file containing number and percentage of basepairs processed and trimmed.

</details>

[Cutadapt](https://cutadapt.readthedocs.io/en/stable/) finds and removes adapter sequences, primers, poly-A tails and other types of low quality sequence from your high-throughput sequencing reads.

### FastQC

<details markdown="1">
<summary>Output files</summary>

* `QC/fastqc_trimmed/`
  * `*_fastqc.html`: FastQC report containing quality metrics for your trimmed fastq files.
  * `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

**NB:** The FastQC plots in this directory are generated relative to the trimmed reads. The regions of low quality have been removed.
</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads.

## Alignment

### BWA-mem

<details markdown="1">
<summary>Output files</summary>

* `alignment/bwa/index/bwamem2/`
  * `*.{0123,amb,ann,bwt.2bit.64,pac}`: BWA genome index files
* `alignment/bwa/`
  * `*.bam`: bam file for each sample.

</details>

[BWA-mem](https://github.com/lh3/bwa) is a software package for mapping DNA sequences against a large reference genome, such as the human genome.

## Alignment post-processing

### SAMtools

<details markdown="1">
<summary>Output files</summary>

**NB:** Please note that the SAMtools' sorted and indexed files are NOT published in the pipeline. Therefore, you won't find them.

* `alignment/bwa/samtools_stats/`
  * SAMtools `<SAMPLE>.sorted.bam.flagstat`, `<SAMPLE>.sorted.bam.idxstats` and `<SAMPLE>.sorted.bam.stats` files generated from the alignment files.

</details>

### LoFreq indelqual

<details markdown="1">
<summary>Output files</summary>

* `variants/lowfreq/indelqual/`
  * `*.indelqual.bam.`: bam file with indel qualities

</details>

[LoFreq_indelqual](https://csb5.github.io/lofreq/commands/) adds indel qualities to bam files.

### SAMtools

<details markdown="1">
<summary>Output files</summary>

**NB:** Please note that the SAMtools' indexed *indelqual.bam files and the genome fasta index are NOT published in the pipeline. Therefore, you won't find them.

</details>

## Variant calling with LoFreq

### LoFreq call-parallel

<details markdown="1">
<summary>Output files</summary>

* `variants/lowfreq/call_parallel/`
  * `*_lofreq.vcf`: VCF file containing variant calls

</details>

[LoFreq_call-parallel](https://csb5.github.io/lofreq/commands/) calls variants using multiple processors.

## Variant calling with iVAR

### iVAR trim

<details markdown="1">
<summary>Output files</summary>

* `variants/ivar/trim/`
  * `*log.`: log file containing trimmed read information

</details>

[iVAR trim](https://github.com/andersen-lab/ivar) trims reads in aligned bam.

### iVAR

<details markdown="1">
<summary>Output files</summary>

* `variants/ivar/tsv/`
  * `*tsv.`: TSV files containing variants calls

</details>

[iVAR](https://github.com/andersen-lab/ivar) is a computational package that contains functions broadly useful for viral SNV identification.

### TSV2VCF

<details markdown="1">
<summary>Output files</summary>

* `variants/ivar/tsv2vcf/`
  * `*_ivar.vcf`: VCF files containing variants calls

</details>

## Variant annotation

### SnpEff

<details markdown="1">
<summary>Output files</summary>

* `variants/snpeff/`
  * SnpEff `<SAMPLE>_<caller>.snpeff.csv`, `<SAMPLE>_<caller>.snpeff.genes.txt` and `<SAMPLE>_<caller>.snpeff.summary.html` files
* `variants/snpeff/vcf`
  *`<SAMPLE>_<caller>_annotated.vcf`: Annotated VCF files

</details>

[SnpEff](http://pcingola.github.io/SnpEff/) is a genome variant annotation tool.


## Variant table generation

<details markdown="1">
<summary>Output files</summary>

* `variants/vartable`
  * `varianttable.csv`: Variant table generated from annotated VCF files 

</details>

## Filtered vcf generation

<details markdown="1">
<summary>Output files</summary>

* `variants/vartable/filteredvars/`
  * `<SAMPLE>_<caller>_filtered.vcf`: Filtered vcf file for each input VCF file

</details>

## Quality control
 
### MultiQC

<details markdown="1">
<summary>Output files</summary>

* `QC/multiqc/`
  * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  * `multiqc_plots/`: directory containing static images from the report in various formats.

</details>


[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

* `pipeline_info/`
  * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  * Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.


