shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml" 

FILES = json.load(open(config['SAMPLES_JSON']))
CLUSTER = json.load(open(config['CLUSTER_JSON']))
white_list = config['white_list']
kallisto_INDEX = config['kallisto_INDEX']
SAMPLES = sorted(FILES.keys())

TARGETS = []

## constructe the target if the inputs are fastqs

# ALL_QC = ["07_multiQC/multiQC_log.html"]


TARGETS.extend( expand("02_barcode_info/{sample}_raw_barcode_count.txt", sample = SAMPLES))
# TARGETS.extend( expand("03_corrected/{sample}_corrected_L001_R1_001.fastq", sample = SAMPLES))
TARGETS.extend( expand("04_count/{sample}/output.bus", sample = SAMPLES))
localrules: all
# localrules will let the rule run locally rather than submitting to cluster
# computing nodes, this is for very small jobs

rule all:
	input: TARGETS


rule polyT_selection: # only keep the read2 with polyT structure
	input:
		r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
		r2 = lambda wildcards: FILES[wildcards.sample]['R2']
	output:
		"01_polyA_seq/{sample}_1.fastq.gz" ,
		"01_polyA_seq/{sample}_2.fastq.gz",
		"01_polyA_seq/{sample}_withoutPolyT_1.fastq.gz",
		"01_polyA_seq/{sample}_withoutPolyT_2.fastq.gz", "00_log/{sample}_polyA_fre_cutadapt.log"
	params:
		jobname = "{sample}"
	threads : 12
	# group: "mygroup"
	message: "trim fastqs {input}: {threads} threads"
	shell: # -g 5' regular polyT 10 length required in second reads 
		"""
		cutadapt  -Z -j {threads}    \
        -G T{{10}}  --action none \
           -o {output[0]} -p {output[1]} \
           --untrimmed-output {output[2]} --untrimmed-paired-output {output[3]} \
            {input[0]} {input[1]}  > {output[4]} 
		"""

rule extract_barcode_umi: # trim_reads make it compariable with kallisto pipeline, keep barcode and umi in reads2, and genomic in read1
    input:
        r1 = "01_polyA_seq/{sample}_1.fastq.gz" ,
        r2 = "01_polyA_seq/{sample}_2.fastq.gz"
    output: 
        r1 = ("01_polyA_trimmed/{sample}_index_L001_R1_001.fastq"),
        r2 = ("01_polyA_trimmed/{sample}_index_L001_R2_001.fastq")
    script:
        "script/raw_fq_update_scRNA.py"


rule barcode_QC: ## examine the barcode , first 18 bp of R2 reads
    input: "01_polyA_trimmed/{sample}_index_L001_R2_001.fastq"  ## contain the umi
    output: "02_barcode_info/{sample}_raw_barcode_count.txt"
    shell:
        """
         awk  '{{if(NR%4==2) print substr($0,1,18)}}' {input} | sort | uniq -c | sort -nr   > {output} 
        """

rule find_right_barcodes: ## update the barcode
    input: "02_barcode_info/{sample}_raw_barcode_count.txt"
    output: sum = "02_barcode_info/{sample}.barcode_final_summary",
            map = "02_barcode_info/{sample}.barcode_final_map",
            log = "00_log/{sample}.barcode_log",
    script:
        "script/barcode_hash_RNA.py"


rule read2_barcode_correction:
    input :
        r1 = "01_polyA_trimmed/{sample}_index_L001_R1_001.fastq",
        r2 = "01_polyA_trimmed/{sample}_index_L001_R2_001.fastq",
        map = "02_barcode_info/{sample}.barcode_final_map"
    output :
        r1 = "03_corrected/{sample}_corrected_L001_R1_001.fastq",
        r2 = "03_corrected/{sample}_corrected_L001_R2_001.fastq",
        log = "00_log/{sample}.R2_reads_correction_log"
    script:
        "script/fq_barcode_correction.py"

rule r1_zip:
    input  : "03_corrected/{sample}_corrected_L001_R1_001.fastq"
    output : "03_corrected/{sample}_corrected_L001_R1_001.fastq.gz"
    threads: 11
    shell:
        "pigz -p {threads} {input}"

rule r2_zip:
    input  : "03_corrected/{sample}_corrected_L001_R2_001.fastq"
    output : "03_corrected/{sample}_corrected_L001_R2_001.fastq.gz"
    threads: 11
    shell:
        "pigz -p {threads} {input}"

rule scRNA_count:  ## use the R1 and corrected R2 reads.
    input:
        r2 = "03_corrected/{sample}_corrected_L001_R2_001.fastq.gz",
        r1 = "03_corrected/{sample}_corrected_L001_R1_001.fastq.gz"
    output: directory("04_count/{sample}"), "04_count/{sample}/output.bus"
    threads: 11
    message: "kallisto {input.r1}: {threads} threads"
    log:
         "00_log/{sample}.kallisto"
    shell:
        """
        kallisto bus  -t {threads} -n -i {kallisto_INDEX}  \
        -x 0,0,18:0,18,30:1,0,0 \
        -o {output[0]} \
        {input[r2]}  {input[r1]}   > {log} 2>&1
        """

 # kallisto bus -i /datacommons/ydiaolab/genome_ref/kallisto_index/homo_sapiens_mRNA/transcriptome.idx -n \
 # -x 0,0,18:0,18,30:1,0,0 -o test RNA_1_sci_HiCAR_200925_corrected_L001_R2_001.fastq.gz RNA_1_sci_HiCAR_200925_corrected_L001_R1_001.fastq.gz
