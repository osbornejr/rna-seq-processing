shell.executable
("/bin/bash")
SAMPLES=expand("{TYPE}/{TYPE}_sample_{N}/{TYPE}_sample_{N}",TYPE=["polyA+","polyA-"],N=[str(n) for n in range(13)][1:])
SAMPLEDIRS=expand("{TYPE}/{TYPE}_sample_{N}",TYPE=["polyA+","polyA-"],N=[str(n) for n in range(13)][1:])
IDS=expand("{TYPE}_sample_{N}",TYPE=["polyA+","polyA-"],N=[str(n) for n in range(13)][1:])

basedir=workflow.basedir+'/'
#TODO put long shell commands into scripts and source?
#TODO split Snakefile into two separate sub workflows, de novo and reference assembly 

rule all:
	input:					
		expand('output-data/{X}/{SAMPLE}_RSEM.{X}.results',SAMPLE=IDS,X=["isoforms","genes"])

include: "de-novo-assembly/Snakefile"
include: "reference-guided-assembly/Snakefile" 		
include: "transcript-classification/Snakefile"	
	
rule clean:
# all in one cleaning of fastq files with fastp. Could alternatively do this with trimmomatic? But need specific adapter sequences etc. TODO double check quality independently after running fastp with fastqc
	input:
		r1="/g/data3/ra94/raw-reads/{SAMPLE}_read_1.fastq.gz",
		r2="/g/data3/ra94/raw-reads/{SAMPLE}_read_2.fastq.gz"
	output:
		r1=temp("clean-reads/{SAMPLE}_read_1_fastp.fastq.gz"),
		r2=temp("clean-reads/{SAMPLE}_read_2_fastp.fastq.gz"),
		html="reports/{SAMPLE}_fastp.html",
		json="reports/{SAMPLE}_fastp.json"
		#other fastp outputs?
		# threads?
	#conda: "environment.yml"
	shell:
		'set +u && '
		'eval "$(conda shell.bash hook)" && '
		'conda activate rna-seq && '
		'set -u && '
		'fastp '
		'-i {input.r1} -I {input.r2} '
		'-o {output.r1} -O {output.r2} '
		'-h {output.html} -j {output.json} '
		'-w {threads} '

