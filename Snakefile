shell.executable
("/bin/bash")

configfile: "config.yaml" 
basedir=workflow.basedir+'/'
#TODO put long shell commands into scripts and source?

wildcard_constraints: 
	sample = '|'.join(config["samples"]) 

include: "rules/"+config["assembly_method"]+"-assembly-rules"
	
rule all:
	input:					
		expand('output-data/{read_type}/{sample}_rsem.{read_type}.results',sample=config["samples"],read_type=["isoforms","genes"])
	
rule clean:
# all in one cleaning of fastq files with fastp. Could alternatively do this with trimmomatic? But need specific adapter sequences etc. TODO double check quality independently after running fastp with fastqc
	input:
		r1=config["path_to_files"]+"/{sample}/{sample}_read_1.fastq.gz",
		r2=config["path_to_files"]+"/{sample}/{sample}_read_2.fastq.gz"
	output: 
		r1="clean-reads/{sample}/{sample}_read_1_fastp.fastq.gz",
		r2="clean-reads/{sample}/{sample}_read_2_fastp.fastq.gz",
		html="reports/{sample}/{sample}_fastp.html",
		json="reports/{sample}/{sample}_fastp.json"
		#other fastp outputs?
		# threads?
	shell: 
		'fastp '
		'-i {input.r1} -I {input.r2} '
		'-o {output.r1} -O {output.r2} '
		'-h {output.html} -j {output.json} '
		'-w {threads} '

