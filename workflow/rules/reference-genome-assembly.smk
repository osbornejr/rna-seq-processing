#Transcriptome assembly using reference genome

#This sets up the reference geneome that STAR will align to. Only needs to be done once for all subsequent STAR alignment jobs 
#NOTE: irrelevant at present, but might be necessary for quantification that doesn't use RSEM (i.e. per gene rather than per-transcript)
#rule star_reference:
#	input:
#		fasta="reference-input/'+config["ref_genome"]+'/'+config["ref_genome"]+'.fasta",
#		gtf="reference-input/'+config["ref_genome"]+'/'+config["ref_genome"]+'.gtf"
#	output:
#		'reference-index/'+config["ref_genome"]+'/SAindex'
#	#conda: "environment.yml"
#	shell:
#		'set +u && '
#		'eval "$(conda shell.bash hook)" && '
#		'conda activate rna-seq && '
#		'set -u && '
#		'mkdir -p {output} && '
#		'STAR --runThreadN {threads} '
#		'--runMode genomeGenerate '
#		'--genomeDir {output} '
#		'--genomeFastaFiles {input.fasta} '
#		'--sjdbGTFfile {input.gtf} '
#		'--sjdbOverhang 100 '


rule rsem_reference:
        input:
                fasta = 'reference-input/'+config["ref_genome"]+'/'+config["ref_genome"]+'.fasta',
                gtf = 'reference-input/'+config["ref_genome"]+'/'+config["ref_genome"]+'.gtf'

        params:
                outdir = 'reference-index/'+config["ref_genome"]+'/'+config["ref_genome"]+'',
		threads = 16

        output:
                grp = 'reference-index/'+config["ref_genome"]+'/'+config["ref_genome"]+'.grp',
                ti = 'reference-index/'+config["ref_genome"]+'/'+config["ref_genome"]+'.ti',
                seq = 'reference-index/'+config["ref_genome"]+'/'+config["ref_genome"]+'.seq',
                chrlist = 'reference-index/'+config["ref_genome"]+'/'+config["ref_genome"]+'.chrlist',
                SAindex = 'reference-index/'+config["ref_genome"]+'/SAindex',
		#refind = directory('reference-index/'+config["ref_genome"]+'')

        #conda: "environment.yml"

        shell:
                'set +u && '
                'eval "$(conda shell.bash hook)" && '
                'conda activate rna-seq && '
                'set -u && '
                'rsem-prepare-reference {input.fasta} --star -p {params.threads} --gtf {input.gtf} {params.outdir}'

##2025 update-- have left the single pass align as is with temp dir in case it is needed again (eg on gadi) other rules changed to work without temp dir on aws
rule align:
	input:
		r1 = 'clean-reads/{sample}/{sample}_read_1_fastp.fastq.gz',
		r2 = 'clean-reads/{sample}/{sample}_read_2_fastp.fastq.gz',
		ref_genome = 'reference-index/'+config["ref_genome"]

	params:
		indir = 'reference-index/'+config["ref_genome"],
		tempdir = '$PBS_JOBFS',
		outdir = 'aligned-reads/{sample}/{sample}',
		threads = 16
	output:
		'aligned-reads/{sample}/{sample}_single_pass/Aligned.sortedByCoord.out.bam',
		'aligned-reads/{sample}/{sample}_single_pass/Aligned.toTranscriptome.out.bam'
	shell:
		#'rm -rf {params.outdir} && '
		#'mkdir {params.outdir} && '
		#'cd {params.outdir} && '
		'STAR --runThreadN {params.threads} '
		'--genomeDir {params.indir} '
		'--outFileNamePrefix {params.tempdir}/ '
		'--readFilesIn {input.r1} {input.r2} '
		'--readFilesCommand zcat '
		'--outSAMtype BAM SortedByCoordinate '
		'--quantMode TranscriptomeSAM GeneCounts && '
		'mv {params.tempdir}/* {params.outdir}'

rule align_pass_1:
	input:
		r1 = 'clean-reads/{sample}/{sample}_read_1_fastp.fastq.gz',
		r2 = 'clean-reads/{sample}/{sample}_read_2_fastp.fastq.gz',
		ref_genome = 'reference-index/'+config["ref_genome"]+'/SAindex'
	params:
		indir  = 'reference-index/'+config["ref_genome"],
		tempdir = '$PBS_JOBFS',
		outdir = 'aligned-reads/{sample}/{sample}_pass_1',
		rmbam = 'aligned-reads/{sample}/{sample}_pass_1/Aligned.out.bam',
		threads = 16
	output:
		'aligned-reads/{sample}/{sample}_pass_1/SJ.out.tab'
	#conda: "environment.yml"
	shell:		
		#'set +u && '
		#'eval "$(conda shell.bash hook)" && '
		#'conda activate rna-seq && '
		#'set -u && '
		#rm -rf {params.outdir} && \
		#mkdir {params.outdir} && \
		#cd {params.outdir} && \
		'mkdir -p {params.outdir} &&'
		'STAR --runThreadN {params.threads} '
		'--outFileNamePrefix {params.outdir}/ '
		'--genomeDir {params.indir} '
		'--readFilesIn {input.r1} {input.r2} '
		'--readFilesCommand zcat '
		'--outSAMtype BAM Unsorted && rm {params.rmbam} '
		#'mv {params.tempdir}/* {params.outdir}' 

rule filter:
	input:
		'aligned-reads/{sample}/{sample}_pass_1/SJ.out.tab',
	output:
		'splice-junctions/{sample}/{sample}_pass_1_SJ.filtered.tab'
	shell:
		'awk "{{if (\$7 >= 3) print \$0 }}" {input[0]} > {input[0]}.filtered && '
		'mv {input[0]}.filtered {output}'

rule align_pass_2:
	input:
		r1 = 'clean-reads/{sample}/{sample}_read_1_fastp.fastq.gz',
		r2 = 'clean-reads/{sample}/{sample}_read_2_fastp.fastq.gz',
		sj_files =expand('splice-junctions/{sample}/{sample}_pass_1_SJ.filtered.tab',sample=config["samples"]),
		ref_genome ='reference-index/'+config["ref_genome"]+'/SAindex'

	params:
		indir  = 'reference-index/'+config["ref_genome"],
		tempdir = '$PBS_JOBFS',
		outdir = 'aligned-reads/{sample}/{sample}_pass_2',
		threads = 16
	output:
		#temp('aligned-reads/{sample}/{sample}_pass_2/Aligned.sortedByCoord.out.bam'),
		'aligned-reads/{sample}/{sample}_pass_2/Aligned.toTranscriptome.out.bam'
	shell:
		#'rm -rf {params.outdir} && '
		#'mkdir {params.outdir} && '
		#'cd {params.outdir} && '
		'mkdir -p {params.outdir} &&'
		'STAR --runThreadN {params.threads} '
		'--genomeDir {params.indir} '
		'--outFileNamePrefix {params.outdir}/ '
		'--readFilesIn {input.r1} {input.r2} '
		'--readFilesCommand zcat '
		'--outSAMtype BAM Unsorted '
		#'--outSAMtype BAM SortedByCoordinate '
		'--sjdbFileChrStartEnd {input.sj_files} '
		'--quantMode TranscriptomeSAM GeneCounts '
		#'mv {params.tempdir}/* {params.outdir}' 


rule rsem:
	input: 
	 	bam = 'aligned-reads/{sample}/{sample}_pass_2/Aligned.toTranscriptome.out.bam',
		ref = expand('reference-index/'+config["ref_genome"]+'/'+config["ref_genome"]+'.{ext}',ext=['grp','ti','seq','chrlist'])

	params: 
		indir  = 'reference-index/'+config["ref_genome"]+'/'+config["ref_genome"], 
		outdir = 'transcript-counts/{sample}',
                outpre = '{sample}_RSEM',
		threads = 16

	output: 'output-data/isoforms/{sample}_RSEM.isoforms.results',
		'output-data/genes/{sample}_RSEM.genes.results'
	
	log:	'logs/rsem-quantification/{sample}_RSEM.log'

	shell:
		'mkdir -p {params.outdir} &&'
		'rsem-calculate-expression --bam --paired-end --calc-pme -p {params.threads} {input.bam} {params.indir} {params.outdir}/{params.outpre} && '
		'mkdir -p output-data/genes && '
		'mkdir -p output-data/isoforms && '
		'cp {params.outdir}/{params.outpre}.genes.results output-data/genes/ && '
		'cp {params.outdir}/{params.outpre}.isoforms.results output-data/isoforms/ '
		#'find ./transcript-counts/ -name "*.genes.results" -exec cp {{}} output-data/genes/ \; && '
		#'find ./transcript-counts/ -name "*.isoforms.results" -exec cp {{}} output-data/isoforms/ \; '

#Use htseq rule to quantify GENE counts. 2025 Note: this is the only rule that requires sorted BAM, so only turn on for this here.
#rule htseq:
#	input:
#		bam = 'aligned-reads/{sample}/{sample}_pass_2/Aligned.sortedByCoord.out.bam',
#		gff = expand('reference-input/'+config["ref_genome"]+'/'+config["ref_genome"]+'.gff3',REF_GENOME=REF_GENOME)
#	output:
#		'gene-counts/{sample}/{sample}_HTSeq_union_gff3_no_gene_ID.log',
#		'gene-counts/{sample}/{sample}_HTSeq.csv'
#	#conda: "environment.yml"
#	shell:
#		'htseq-count -m union -s no -t gene -i ID -r pos -f bam {input.bam} {input.gff} &> {output[0]} && '
#		'grep {GENE_CODE} {output[0]} | sed "s/gene://g" > {output[1]}'
