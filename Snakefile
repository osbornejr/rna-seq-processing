shell.executable("/bin/bash")
SAMPLES=expand("{TYPE}/{TYPE}_sample_{N}/{TYPE}_sample_{N}",TYPE=["polyA+","polyA-"],N=[str(n) for n in range(13)][1:])
IDS=expand("{TYPE}_sample_{N}",TYPE=["polyA+","polyA-"],N=[str(n) for n in range(13)][1:])
REF_GENOME="Kabuli_UWA-v2.6.3"
GENE_CODE="Ca"


rule all:
	input:
		expand("clean-reads/{SAMPLE}_read_{N}_fastp.fastq.gz",SAMPLE=SAMPLES,N=["1","2"]),
		expand("reports/{SAMPLE}_fastp.{EXT}",SAMPLE=SAMPLES,EXT=["html","json"]),
		#expand('reference-index/{REF_GENOME}_star',REF_GENOME=REF_GENOME),
		expand('aligned-reads/{SAMPLE}_pass_1/SJ.out.tab',SAMPLE=SAMPLES),
		expand('splice-junctions/{SAMPLE}_pass_1_SJ.filtered.tab',SAMPLE=SAMPLES),
		expand('aligned-reads/{SAMPLE}_pass_2/Aligned.sortedByCoord.out.bam',SAMPLE=SAMPLES),
		expand('aligned-reads/{SAMPLE}_pass_2/Aligned.toTranscriptome.out.bam',SAMPLE=SAMPLES),
		expand('reference-index/{REF_GENOME}/{REF_GENOME}.grp',REF_GENOME=REF_GENOME),
		expand('reference-index/{REF_GENOME}/{REF_GENOME}.ti',REF_GENOME=REF_GENOME),
		expand('reference-index/{REF_GENOME}/{REF_GENOME}.seq',REF_GENOME=REF_GENOME),
		expand('reference-index/{REF_GENOME}/{REF_GENOME}.chrlist',REF_GENOME=REF_GENOME),
		expand('reference-index/{REF_GENOME}/SAindex',REF_GENOME=REF_GENOME),
		expand('transcript-counts/{SAMPLE}_rsem.isoforms.results',SAMPLE=SAMPLES),
		expand('transcript-counts/{SAMPLE}_rsem.genes.results',SAMPLE=SAMPLES),		
		#expand('gene-counts/{SAMPLE}_HTSeq_union_gff3_no_gene_ID.log',SAMPLE=SAMPLES),
		#expand('gene-counts/{SAMPLE}_HTSeq.csv',SAMPLE=SAMPLES)


rule clean:
	input:
		r1="raw-reads/{SAMPLE}_read_1.fastq.gz",
		r2="raw-reads/{SAMPLE}_read_2.fastq.gz"
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
#This sets up the reference geneome that STAR will align to. Only needs to be done once for all subsequent STAR alignment jobs 
#rule star_reference:
#	input:
#		fasta="reference-input/{REF_GENOME}/{REF_GENOME}.fasta",
#		gtf="reference-input/{REF_GENOME}/{REF_GENOME}.gtf"
#	output:
#		'reference-index/{REF_GENOME}/SAindex'
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
                fasta = 'reference-input/{REF_GENOME}/{REF_GENOME}.fasta',
                gtf = 'reference-input/{REF_GENOME}/{REF_GENOME}.gtf'

        params:
                outdir = 'reference-index/{REF_GENOME}/{REF_GENOME}'

        output:
                grp = 'reference-index/{REF_GENOME}/{REF_GENOME}.grp',
                ti = 'reference-index/{REF_GENOME}/{REF_GENOME}.ti',
                seq = 'reference-index/{REF_GENOME}/{REF_GENOME}.seq',
                chrlist = 'reference-index/{REF_GENOME}/{REF_GENOME}.chrlist',
                SAindex = 'reference-index/{REF_GENOME}/SAindex',
		#refind = directory('reference-index/{REF_GENOME}')

        #conda: "environment.yml"

        shell:
                'set +u && '
                'eval "$(conda shell.bash hook)" && '
                'conda activate rna-seq && '
                'set -u && '
                'rsem-prepare-reference {input.fasta} --star -p {threads} --gtf {input.gtf} {params.outdir}'

rule align:
	input:
		r1 = 'clean-reads/{SAMPLE}_read_1_fastp.fastq.gz',
		r2 = 'clean-reads/{SAMPLE}_read_2_fastp.fastq.gz',
		ref_genome = expand('reference-index/{REF_GENOME}',REF_GENOME=REF_GENOME)

	params:
		indir  = expand('reference-index/{REF_GENOME}',REF_GENOME=REF_GENOME),
		tempdir = '$PBS_JOBFS',
		outdir = 'aligned-reads/{SAMPLE}',
	output:
		temp('aligned-reads/{SAMPLE}_single_pass/Aligned.sortedByCoord.out.bam'),
		temp('aligned-reads/{SAMPLE}_single_pass/Aligned.toTranscriptome.out.bam')
	#conda: "environment.yml"
	shell:
		'set +u && '
		'eval "$(conda shell.bash hook)" && '
		'conda activate rna-seq && '
		'set -u && '
		#'rm -rf {params.outdir} && '
		#'mkdir {params.outdir} && '
		#'cd {params.outdir} && '
		'STAR --runThreadN {threads} '
		'--genomeDir {params.indir} '
		'--outFileNamePrefix {params.tempdir}/ '
		'--readFilesIn {input.r1} {input.r2} '
		'--readFilesCommand zcat '
		'--outSAMtype BAM SortedByCoordinate '
		'--quantMode TranscriptomeSAM GeneCounts && '
		'mv {params.tempdir}/* {params.outdir}'

rule align_pass_1:
	input:
		r1 = 'clean-reads/{SAMPLE}_read_1_fastp.fastq.gz',
		r2 = 'clean-reads/{SAMPLE}_read_2_fastp.fastq.gz',
		ref_genome = expand('reference-index/{REF_GENOME}/SAindex',REF_GENOME=REF_GENOME)
	params:
		indir  = expand('reference-index/{REF_GENOME}',REF_GENOME=REF_GENOME),
		tempdir = '$PBS_JOBFS',
		outdir = 'aligned-reads/{SAMPLE}_pass_1',
		rmbam = 'aligned-reads/{SAMPLE}_pass_1/Aligned.out.bam'
	output:
		temp('aligned-reads/{SAMPLE}_pass_1/SJ.out.tab')
	#conda: "environment.yml"
	shell:		
		'set +u && '
		'eval "$(conda shell.bash hook)" && '
		'conda activate rna-seq && '
		'set -u && '
		#rm -rf {params.outdir} && \
		#mkdir {params.outdir} && \
		#cd {params.outdir} && \
		'STAR --runThreadN {threads} '
		'--outFileNamePrefix {params.tempdir}/ '
		'--genomeDir {params.indir} '
		'--readFilesIn {input.r1} {input.r2} '
		'--readFilesCommand zcat '
		'--outSAMtype BAM Unsorted && rm {params.rmbam} &&'
		'mv {params.tempdir}/* {params.outdir}' 

rule filter:
	input:
		'aligned-reads/{SAMPLE}_pass_1/SJ.out.tab',
	output:
		temp('splice-junctions/{SAMPLE}_pass_1_SJ.filtered.tab')
	#conda: "environment.yml"
	shell:
		'set +u && '
		'eval "$(conda shell.bash hook)" && '
		'conda activate rna-seq && '
		'set -u && '
		'awk "{{if (\$7 >= 3) print \$0 }}" {input[0]} > {input[0]}.filtered && '
		'mv {input[0]}.filtered {output}'

rule align_pass_2:
	input:
		r1 = 'clean-reads/{SAMPLE}_read_1_fastp.fastq.gz',
		r2 = 'clean-reads/{SAMPLE}_read_2_fastp.fastq.gz',
		sj_files =expand('splice-junctions/{SAMPLE}_pass_1_SJ.filtered.tab',SAMPLE=SAMPLES),
		ref_genome = expand('reference-index/{REF_GENOME}/SAindex',REF_GENOME=REF_GENOME)

	params:
		indir  = expand('reference-index/{REF_GENOME}',REF_GENOME=REF_GENOME),
		tempdir = '$PBS_JOBFS',
		outdir = 'aligned-reads/{SAMPLE}_pass_2'
	output:
		temp('aligned-reads/{SAMPLE}_pass_2/Aligned.sortedByCoord.out.bam'),
		temp('aligned-reads/{SAMPLE}_pass_2/Aligned.toTranscriptome.out.bam')
	#conda: "environment.yml"
	shell:
		'set +u && '
		'eval "$(conda shell.bash hook)" && '
		'conda activate rna-seq && '
		'set -u && '
		#'rm -rf {params.outdir} && '
		#'mkdir {params.outdir} && '
		#'cd {params.outdir} && '
		'STAR --runThreadN {threads} '
		'--genomeDir {params.indir} '
		'--outFileNamePrefix {params.tempdir}/ '
		'--readFilesIn {input.r1} {input.r2} '
		'--readFilesCommand zcat '
		'--outSAMtype BAM SortedByCoordinate '
		'--sjdbFileChrStartEnd {input.sj_files} '
		'--quantMode TranscriptomeSAM GeneCounts && '
		'mv {params.tempdir}/* {params.outdir}' 


rule rsem:
	input:
		bam = 'aligned-reads/{SAMPLE}_pass_2/Aligned.toTranscriptome.out.bam',
		ref = expand('reference-index/{REF_GENOME}/{REF_GENOME}.{ext}',REF_GENOME=REF_GENOME,ext=['grp','ti','seq','chrlist'])

	params: 
		indir  = expand('reference-index/{REF_GENOME}/{REF_GENOME}',REF_GENOME=REF_GENOME), 
		outdir = 'transcript-counts/{SAMPLE}_rsem'

	output: 'transcript-counts/{SAMPLE}_rsem.isoforms.results',
		'transcript-counts/{SAMPLE}_rsem.genes.results'
	#log:	'logs/{SAMPLE}_rsem.log'
	#conda: "environment.yml"

	shell:
		'set +u && '
		'eval "$(conda shell.bash hook)" && '
		'conda activate rna-seq && '
		'set -u && '
		'rsem-calculate-expression --bam --paired-end -p 1 {input.bam} {params.indir} {params.outdir}'
		

#Use htseq rule to quantify GENE counts.  
rule htseq:
	input:
		bam = 'aligned-reads/{SAMPLE}_pass_2/Aligned.sortedByCoord.out.bam',
		gff = expand('reference-input/{REF_GENOME}/{REF_GENOME}.gff3',REF_GENOME=REF_GENOME)
	output:
		'gene-counts/{SAMPLE}_HTSeq_union_gff3_no_gene_ID.log',
		'gene-counts/{SAMPLE}_HTSeq.csv'
	#conda: "environment.yml"
	shell:
		'set +u && '
		'eval "$(conda shell.bash hook)" && '
		'conda activate rna-seq && '
		'set -u && '
		'htseq-count -m union -s no -t gene -i ID -r pos -f bam {input.bam} {input.gff} &> {output[0]} && '
		'grep {GENE_CODE} {output[0]} | sed "s/gene://g" > {output[1]}'
