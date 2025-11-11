#Transcriptome assembly using reference genome
###2025 note: For reference genome runs, very large machines are not necessary for any rules. 8 or 16 CPUS with ~32GB is what each needs at max capacity, but there are still large periods when it will run on single thread (see rsem rule note). The whole thing for all 24 chickpea samples can be run on a machine the size above in a day or so.
## if a faster turn around is needed, then running each job on its own node of that size (i.e. cluster style) is MUCH more efficient than one large AWS style instance, which will only use 5-20% of CPU at any one time.

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

##2025 update-- have left the single pass align as is with temp dir in case it is needed again (eg on gadi) other rules changed to work without temp dir on aws. note all rules have threads manually set as parameters now rather than by cluster config... in long run this would be set now by snakemake profiles.
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
	output:##2025 note: do we want this to be discarded automatically? Depends on directory size(s)
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

##2025 note: if it is ever necessary, this could be optimised much more for a HPC by running each constituent part separately.
# eg. for a given input alignment:
#FIRST COMMAND (only runs on a single process, can be looong dependning on size of alignment)
#rsem-parse-alignments reference-index/<chosen-ref>/<chosen-ref> transcript-counts/polyA_plus_sample_X/polyA-plus_X.temp/polyA_plus_sample_5 transcript-counts/polyA_plus_sample_X/polyA-plus_X.stat/polyA_plus_sample_5
#       aligned-reads/polyA_plus_sample_5/polyA-plus_5_pass_2/Aligned.toTranscriptome.out.bam 3 -tag XM
#SECOND COMMAND (also single process, pretty quick though)
#rsem-build-read-index 32 1 0 transcript-counts/polyA_plus_sample_8/polyA_plus_sample_8.temp/polyA_plus_sample_8_alignable_1.fq transcript-counts/polyA_plus_sample_8/polyA_plus_sample_8.temp/polyA_plus_sample_8_alignable_2.fq
#THIRD COMMAND (runs in parallel for most of run, may trigger gibbs automatically?) 
#rsem-run-em reference-index/CDCFrontier_v2.0/CDCFrontier_v2.0 3 transcript-counts/polyA_plus_sample_1/polyA_plus_sample_1 transcript-counts/polyA_plus_sample_1/polyA_plus_sample_1.temp/polyA_plus_sample_1 transcript-counts/polyA_plus_sample_1/polyA_plus_sample_1.stat/polyA_plus_sample_1 -p 16 -b aligned-reads/polyA_plus_sample_1/polyA_plus_sample_1_pass_2/Aligned.toTranscriptome.out.bam 0 --gibbs-out
#FOURTH COMMAND (also in parallel)
#rsem-run-gibbs reference-index/CDCFrontier_v2.0/CDCFrontier_v2.0 transcript-counts/polyA_minus_sample_3/polyA_minus_sample_3.temp/polyA_minus_sample_3 transcript-counts/polyA_minus_sample_3/polyA_minus_sample_3.stat/polyA_minus_sample_3 200 1000 1 -p 16
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
