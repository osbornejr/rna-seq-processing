shell.executable
("/bin/bash")
SAMPLES=expand("{TYPE}/{TYPE}_sample_{N}/{TYPE}_sample_{N}",TYPE=["polyA+","polyA-"],N=[str(n) for n in range(13)][1:])
IDS=expand("{TYPE}_sample_{N}",TYPE=["polyA+","polyA-"],N=[str(n) for n in range(13)][1:])
REF_GENOME="Kabuli_UWA-v2.6.3"
GENE_CODE="Ca"

trinitydir="trinity-transcriptome-assembly/"
basedir=workflow.basedir+'/'
#Transcriptome assembly using reference genome
#rule all:
#	input:
#		expand("clean-reads/{SAMPLE}_read_{N}_fastp.fastq.gz",SAMPLE=SAMPLES,N=["1","2"]),
#		expand("reports/{SAMPLE}_fastp.{EXT}",SAMPLE=SAMPLES,EXT=["html","json"]),
#		#expand('reference-index/{REF_GENOME}_star',REF_GENOME=REF_GENOME),
#		expand('aligned-reads/{SAMPLE}_pass_1/SJ.out.tab',SAMPLE=SAMPLES),
#		expand('splice-junctions/{SAMPLE}_pass_1_SJ.filtered.tab',SAMPLE=SAMPLES),
#		expand('aligned-reads/{SAMPLE}_pass_2/Aligned.sortedByCoord.out.bam',SAMPLE=SAMPLES),
#		expand('aligned-reads/{SAMPLE}_pass_2/Aligned.toTranscriptome.out.bam',SAMPLE=SAMPLES),
#		expand('reference-index/{REF_GENOME}/{REF_GENOME}.grp',REF_GENOME=REF_GENOME),
#		expand('reference-index/{REF_GENOME}/{REF_GENOME}.ti',REF_GENOME=REF_GENOME),
#		expand('reference-index/{REF_GENOME}/{REF_GENOME}.seq',REF_GENOME=REF_GENOME),
#		expand('reference-index/{REF_GENOME}/{REF_GENOME}.chrlist',REF_GENOME=REF_GENOME),
#		expand('reference-index/{REF_GENOME}/SAindex',REF_GENOME=REF_GENOME),
#		expand('transcript-counts/{SAMPLE}_rsem.isoforms.results',SAMPLE=SAMPLES),
#		expand('transcript-counts/{SAMPLE}_rsem.genes.results',SAMPLE=SAMPLES),		
#		#expand('gene-counts/{SAMPLE}_HTSeq_union_gff3_no_gene_ID.log',SAMPLE=SAMPLES),
#		#expand('gene-counts/{SAMPLE}_HTSeq.csv',SAMPLE=SAMPLES)

#De novo transcriptome assembly
rule all:
	input:					
		expand("clean-reads/{SAMPLE}_read_{N}_fastp.fastq.gz",SAMPLE=SAMPLES,N=["1","2"]),
		expand("reports/{SAMPLE}_fastp.{EXT}",SAMPLE=SAMPLES,EXT=["html","json"]),
		basedir+"normalised-reads/left.norm.fq",
		basedir+"normalised-reads/right.norm.fq",
		trinitydir+"Trinity.fasta"
		
rule trinity_normalisation:
	input: 
		left=expand(basedir+"clean-reads/{SAMPLE}_read_1_fastp.fastq.gz",SAMPLE=SAMPLES),
                right=expand(basedir+"clean-reads/{SAMPLE}_read_2_fastp.fastq.gz",SAMPLE=SAMPLES)
	output: 
		left=basedir+"normalised-reads/left.norm.fq",
		right=basedir+"normalised-reads/right.norm.fq"
	
	run:
                left=",".join(map(str,input.left))
                right=",".join(map(str,input.right))
                shell(
                'set +u && '
                'eval "$(conda shell.bash hook)" && '
                'conda activate rna-seq && '
                'set -u && '
		'mkdir -p normalised-reads && '
		'cd normalised-reads && '
		'insilico_read_normalization.pl '
		'--seqType fq '
		'--JM 100G '
		'--CPU 16 '
		'--min_cov 1 '
		'--max_cov 200 '
		'--max_CV 10000 '
		'--left {left} '
		'--right {right} '
		'--pairs_together '
		'--PARALLEL_STATS > ../logs/trinity/trinity_norm.out')


##This method works, but is slow because the grid execution will not work when using temporary node storage.
#rule trinity_assembly:
#	input:
#		left=basedir+"/normalised-reads/left.norm.fq",
#		right=basedir+"/normalised-reads/right.norm.fq"
#
#	params:
#		tempdir="$PBS_JOBFS"
#	
#	output:
#		trinitydir+"/Trinity.fasta"
#	
#	shell:
#		'set +u && '
#                'eval "$(conda shell.bash hook)" && '
#                'conda activate rna-seq && '
#                'set -u && '
#		#run trinity in grid mode from root of temp drive
#		'cd {params.tempdir} && '
#		'Trinity '
#		#'--grid_exec "'+basedir+'/HpcGridRunner/hpc_cmds_GridRunner.pl --grid_conf '+basedir+'/HpcGridRunner/hpc_conf/nci_trinity.conf -c" '
#		'--seqType fq '
#		'--left {input.left} '
#		'--right {input.right} '
#		'--CPU 16 '
#		'--max_memory 950G '
#		'--no_normalize_reads '
#		'--output '+trinitydir+' > '+basedir+'/logs/trinity/trinity.out; '
#		# rescue output files from temp drive before exit
#		'mv farmit* '+basedir+' && '
#		'mv '+trinitydir+'/read_partitions/Fb_0/CBin_41 '+basedir+' && '   
#		'mv {params.tempdir}/'+trinitydir+'/Trinity.fasta {params.tempdir}/'+trinitydir+'/Trinity.fasta.gene_trans_map '+basedir+trinitydir+'/ '
	
#This rule runs Trinity phase 1 using the storage on a highmemory node in Gadi. This allows unlimited numbers of files to be created  (and could theoretically mean running without any normalisation?). The resulting read partitions are zipped up and re-stored locally, along with the recursive commands required in phase 2. (On other systems this method might not be necessary depending on file count quotas. In any case, the below rule is currently context specific to a PBS/Gadi/NCI framework) TODO remove this reliance?  
rule trinity_assembly_phase_1:
	input:
		left=basedir+"normalised-reads/left.norm.fq",
		right=basedir+"normalised-reads/right.norm.fq"
	
	params:	
		tempdir="$PBS_JOBFS"

	output:
		trinitydir+"phase_1.tar.gz"
	
	run:			
		shell(
                'set +u && '
                'eval "$(conda shell.bash hook)" && '
                'conda activate rna-seq && '
                'set -u && '
		#run trinity in grid mode from root of temp drive
		'cd {params.tempdir} && '
		'Trinity '
		'--no_distributed_trinity_exec '
		'--seqType fq '
		'--left {input.left} '
		'--right {input.right} '
		'--CPU 16 '
		'--max_memory 950G '
		'--no_normalize_reads '
		'--output '+trinitydir+' > '+basedir+'logs/trinity/trinity_phase_1.out && ' 
		"""sed -i 's~"{params.tempdir}"""+trinitydir+""""~PHASE_2_PREFIX~g' """+trinitydir+"""recursive_trinity.cmds && """
		'tar -cvzf '+basedir+trinitydir+'phase_1.tar.gz '+trinitydir)  

#TODO this rule	could either be run outside of trinity, manually pasting together the segment runs at the end, or the parallel command could possibly be passed as a --grid_exec command to trinity phase 2. The latter is neater, but will require all of the phase one .ok files to be passed as well. perhaps as simple as zipping EVERYTHINg up after phase 1.	
rule trinity_assembly_phase_2:
	input:
		trinitydir+"phase_1.tar.gz",
		left=basedir+"normalised-reads/left.norm.fq",
		right=basedir+"normalised-reads/right.norm.fq"
	
	params:
		tempdir="$PBS_JOBFS/",
		n_cpus="$PBS_NCPUS"
	
	output:
		trinitydir+"phase_2.tar.gz",
		trinitydir+"Trinity.fasta"
	
	run: 
		shell(
		'cd {params.tempdir} && '
		'tar -xvzf '+basedir+trinitydir+'phase_1.tar.gz && '
		"""sed -i 's~PHASE_2_PREFIX~"{params.tempdir}"""+trinitydir+""""~g' """+trinitydir+"""recursive_trinity.cmds && """
		#run trinity in grid mode from root of temp drive
		'Trinity '
		"""--grid_exec "parallel -j {params.n_cpus} pbsdsh -n {{%}} -- bash -l -c '{{}}'< " """
		'--seqType fq '
		'--left {input.left} '
		'--right {input.right} '
		'--CPU 16 '
		'--max_memory 10G '
		'--no_normalize_reads '
		'--output '+trinitydir+' > '+basedir+'/logs/trinity/trinity_phase_2.out && '
		'tar -cvzf '+basedir+trinitydir+'phase_2.tar.gz '+trinitydir+' && '
		'cd '+basedir+trinitydir+' && '
		'tar -xvzf phase_2.tar.gz && '
		'mv '+trinitydir+'Trinity.fasta ./') 
		
		
	
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
