shell.executable
("/bin/bash")
SAMPLES=expand("{TYPE}/{TYPE}_sample_{N}/{TYPE}_sample_{N}",TYPE=["polyA+","polyA-"],N=[str(n) for n in range(13)][1:])
IDS=expand("{TYPE}_sample_{N}",TYPE=["polyA+","polyA-"],N=[str(n) for n in range(13)][1:])
REF_GENOME="Kabuli_UWA-v2.6.3"
GENE_CODE="Ca"

trinitydir="trinity-transcriptome-assembly/"
basedir=workflow.basedir+'/'
#TODO put long shell commands into scripts and source?
#TODO split Snakefile into two separate sub workflows, de novo and reference assembly 
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
	
	log: 	
		basedir+"logs/trinity/trinity_assembly_normalisation.out"	
	run:
                left=",".join(map(str,input.left))
                right=",".join(map(str,input.right))
                shell(
		'mkdir -p normalised-reads >> {log} && '
		'cd normalised-reads >> {log} && '
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
		'--PARALLEL_STATS >> {log}')


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
	
rule trinity_assembly_phase_1:
#This rule runs Trinity phase 1 using the storage on a highmemory node in Gadi. This allows unlimited numbers of files to be created  (and could theoretically mean running without any normalisation?). The resulting read partitions are zipped up and re-stored locally, along with the recursive commands required in phase 2. (On other systems this method might not be necessary depending on file count quotas. In any case, the below rule is currently context specific to a PBS/Gadi/NCI framework) TODO remove this reliance?  
	input:
		left=basedir+"normalised-reads/left.norm.fq",
		right=basedir+"normalised-reads/right.norm.fq"
	
	params:	
		tempdir="$PBS_JOBFS/"

	output:
		trinitydir+"phase_1.tar.gz"
	
	log: 	
		basedir+"logs/trinity/trinity_assembly_phase_1.out"	
	run:			
		shell(
		#run trinity in grid mode from root of temp drive
		'cd {params.tempdir} >> {log} && '
		'Trinity '
		'--no_distributed_trinity_exec '
		'--seqType fq '
		'--left {input.left} '
		'--right {input.right} '
		'--CPU 16 '
		'--max_memory 950G '
		'--no_normalize_reads '
		'--output '+trinitydir+'  >> {log} && ' 
		'sed -i "s~{params.tempdir}'+trinitydir+'~PHASE_2_PREFIX~g" '+trinitydir+'recursive_trinity.cmds >> {log} && '
		'head '+trinitydir+'recursive_trinity.cmds >> {log} && '
		'tar -czf '+basedir+trinitydir+'phase_1.tar.gz '+trinitydir+'recursive_trinity.cmds '+trinitydir+'/read_partitions >> {log} ' )  



rule trinity_assembly_phase_2:
#phase 2 is executed outside of Trinity by simply running the recursive_trinity.cmds in parallel. 
#This requires a large number of CPUs, as well as the storage of a very high file count.
# Running these commands separately on 1000s of CPUs allows us to fully exploit the parallelisation speed up offered at this step, and is essentially required to keep the job under the 48 hour time limit.
# The phase 1 tarball is unzipped onto the primary node (Node 0), and 48 cores of this node is dedicated to "housekeeping".
# They are responsible for bouncing down the required input file for each command to /scratch/, and retrieving the output upon command completion to bring it back to Node 0.
# All other cores are dedicated to running the commands within /scratch/ space.
# This procedure is necessary because compute cores can only access storage on their own node and on /scratch/.
# The bounce down and up is coupled with a bit of cleaning by the housekeepers to make sure that the file count does not get too large on scratch.
# Once all commands are complete, the node directory is tarballed to /scratch/ (precious cargo!) and then the Trinity aggregate commands are run on the outputs.
# There may be further commands to run once the annotation process is clearer (this is why we tarball the directory).
# Outputs of this step will be the phase 2 tarball along with the final Trinity.fasta transcriptome file.
#TODO make separate path/switch that allows the non- parallel/slight parallel version to run if so desired (will require phase 1 tarball to give full trinity directory. (i.e. SLOW)
	input:
		trinitydir+"phase_1.tar.gz",
		left=basedir+"normalised-reads/left.norm.fq",
		right=basedir+"normalised-reads/right.norm.fq"
	
	params:
		tempdir="$PBS_JOBFS/",
		basedir=basedir,
		n_cpus="$PBS_NCPUS",
		threads_per_node=48
	
	output:
		trinitydir+"phase_2.tar.gz"
	log: 	
		basedir+"logs/trinity/trinity_assembly_phase_2.out"	
	run: 
		shell(
		#unzip phase 1 to housekeeping node and change prefix to suit bouncedown
		'tar -xzf '+trinitydir+'phase_1.tar.gz -C {params.tempdir} >> {log} && '
		'sed -i "s~PHASE_2_PREFIX~'+trinitydir+'~g" {params.tempdir}'+trinitydir+'recursive_trinity.cmds && '
		'head {params.tempdir}'+trinitydir+'recursive_trinity.cmds >> {log} && '
		#write bounce commands to be run on each cmd. The last CPU will always belong to Node 0, where phase 1 should be extracted to. This node is used as a housekeeper to move inputs and outputs between Node 0 and the scratch drive. All other CPUS will be assigned to 			process a specific command from recursive_trinity.cmds. TODO may need more than one housekeeper?
		'trinity_bounce() {{ '
		'basedir='+basedir+'; '
		'mkdir -p $6 &>> {log}; mv {params.tempdir}$2 $6 &>> {log}; '
		'pbsdsh -n $(($7+{params.threads_per_node}-1)) -- bash -l -c "cd $basedir;$1 $2 $3 $4 $5 &>> {log};" &>> {log} ; '
		'mv $4.Trinity.fasta {params.tempdir}$6;rm $2 &>> {log};rmdir -p --ignore-fail-on-non-empty $6 &>> {log}; }} && '
		'export -f trinity_bounce && '
		# remove any relics of past runs
		'rm -rf '+trinitydir+'read_partitions && '
		#run GNU parallel to distribute and bounce cmds to available nodes. Note that parallel is given one node less than the total (Node 0 will be for housekeeping, as prescribed in trinity_bounce
		"""cat {params.tempdir}"""+trinitydir+"""recursive_trinity.cmds | parallel --colsep '"' --env trinity_bounce -j $(({params.n_cpus}-{params.threads_per_node})) trinity_bounce {{1}} {{2}} {{3}} {{4}} {{5}} {{2//}} {{%}}  && """
		#save job directory to zip (very precious; must be run from job directory)
		'cd {params.tempdir} && '
		'tar -cvzf '+basedir+trinitydir+'phase_2.tar.gz '+trinitydir+' >> {log} ')  

rule trinity_assembly_finalise:
	input:
		trinitydir+"phase_2.tar.gz"
	
	params:
		tempdir="$PBS_JOBFS/",
		basedir=basedir
	
	output:
		#trinitydir+"Trinity.fasta"
	log: 	
		basedir+"logs/trinity/trinity_assembly_finalise.out"	
	run: 
		shell(
		#unzip phase 2 to job node for final collation steps
		'tar -xzf '+trinitydir+'phase_2.tar.gz -C {params.tempdir} >> {log} && '
		#aggregate found reads to one transcriptome TODO remove reference to specific version of Trinity TODO need more than just fasta file to annotate?
		'find {params.tempdir}'+trinitydir+'read_partitions/ -name "*Trinity.fasta" | $CONDA_PREFIX/opt/trinity-2.9.1/util/support_scripts/partitioned_trinity_aggregator.pl --token_prefix TRINITY_DN --output_prefix Trinity >'+trinitydir+'Trinity.fasta')

rule trinity_assembly_phase_1_full:
#This rule runs Trinity phase 1 using the storage on a highmemory node in Gadi. This allows unlimited numbers of files to be created  (and could theoretically mean running without any normalisation?). The resulting trinity directory is tarballed and required in phase 2 full. (On other systems this method might not be necessary depending on file count quotas. In any case, the below rule is currently context specific to a PBS/Gadi/NCI framework) TODO remove this reliance?  
	input:
		left=basedir+"normalised-reads/left.norm.fq",
		right=basedir+"normalised-reads/right.norm.fq"
	
	params:	
		tempdir="$PBS_JOBFS/"

	output:
		trinitydir+"phase_1_full.tar.gz"
	
	log: 	
		basedir+"logs/trinity/trinity_assembly_phase_1_full.out"	
	run:			
		shell(
		'cd {params.tempdir} >> {log} && '
		'Trinity '
		'--no_distributed_trinity_exec '
		'--seqType fq '
		'--left {input.left} '
		'--right {input.right} '
		'--CPU 16 '
		'--max_memory 950G '
		'--no_normalize_reads '
		'--output '+trinitydir+'  >> {log} && ' 
		'sed -i "s~{params.tempdir}'+trinitydir+'~PHASE_2_PREFIX~g" '+trinitydir+'recursive_trinity.cmds >> {log} && '
		'head '+trinitydir+'recursive_trinity.cmds >> {log} && '
		'tar -czf '+basedir+trinitydir+'phase_1_full.tar.gz '+trinitydir+' >> {log} ' )  

rule trinity_assembly_phase_2_full:
	input:
		trinitydir+"phase_1_full.tar.gz",
		left=basedir+"normalised-reads/left.norm.fq",
		right=basedir+"normalised-reads/right.norm.fq"
	
	params:
		tempdir="$PBS_JOBFS/",
		basedir=basedir,
		n_cpus="$PBS_NCPUS",
		threads_per_node=48
	
	output:
		trinitydir+"phase_2_full.tar.gz",
		trinitydir+"Trinity.fasta"
	log: 	
		basedir+"logs/trinity/trinity_assembly_phase_2_full.out"	
	run: 
		shell(
		#unzip phase 1 full to job node and change prefix to suit node
		'tar -xzf '+trinitydir+'phase_1.tar.gz -C {params.tempdir} >> {log} && '
		'sed -i "s~PHASE_2_PREFIX~'+trinitydir+'~g" {params.tempdir}'+trinitydir+'recursive_trinity.cmds && '
		'head {params.tempdir}'+trinitydir+'recursive_trinity.cmds >> {log} && '
		#run trinity in grid mode from job node
		'cd {params.tempdir} && '
		'Trinity '
		"""--grid_exec "parallel -j {params.n_cpus} pbsdsh -n {{%}} -- bash -l -c '{{}}'< " """
		'--seqType fq '
		'--left {input.left} '
		'--right {input.right} '
		'--CPU 16 '
		'--max_memory 10G '
		'--no_normalize_reads '
		'--output '+trinitydir+' >> {log} && '
		#move final transcriptome to base
		'mv '+trinitydir+'Trinity.fasta '+basedir+trinitydir+' >> {log} && '
		'tar -cvzf '+basedir+trinitydir+'phase_2_full.tar.gz '+trinitydir+' >> {log}') 
		
	
rule clean:
# all in one cleaning of fastq files with fastp. Could alternatively do this with trimmomatic? But need specific adapter sequences etc. TODO double check quality independently after running fastp with fastqc
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
