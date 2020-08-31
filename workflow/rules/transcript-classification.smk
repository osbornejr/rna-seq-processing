import math
def batches(file_path,transcripts_per_batch):
	num_lines = sum(1 for line in open(file_path))
	bat = math.ceil(num_lines/(transcripts_per_batch*2))
	return(bat)
#TODO fix blast process so that it is generalisable to any dataset. This might involve using subworkflow concept which produces the transcriptome before blast etc are run (could be messy).
bat=batches(trinitydir+'Trinity.fasta',100000)
rule pre_blastx:
# remove transcripts below certain length, then split transcriptome up into 
	input:
		transcriptome=trinitydir+'Trinity.fasta'
	params:
		outdir="blastx",
		transcripts_per_batch=100000,
		min_length=config["min_transcript_length"]
	output:
		expand('blastx/blast_input_{X}.fasta',X=list(range(1+bat))[1:bat+1])
	shell:
		'''awk 'BEGIN {{FS = "[ ,=]" ;  RS = ">" ; ORS = ""}} $3 >= {params.min_length} {{print ">"$0}}' '''
		'''{input.transcriptome} > {params.outdir}/Trinity{params.min_length}.fasta && '''
		'bash '+flowdir+'scripts/split_transcriptome_for_blast.sh {params.outdir}/Trinity{params.min_length}.fasta {params.outdir} {params.transcripts_per_batch}'

BAT=list(range(1+bat))[1:bat+1]
rule blastx:
#TODO maybe switch if possible to running each batch on the same node. Might require some form of MPI parallelism?
	input:
		query='blastx/blast_input_{batch}.fasta',
		subject=input_path+config["blastx_subject"]
	params:
		db_location='blastx/uniprot_db'
	output: 
		'blastx/output_{batch}.blast.txt'
	log: 	
		'logs/transcript-classification/blastx_{batch}.out'
	shell:
		'makeblastdb '
		'-in {input.subject} '
		'-out {params.db_location} '
		'-dbtype prot '
		'-parse_seqids && '
		'blastx '
		'-query {input.query} '
		'-db {params.db_location} '
		'-num_threads 48 '
		'-out {output} >> {log}'

rule post_blastx:
	input:
		expand('blastx/output_{batch}.blast.txt',batch=list(range(1+bat))[1:bat+1]),	
		'blastx/Trinity'+config["min_transcript_length"]+'.fasta'
	params:
		min_length=config["min_transcript_length"]
	output:
		'blastx/blast-coding.fasta',
		'blastx/blast-non-coding.fasta'
	shell:
		'''rm -f blastx/blast-no-hit.list && '''	
		'''for file in blastx/output*; do awk 'BEGIN {{RS = "Query=";ORS=""}} /***** No hits found *****/ {{print ">"$1"\\n"}} ' $PWD/$file >> blastx/blast-no-hit.list ;done &&'''
		'''grep -A1 -Ff blastx/blast-no-hit.list blastx/Trinity{params.min_length}.fasta | grep -v '^--' > blastx/blast-non-coding.fasta && '''
		'''grep -vFf blastx/blast-non-coding.fasta blastx/Trinity{params.min_length}.fasta > blastx/blast-coding.fasta'''	
rule rnasamba_classify:
	input:
		noncoding='blastx/blast-coding.fasta',
		coding='blastx/blast-non-coding.fasta',
		training_set=input_path+config["RNAsamba_training_set"]
	output:
		noncoding='rnasamba/non-code-hits.list',
		coding='rnasamba/code-hits.list'
	shell:
		'rnasamba classify rnasamba/non-code-classification.tsv {input.noncoding} {input.training_set} && '
		'rnasamba classify rnasamba/code-classification.tsv {input.coding} {input.training_set} && '
		'''awk '/noncoding/ {{print $1}}' non-code-classification.tsv > {output.noncoding} && '''
		'''awk '/ coding/ {{print $1}}' code-classification.tsv > {output.coding} '''
