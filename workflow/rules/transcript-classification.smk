import math
def batches(file_path,transcripts_per_batch):
	num_lines = sum(1 for line in open(file_path))
	bat = math.ceil(num_lines/(transcripts_per_batch*2))
	return(bat)
#TODO fix blast process so that it is generalisable to any dataset. This might involve using subworkflow concept which produces the transcriptome before blast etc are run (could be messy).
bat=batches(trinitydir+'Trinity.fasta',100000)
rule blast_prep:
# remove transcripts below certain length, then split transcriptome up into 
	input:
		transcriptome=trinitydir+'Trinity.fasta'
	params:
		outdir="blastx",
		transcripts_per_batch=100000,
		min_length=150
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
