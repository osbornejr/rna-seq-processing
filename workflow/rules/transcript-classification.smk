import math
def batches(file_path,transcripts_per_batch):
	num_lines = sum(1 for line in open(file_path))
	bat = math.ceil(num_lines/(transcripts_per_batch*2))
	return(bat)

bat=batches(trinitydir+'Trinity.fasta',100000)
rule blast_batch:
	input:
		transcriptome=trinitydir+'Trinity.fasta'
	params:
		outdir="blastx",
		transcripts_per_batch=100000
	output:
		expand('blastx/blast_input_{X}.fasta',X=list(range(1+bat))[1:bat+1])
	shell:
		'bash '+flowdir+'scripts/split_transcriptome_for_blast.sh {input.transcriptome} {params.outdir} {params.transcripts_per_batch}'

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
