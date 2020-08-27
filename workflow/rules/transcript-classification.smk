rule blastx:
	input:
		query=trinitydir+'Trinity.test.fasta',
		subject=input_path+config["blastx_subject"]
	params:
		db_location='blastx/uniprot_db'
	output:
		blastx='blastx/output.blast.txt'
	log: 	
		'logs/transcript-classification/blastx.out'
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
		'-out {output.blastx} >> {log}'
