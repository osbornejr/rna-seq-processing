rule blastx:
	input:
		query=trinitydir+'Trinity.test.fasta',
		subject=input_path+config["blastx_subject"]
	output:
		db_location='blastx/uniprot_db',
		blastx='blastx/output.blast.txt'
	log: 	
		'logs/transcript-classification/blastx.out'
	shell:
		'makeblastdb '
		'-in {input.subject} '
		'-out {output.db_location} '
		'-dbtype prot '
		'-parse_seqids && '
		'blastx '
		'-query {input.query} '
		'-db {output.db_location} '
		'-num_threads 16 '
		'-out {output.blastx} >> {log}'
