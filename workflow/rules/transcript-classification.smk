rule blastx:
	input:
		query=trinitydir+'Trinity.fasta',
		subject=config["path_to_files"]+'/../uniprot_sprot.fasta'
	output:
		'blastx/output.blast.txt'
	log: 	
		'logs/transcript-classification/blastx.out'
	shell:
		'blastx '
		'-query {input.query} '
		'-subject {input.subject} '
		'-out {output} >> {log}'
