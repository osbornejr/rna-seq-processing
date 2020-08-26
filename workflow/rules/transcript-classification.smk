rule blastx:
	input:
		query=trinitydir+'Trinity.fasta',
		subject=input_path+config["blastx_subject"]
	output:
		'blastx/output.blast.txt'
	log: 	
		'logs/transcript-classification/blastx.out'
	shell:
		'blastx '
		'-query {input.query} '
		'-subject {input.subject} '
		'-num_threads {threads} '
		'-out {output} >> {log}'
