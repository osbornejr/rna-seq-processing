##make blast core_nt database on aws instance
### you can also get the database directly from the NCBI ftp  using ``update_blastdb.pl --decompress nt``. note that core_nt is a condensed version of nt that ncbi now recommends, and it is approx 1/3 of the size. Either way, a huge download, so you need a lot of space (including to decompress, where it gets even larger). If operating on AWS, the s3 bucket is much preferable; not only does it download a lot quicker than the glacial ncbi server, but it comes already decompressed. you can also get the full nt on the s3 bucket, or other subsets of it if required.    
mkdir blast_db
aws s3 cp --no-sign-request s3://ncbi-blast-databases/2025-09-11-01-05-02/ blast_db/ --exclude "*" --include "core_nt*" --recursive
##run on large instance with multitrheading
#NOTE: blast is not great at providing documentation or feedback as it goes, but it can be run very efficiently if given the right resources. How do you know what's right though? See below...
#some number of threads is good (more the merrier, it will use them all full tilt, but there is some single threaded downtime between each batch that it runs, so don't pay for too many, diminishing returns ;) ) but the main thing is to ensure there is enough RAM to allow the database you are blasting to be cached.
#You find this somewhat cryptically using the core_nt.nal (or *.nal for whatever db you use), which is just a json file with an entry "bytes-to-cache=". 
## In the case of core_nt, this is currently approx 235GB, so a high memory instance with 256GB RAM is perfect-- a bit of head room is required for other processes on top of the shared cache.
blastn -num_threads 32 -query query.fasta -db /home/ec2-user/rna-seq-processing/blast_db/core_nt -out results.txt -task megablast

##after running, we want to condense these down to one match (or none if no hits) per transcript. there are python scripts to do this. because there are several potential options on how to select, we grab each of
##choose best Cicer arietinum non-predicted, if not best cicer arietinum, if not top match of any species
python3 scripts/filter_blastn.py results.txt filtered_results.txt
##choose best Cicer arietinum, if not top match of any species
python3 scripts/filter_blastn_no_nonpredict.py results.txt filtered_results_no_nonpredict.txt
##just top match of any species
python3 scripts/filter_blastn_just_top_any_species.py results.txt filtered_results_just_top_any_species.txt

##typically the first two more selective ones will fail on some small number of transcripts. this could be debugged, but easiest to do so on subset so for now we form small fasta files for those missing 
#for each filtered_*.txt file:
cut -d ',' -f 1 filtered_output_just_top_any_species.txt > filtered_names_just_top_any_species.txt
## (this command needs full transcriptome names derived from Trinity.fasta
comm -23 <(sort trinity_names.txt) <(sort filtered_names.txt) > missing.txt
while IFS= read -r p; do rg -A1 --no-context-separator "^>${p}\b" Trinity.fasta; done < missing_ref-seq.txt >missing_ref-seq.fasta
