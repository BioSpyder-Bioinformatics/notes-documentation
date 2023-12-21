# How to make a new local blast genome
First of all find the genome fasta in Ensembl db
Download the zipped fasta. Most likely you'll need the DNA toplevel fasta

## Create DB
makeblastdb -in ./dna_fasta/Oar_v4.0.dna.toplevel.fa -dbtype nucl -parse_seqids -out Oar_v4_blastdb -title Oar_v4_blastdb
 


## Run BLAST on this genome
blastn -query ./probe_fasta/sheep_udo.fa  -db Oar_v4_blastdb -out sheep_do_blast_out  -outfmt 6 -task blastn-short



### Read in with pandas
blast_df=pd.read_csv(blastoutput,sep='\t', header=None, names=['QUERY','TARGET','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'])



