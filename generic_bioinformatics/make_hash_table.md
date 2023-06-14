# Make hash table

To make hash table:
- Prepare fasta file with all the sequences of interest
- Using python2 run makeGtf.py (kevin)
    - `python /srv/ssrefs/makeGtf.py input.fa > input.gtf`
- Run STAR to make actual hash table
    - `opt/star/STAR --runMode genomeGenerate --runThreadN 7 --genomeDir ./ --genomeFastaFiles input.fa --genomeSAindexNbases 4 --sjdbGTFfile input.gtf --sjdbGTFfeatureExon exon`