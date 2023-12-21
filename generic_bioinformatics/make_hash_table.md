# Make hash table

To make hash table:
- Prepare fasta file with all the sequences of interest
- Using python2 run makeGtf.py (kevin)
    - `python /srv/ssrefs/makeGtf.py input.fa > input.gtf`
- Run STAR to make actual hash table
    - `/opt/star/STAR --runMode genomeGenerate --runThreadN 7 --genomeDir ./ --genomeFastaFiles input.fa --genomeSAindexNbases 4 --sjdbGTFfile input.gtf --sjdbGTFfeatureExon exon`


Make sure folder, fa and gtf are named the same. References are in /srv/ssrefs


From htop I understood that the main script that manages the alignments is found in /home/shiny/tsqrScripts/queue_alignment.py (I know, it's ironic that it is a python script ahahah).
With a bit of 'find' and 'grep' piping I've realised that the bams and log files are in /srv/data/rna8/, from there I managed to access the actual alignment logs of RNA8, and the error message was....... 
IOError: [Errno 2] No such file or directory: '/srv/ssrefs/RNA8Refs/pharmaq_probes_new_gtf/pharmaq_probes_new_gtf.fa'
This was because I did not name the fasta file and gtf file the same as the directory name... 😫
After renaming the files the alignment failed again, this was because the user 'shiny' did not have the permissions to write in my reference directory. After fixing the permissions, now everything looks and works just as it should! So this should mean that we can now troubleshoot RNA8/TempoSeqR alignments ssa bit easier.

FAILED ALIGNMENTS ARE IN /data/tsqr/Username