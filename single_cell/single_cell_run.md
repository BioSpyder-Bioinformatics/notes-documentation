###

Start connecting to teamviewer

- Very careful cus you have limited time!!

- They leave the sequencing software on
	+ Check that top right has a V green, that means the run is done
	
- Open filezilla 
	+ Go in folder'output' 
	+ Last folder at the bottom is the last run!
	+ Open connection to kevin
		* Create new folder (BIOS#### - given in email)
	+ Select all in their end
		* Drag and drop to kevin
	+ Done



- Connect to kevin
	+ Go to folder with new data
	+ Perform qc on reads
		* `fastqc BIOS*`
			- This produces html reports



- Run BIOS_single_cell on it 
	+ Get fasta file for BWA index (in this case human WT 2.0)
		* /home/salvo/bwaIndices/TempO-Seq_Human_Whole_Transcriptome_2.0_ManifestV2.fa
	+ Get annotation file (GTF)
		* /srv/ssrefs/homo_sapiens/....2.0_manifestV2
	+ These were copied in home dir










### Get stats on Bar code 1
*script* get_stats_bc1.py
usage:
`python3 get_stats_bc1.py -i ./pair01/,./pair02/ -o bios...pc1.csv -bc1 bc1_seqs.txt`
- i:
	+ Comma separated list of files where the reads are
- o:
	+ Output file in csv
- bc1:
	+ tsv file with Condition - Sequence as header






































