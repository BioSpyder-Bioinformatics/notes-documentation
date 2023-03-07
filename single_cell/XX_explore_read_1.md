# Explore read 1

It can happen that the read 2 have much lower read quality than read 1! This only leaves you with the ability to analyse read 1s.


Align everything on BWA for read 1
bwa from conda bin!!
`bwa mem -v 1 -c 2 -L 100 {manifest file for BWA.fa} {file of read 1s.fastq.gz} > alignment.sam`
The flags are essential to keep a strict alignment and not break the reads across the genome


The obtained alignment.sam can be analysed with samtools


### Samtools view

- Convert sam to bam (takes less space)
	+ `samtools view -bS alignment.sam > alignment.bam`
- Obtain ALL aligned reads
	+ `samtools view -F 4 alignment.sam > mapped.bam`
- Obtain all NOT aligned reads
	+ `samtools view -f 4 alignment.sam > unmapped.bam`


### Search for motif
Basically we want to have a look for a recurrent motif on the read *2*, so to make sure that there was a certain level of reproducibility in the sequencing itself.

On website meme-suite.org, using the tool Meme Motif Search:

- Prerequisite
	+ Got a fasta file by taking the aligned reads with Samtools, getting the head, and searching the identifier in the read 2 files
	+ Capture the identifier and part of the read using grep (first line only)
	+ Append all this to file

- On motif page, upload the fasta file (or copy paste text) and start the search

