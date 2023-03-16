# This is the script to process fasta files

# !!!FOR EACH FILE

## Make the temporary directory and copy the fasta file in!

## Produce the GTF file (THIS IS ALWAYS THE SAME BUT IT'S required in each folder)
Script makeGtf.py only takes in standard input 1 argument, eg reference genome
`python makeGtf.py ../bwaIndices/TempO-Seq_Mouse_Whole_Transcriptome_1.0.fa > input.gtf`

### Produce the Hash Table (Inside each directory)
Aligner STAR produces it
`STAR --runMode genomeGenerate --runThreadN 7 --genomeDir ./ --genomeFastaFiles $1 --genomeSAindexNbases 4 --sjdbGTFfile $2 --sjdbGTFfeatureExon exon`
Where
- $1 -> reference_genome.fa
- $2 -> Newly produced GTF file


Example:
`STAR --runMode genomeGenerate --runThreadN 7 --genomeDir ./ --genomeFastaFiles ../../bwaIndices/TempO-Seq_Mouse_Whole_Transcriptome_1.0.fa --genomeSAindexNbases 4 --sjdbGTFfile input.gtf --sjdbGTFfeatureExon exon`



## Run the STAR alignment
`STAR --genomeDir $1 --readFilesIn $2 --runThreadN 8 --outSAMtype BAM SortedByCoordinate --scoreDelOpen -10000 --scoreInsOpen -10000 --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMunmapped Within --outFileNamePrefix $3`
Where 
- $1 -> Path to directory where the hash table is 
- $2 -> Path to fastq file
- $3 -> Prefix of output files

Example
`STAR --genomeDir ./ --readFilesIn <(gunzip -c ./rat_DMF_hi_RPH_6_S54_L002_R1_001.fastq.gz) --runThreadN 8 --outSAMtype BAM SortedByCoordinate --scoreDelOpen -10000 --scoreInsOpen -10000 --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMunmapped Within --outFileNamePrefix my_alignment_`
The gunzip is required for .gz files


## Produce the count table with featureCount
`featureCounts -a annotation_file.gtf -o read_count.txt alignment.bam`

Example (relative to all this)
`featureCounts -a input.gtf -o read_count.txt alingnment_tryAligned.sortedByCoord.out.bam`





































-------- 
Make directory and copy file in
`mkdir temp`
`cp my_reads.fa temp/`
`cd temp`
Make GTF reference
`python ../makeGtf.py /home/gioele/bwaIndices/TempO-Seq_Rat_Whole_Transcriptome_1.0.fa > input.gtf`
Create hash table for star
`STAR --runMode genomeGenerate --runThreadN 7 --genomeDir ./ --genomeFastaFiles /home/gioele/bwaIndices/TempO-Seq_Rat_Whole_Transcriptome_1.0.fa --genomeSAindexNbases 4 --sjdbGTFfile input.gtf --sjdbGTFfeatureExon exon`
Run star alignment 
`STAR --genomeDir ./ --readFilesIn <(gunzip -c ./rat_DMF_hi_RPH_6_S54_L002_R1_001.fastq.gz) --runThreadN 8 --outSAMtype BAM SortedByCoordinate --scoreDelOpen -10000 --scoreInsOpen -10000 --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMunmapped Within --outFileNamePrefix my_alignment_`
Produce count table 
`featureCounts -a input.gtf -o read_count.txt my_alignment_Aligned.sortedByCoord.out.bam`
Clean everything from temp directory other than file of interest (read_count.txt)
Concatenate 7th col of tables