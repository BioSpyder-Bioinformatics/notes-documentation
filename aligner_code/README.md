# Code for the aligner part of TempOSeqPortal

### Passages
For each FastQ(s) OR Fasta
- Make Temp folder
    + Copy in folder
        * Make GTF file
        * Make HashTable
    + Use Fastq + Hash Table to perform STAR alignment
        * Produces BAM file
        * Get the Count table 
    + Use BAM file + FeatureCount
        * Get count table
        * From this count table concatenate each 5th column to produce the overall count table




