# BIOS_single_cell analysis
This is the entry point of the single cell analyisis!

In order to perform this analysis, there a prerequisite passage, which is to concatenate the sequencing files into one per read per sample.
`zcat filename | gzip` (THIS IS ONLY MADE IF THERE IS 1+ RUNS TO CONCATENATE)



##### Installation
If not installed yet, you can either copy the file from Salvo's home (/home/salvo/SingleCell/BioS_single_cell), or clone it from GitHub (better tbh, so you can pull when needed instead of copying again).
Install the script use the helper script `bash install.sh`

Check that the installation was successful with 
`./BioS_single_cell -h`


### Usage
The software takes in several inputs:
- -1
	+ Read 1 of the FastQ files
- -2
	+ Read 1 of the FastQ files
- -ht
	+ Reference for the probe pannel (in this case the BWA Index)
- -a
	+ Annotation file of the probe pannel (GTF extension)
- -o
	+ Output folder
- -t
	+ Number of processors to use
- -n
	+ Minimum number of reads a cell should have (cutoff)
- -species
	+ Comma separated list of species

Extra resources:
- Hash tables for BWA
	+ /home/salvo/bwaIndices/
- GTF annotation files
	+ /home/salvo/annotation_files/

*Command Example*
`./BIOS_single_cell -1 read1.fastq.gz -2 read2.fastq.gz -o output_folder -n 10000 -ht /home/salvo/bwaIndices/TempO-Seq_Human_Whole_Transcriptome_2.1.fa -a /home/salvo/annotation_files/TempO-Seq_Human_Whole_Transcriptome_2.1.gtf -species human -t 50`


### Output
- `logfile.txt`
- `sequencing_stats.txt`
- `{species}_variable_threshold_plot.jpg`
- `{species}_mapped_umi_vs_detected_genes.jpg`
- `{species}_umi_count.csv`
- `{species}_read_count.csv`
- `barcode_quantification.csv`









