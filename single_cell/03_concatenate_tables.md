# Concatenate Tables

This step is not dependent on 02_bc1_analysis, but it is dependent on 01_bios_single_cell_analysis. 

It basically concatenates all the read_count and umi_count which were produced by BioS_Single_Cell into one CSV file (which is the complete gene count table)

The script is `concatenate_index_alignments.py`
It takes in 
- f
	+ Comma separated list of folders where the outputs of BioS_Single_Cell are
- o
	+ Output file name, no need to add the csv
- s 
	+ Species of interest (likely human)


You get the concatenated reads and umi tables

- TBL2MATRIX for UMI!!!! To then be able to run it on seurat!!













