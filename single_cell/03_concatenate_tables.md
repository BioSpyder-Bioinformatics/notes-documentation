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
IF ISSUES, remove the bottom lines of sequencing_stats.txt (>14)


--------



--------

This is mainly for internal tracking, the script basically compiles all the sequencing stats into a single CSV, for it to be easier

The script is `get_summary_stats.py` and you find it in /home/salvo/A_project/SingleCell/BioRad/data_analysis/get_summary_statistics.py

The script requires the changing of the 'final_table' variable's directory, to the first directory of interest (eg for S01 to S20, you want to indicate S01)

Usage of the script (no flags, positional arguments):
- List of files (comma separated)
- Cutoff value used for the analysis
- Outfile (with .csv)

This produces a table that Dennis is interested in!


----------


- TBL2MATRIX for _UMI_!!!! To then be able to run it on seurat!!













