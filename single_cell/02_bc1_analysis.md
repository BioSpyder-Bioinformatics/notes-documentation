# Barcode 1 Analysis

*Barcode 1* (BC1) is the first barcode ligated to the cells during the assay. Being it the first passage, the BC1 sequences can be associated to the well of reference, and consequently to a unique condition given to the well.

The BC1 analysis can be performed after doing the alignment with BIOS_single_cell, given that we need to input a list of barcodes! The list of barcodes is extracted from the count table header!

BC1 analysis is useful to check the quality of the overall experiment and how homogeneus is the distribution of different conditions (as it looks at the averages per condition)


The *TOOL* is: get_stats_bc1.py
Takes in:
- List of output folders of BIOS_single_cell run 
	+ Flag -i 
	+ /home/gioele/epa/epa_01,/home/gioele/epa/epa_02,...
- Output file (add .csv at the end!)
	+ Flag -o
- Sequences of bar code 1 and associated conditions (need to get this from people who planned he experiment)
	+ Flag -bc1



To perform BC1 analysis:
- Make the appropriate number of folders (depending on how many experiments, could be one!) (You do this bc the script outputs in the working directory)
- Copy the script into the folder (you can call it from main if you prefer, but it produces the counts into main)
- Make the file with the condition and associated BC1 sequence
	+ Get this from whoever planned the experiment (One Note or email!)
	+ VERY IMPORTANT that the conditions' name have an underscore-number at the end (eg xxx_#). This is essential to differentiate the different wells, which otherwise, as they're named the same, would over-write one another. (This ending number gets then removed for the remainder of the analysis)
	+ CHANGE NAME OF column Sample to Condition
- Ready to start!



##### How To
*Create the TSV file for -bc1*
- not sure, check with Salvo



##### If failures
If BIOS_single_cell for any reason did not produce the sequencing_stats.txt file, you can skip the passages and produce it straight with calculate_stats.py!







Move onto jupyter notebook to make averages!















