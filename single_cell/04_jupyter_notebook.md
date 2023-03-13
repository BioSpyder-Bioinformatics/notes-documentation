# Analysis on the Jupyter notebook
Copy of the notebook in this folder!

- Before actually starting on the notebook
	+ Download locally the file produced by the BC1 analysis (guide 02_bc1_analysis.md)
	+ Get also for the summary stats made from get_summary_stats as Dennis is interested
	+ In BC1_analysis: In the column condition, remove everything fater the underscore (the number you remove basically identifies different wells, by removing it you simply merge all the conditions to analyse together)
		* Easier way:
			- Make new col next to 'Condition' on excel
			- Copy the whole column in sublime
				+ Replace _ with \t
			- Paste back the column, the numbers after the _ are now in the new col you just made
				+ Remove that col
			- In column 'Total number of cells' (or something) make all the values _1_ !!!!!!!!!!!!!!!!!!!!!!! DONT FORGET!!!


- On Jupyter notebook
	+ Change the name of the importing file in the 2nd block!
	+ On block 6, under the graph avg detected genes vs umi, change the name of the csv output file
		* This file is one of the files Dennis is interested in

- Now make the csv file nicer...
	+ Copy the numbers from the 2nd block and paste it on the first
	+ Delete the second block
	+ Format cells to include each block and color them (mean, std, sem)
		* Select cells
			- format
				+ alignment
					* Text alignment
						- Horizontal and vertical center
					* merge cells


