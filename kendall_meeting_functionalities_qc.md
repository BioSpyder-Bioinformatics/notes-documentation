# Functionalities

*From file count table* (COUNT TABLE )
- For controls (prefix: HBrain) either triplicates or duplicates 
	+ Find the triplicates/duplicates based on initial bit of stirng (function for detect prefix!)

- Compute a R2 between replicates (new brain to new brain), show it up
- R2 between new to old controls (new brain to old brain)
	- Scatterplots of replicates 
	- Scatterplot of new vs old!

CHECK R2 value of replicates new <   x < R2 value of relicates old
                     x = average R2 of new vs old 
                     You want x to be in between, else error!!!


- Average number of mapped reads in positive controls,  










(SO there are two naming conventions for internal and customer experiments!)
- Customer experiments always have the _bspc_ and _bsnc_ inside the string!!
- Internal experiments don't always do!!




 *_bspc_*










Need to make a function to match regular expression  of bspc and match the prefixes to combine replicates 

















Other spreadsheet

Input: count table of controls



- Sheet 1:
	+ Sum column of control (basically total reads)
		* Finds average of all total reads
			- Normalisation (this total reads / total total reads)
	+ CV
		* >0 >20 >100 -> You want these numbers to be similar between controls categories
			- Looking at ~ +/-0.1 of each other
		* CV is standard deviation over average
		* If average > 0 then StD / AVG



- Sheet 2:
	+ Makes the count per million metric ( but its count per 6 million )
















# She sent unlocked spreadsheets via email!! Have a look at formulas!

