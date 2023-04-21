# Data analysis in Biospyder (By Kendall)

Steps
- Review traces (for pcr I suppose)
	+ They make sure they overlap if using the same material

- Mapped-Unmapped file
	+ They do the average of the mapped percentage of reads (average because they run triplicates for controls)
		* You want ~>96% in average and not much more of 1% difference between replicates
	+ Usually 2 pos controls but they do 3 to hit as many genes as possible + 1 negative control
	+ Make correlation graphs 
		* Between different control replicates, you want linear graphs (scatter plots)
	+ Check variation between old and new brain positive controls is about as much as the one between the triplicates controls.
		* Check that the proportion of read counts is also close
		* When performing this analysis we do normalise the data based on 6millions for some reason
	
- In temposeqr check the data in graphs
	+ Barplot for read counts (check its not skewed)
	+ Dendrogram to see homology association (we want to see relicates grouped together)
	+ PCA (same reason as before)
	+ Log2 Scatter plot (with R2 for association)

Find more in Biospyder shared docs, Documents -> Manufacturing -> QC Manufacturing -> 2023 -> Find the protol of reference

























