# Transfer data with Teamviewer

- Access teamviewer with the given password
- Open the collapsed icon for the sequencer and have a look at the top right, if there is a green [V]
- If there is the data is ready, collapse the application and open filezilla

In filezilla:
- Get the latest run folder
	+ Starts with 2303##_NN_....
	+ It is the last file of the series named like that
	+ Inside that
		* Folder alignment1
			- Another folder (just 1)
				+ Inside there there is the folder fasta (the one of interest)
- Connect to kevin
	+ File, create connection, kevin, ok 
	+ Make new directory where to store the data
- Select everything in the alignment folder, drag and drop in kevin!





# Preliminary analysis
Run fastqc

Move to the folder of interest
`fastqc *.gz`

Have a look at the HTML report!



