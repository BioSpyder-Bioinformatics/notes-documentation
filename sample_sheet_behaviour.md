# On sample sheet download

*Options*
- Plate layout
	+ 24-48-96 wells
- Number of plates (IF 96 well plate, otherwise just 1)
	+ 1 to 8

- Project name
- Experiment name
- Additional comments
(These three pieces of metadata don't really get reported on the sample sheet)



You then download in CSV format 



## CSV file structure
_When downloading 24 well layout_
You get the classic 96 well plate layout, but with only 4 columns and 6 rows, the other rows are marked as NoIndex (last 2 rows are NoIndex1 and NoIndex2)
!!! Issue: When downloading the 24 well template, there is no space for the pos neg controls, are the last two rows (eg NoIndex1 and NoIndex2) responsible for holding those? If yes, do we want to change the naming convention



_When downloading 48 well layout_
You get the classic 96 well plate layout, but with 8 columns (last 4 are marked as NoIndex) and the classic 8 rows. (Be careful, there are 2 extra rows, I suppose for positive/negative controls, but still, is this run on a 96 well plate?? Otherwise there's no physical space to fit)



_When downloading96 well layout_
You get the classic 96 well plate layout (12 cols, 8 rows) (Be careful, some of these are meant to be controls, you don't want allow the user to use it all!)
If you have selected more than 1 plate, these get stuck one on top of the other with only 1 row of space in between






# On reupload

*For 24 well plate*

It pretty much displays the whole 96-well plate but with whe empty columns greyed out and showing NA. For index selection, this is pre-assigned as only 1 exists (Check the CSV though)
On confirm, it prompts to the page for download in the different instruments (MiSeq-NextSeq-MiniSeq-HiSeq)
It behaves as the 96 well plate, it pretty much melts the table and appends the indexes
!!! Very concerned about pos/neg control position as it doesn't get reported in the sample sheet
!!! Need to ask where the controls fit in here and if I have to enforce some sort of restriction



*For 48 well plate*

Same as the 24-well. It displays the 96-well plate but with the empty columns greyed out. Index selection is preassigned as well as only 1 exists for 48 well plates layouts. Very annoying -> if you fill the last 2 rows (meant for controls) it will not give any error and just output empty CSVs when getting the sample sheets for NGS machinery. Table is melted normally otherwise. (Even more concern for the 24 well plate's spacing for controls :-( )



*For 96 well plate*

Allows you to pick the plate's indexes. The only control there is is not allowing the user to select the same index twice during plate's selection (per experiment). No enforcement of leaving 2 columns empty for controls, this might need to be added.
On confirm you get redirected to the page for the NGS machinery's Sample Sheet download. The downloaded CSV has metadata specific to the machinery and the tables melted one at the time (eg melt table1, melt table 2 etc)








# NGS-specific sample sheet formatting

*For all filenames*
`(selected_project_name)_(instrument_name).csv`

### MiSeq Sample Sheet
`[project_name]_SR_MS.csv`

CSV template, described bellow 
```csv
"[Header]",""," "," "," "," "," "," "," "," "
"IEMFileVersion","4"," "," "," "," "," "," "," "," "
"Investigator","BioSpyder"," "," "," "," "," "," "," "," "
"Project Name","second_project"," "," "," "," "," "," "," "," "
"Experiment Name","second_experiment"," "," "," "," "," "," "," "," "
"Date","01/16/2023"," "," "," "," "," "," "," "," "
"Workflow","GenerateFASTQ"," "," "," "," "," "," "," "," "
"Application","FASTQ Only"," "," "," "," "," "," "," "," "
"Assay","Nextera"," "," "," "," "," "," "," "," "
"Description","MiSeq"," "," "," "," "," "," "," "," "
"Chemistry","Amplicon"," "," "," "," "," "," "," "," "
"Additional Comments","comment_sample","","","","","","","",""
"[Manifests]",""," "," "," "," "," "," "," "," "
"",""," "," "," "," "," "," "," "," "
"[Reads]",""," "," "," "," "," "," "," "," "
"50",""," "," "," "," "," "," "," "," "
"[Settings]",""," "," "," "," "," "," "," "," "
"CustomIndexPrimerMix","C2"," "," "," "," "," "," "," "," "
"",""," "," "," "," "," "," "," "," "
"[Data]",""," "," "," "," "," "," "," "," "
"Sample_ID","Sample_Name","Sample_Plate","Sample_Well","I7_Index_ID","index","I5_Index_ID","index2","Sample_Project","Description"
"1","A1","A","A01","R801","AAGACTCTT","F801","AAGGTGTTT","",""
```

*Breakdown*
(Section) Header
- IEMFileVersion -> 4 
	+ (?)
- Investigator -> BioSpyder 
	+ not sure if it changes processes
- Project Name -> 
	+ Based on user selection
- Experiment Name ->
	+ User-based
- Date -> 
	+ Date in format MM/DD/YYYY
- Workflow -> GenerateFASTQ
- Application -> FASTQ Only
- Assay -> Nextera
- Description -> MiSeq
- Chemistry -> Amplicon
	+ All these above are probably required by the NGS machinery to operate specific protocols
- Additional comments -> 
	+ User's comments

(Section) Manifests
- Empty (not sure why)

(Section) Reads
- 50
	+ For NGSM

(Section) Settings
- CustomIndexPrimerMix -> C2
	+ Seems to be a constant
	
(Section) Data
All the melted tables, with columns as:
- Sample_ID
	+ from 1 to N
- Sample_Name
	+ User-given sample names
- Sample_Plate
	+ User-assigned index of reference (A-D800/A-T900)(Only letter asssigned)
- Sample_Well
	+ Well assigned to sample, in format RCC (row-col-col, eg. A01) 
- I7_Index_ID
	+ ID of I7 index in format R### (# is a number)
- index
	+ I7 index sequence
- I5_Index_ID
	+ ID of I5 index in format F###
- index2 
	+ I5 index sequence
- Sample_Project
	+ Empty by default
- Description
	+ Empty by default






### HiSeq Sample Sheet
Surprisingly this file does not have any metadata, it's exactly the same as the Data section of MiSeq. This means that user's metadata such as experiment name and description is lost. The CSV file is composed by the melted table with the columns as: (copy pasted from above)
(The file starts with the column names so no [stuff in bracket for the machine])

*CSV File Trace*
```csv
"Sample_ID","Sample_Name","Sample_Plate","Sample_Well","I7_Index_ID","index","I5_Index_ID","index2","Sample_Project","Description"
"1","A1","A","A01","R801","AAGACTCTT","F801","AAGGTGTTT","",""
"2","B1","A","B01","R801","AAGACTCTT","F803","AACTACAGC","",""
```

*Breakdown*
(Section) Data
All the melted tables, with columns as:
- Sample_ID
	+ from 1 to N
- Sample_Name
	+ User-given sample names
- Sample_Plate
	+ User-assigned index of reference (A-D800/A-T900)(Only letter asssigned)
- Sample_Well
	+ Well assigned to sample, in format RCC (row-col-col, eg. A01) 
- I7_Index_ID
	+ ID of I7 index in format R### (# is a number)
- index
	+ I7 index sequence
- I5_Index_ID
	+ ID of I5 index in format F###
- index2 
	+ I5 index sequence
- Sample_Project
	+ Empty by default
- Description
	+ Empty by default







### NextSeq Sample Sheet
This is almost identical to the MiSeq sample sheet, only a couple of parameters for the machine are changed (application and description in the header's data). 

*CSV File Trace*
```csv
"[Header]",""," "," "," "," "," "," "," "," "
"IEMFileVersion","4"," "," "," "," "," "," "," "," "
"Investigator","BioSpyder"," "," "," "," "," "," "," "," "
"Project Name","second_project"," "," "," "," "," "," "," "," "
"Experiment Name","second_experiment"," "," "," "," "," "," "," "," "
"Date","01/16/2023"," "," "," "," "," "," "," "," "
"Workflow","GenerateFASTQ"," "," "," "," "," "," "," "," "
"Application","NextSeq FASTQ Only"," "," "," "," "," "," "," "," "
"Assay","Nextera"," "," "," "," "," "," "," "," "
"Description","NextSeq"," "," "," "," "," "," "," "," "
"Chemistry","Amplicon"," "," "," "," "," "," "," "," "
"Additional Comments","comment_sample","","","","","","","",""
"[Manifests]",""," "," "," "," "," "," "," "," "
"",""," "," "," "," "," "," "," "," "
"[Reads]",""," "," "," "," "," "," "," "," "
"50",""," "," "," "," "," "," "," "," "
"[Settings]",""," "," "," "," "," "," "," "," "
"CustomIndexPrimerMix","C2"," "," "," "," "," "," "," "," "
"",""," "," "," "," "," "," "," "," "
"[Data]",""," "," "," "," "," "," "," "," "
"Sample_ID","Sample_Name","Sample_Plate","Sample_Well","I7_Index_ID","index","I5_Index_ID","index2","Sample_Project","Description"
"1","A1","A","A01","R801","AAGACTCTT","F801","AAACACCTT","",""
```


*Breakdown*
(Section) Header
- IEMFileVersion -> 4 
- Investigator -> BioSpyder 
- Project Name -> 
- Experiment Name ->
- Date -> MM/DD/YYYY
- Workflow -> GenerateFASTQ
- Application -> NextSeq FASTQ Only *THIS CHANGES*
- Assay -> Nextera
- Description -> NextSeq *THIS CHANGES*
- Chemistry -> Amplicon
- Additional comments -> 

(Section) Manifests
- Empty 

(Section) Reads
- 50

(Section) Settings
- CustomIndexPrimerMix -> C2

(Section) Data
All the melted tables, with columns as:
- Sample_ID
- Sample_Name
- Sample_Plate
- Sample_Well
- I7_Index_ID
- index
- I5_Index_ID
- index2 
- Sample_Project
- Description





### MiniSeq Sample Sheet
This is the most peculiar one, it returns a zipped file with a TSV *and* a CSV
!!!!!!! *__BE VERY CAREFUL__*, this one uses reverse complement sequences for the I5 INDEXES

###### TSV
TSV file seems to be the same regardless of the plate size. It looks more like a manifest file, with metadata on top, and indices name and sequence after
*Very strange*
It gives 96 'R' indexes and only 64 'F' indexes _CHECK INDEXES WITH NEW ONES_

*TSV Template*
```tsv
[Kit]			
Name	BioSpyder900		
Description	"BioSpyder, Inc. 900 series"		
IndexStrategy	All		
ReadType	All		
DefaultReadLength1	50		
DefaultReadLength2	0		
			
[Resources]			
Name	Type	Format	Value
			
[Indices]			
Name	Sequence	IndexReadNumber	
R901	GTATTATTG	1	
…
R996	GATAGCGGA	1	
			
F901	AACTACAGC	2	
…
F964	CGCTACAAC	2	
			
[SupportedModules]			
GenerateFastQWorkflow	
```


*TSV file breakdown*
(Section) Kit
- Name -> BioSpyder900
- Description -> "BioSpyder, Inc. 900 Series"
- IndexStrategy -> All
- ReadType -> All
- DefaultReadLength1 -> 50
- DefaultReadLength2 -> 0

(Section) Resources
(all on the same line as if they were columns, but nothing under it)
Name	Type	Format	Value

(Section) Indices
- Name -> index name (R/F###)
- Sequence -> Actual index sequence
- IndexReadNumber -> R sequences have 1 assigned, F sequences have 2 assigned

(Section) SupportedModules
- GenerateFastQWorkFlow
	+ This is just there on its own






###### CSV
*CSV Template*
```csv
[Header],,,,,,
Experiment Name,second_experiment,,,,,
Date,01/16/2023,,,,,
Module,GenerateFASTQ - 2.0.1,,,,,
Workflow,GenerateFASTQ,,,,,
Library Prep Kit,BioSpyder900,,,,,
Chemistry,Amplicon,,,,,
[Reads],,,,,,
50,,,,,,
[Settings],,,,,,
adapter,CTGTCTCTTATACACATCT,,,,,
[Data],,,,,,
Sample_ID,Description,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project
A1,,R801,AAGACTCTT,F801,AAGGTGTTT,
```




*CSV Breakdown* (No spacing at all throughout document)
(Section) Header
- Experiment Name -> 
	+ user-defined
- Date -> 
	+ automatic in MM/DD/YYYY format
- Module -> GenerateFASTQ - 2.0.1
- Workflow -> GenerateFASTQ
- Library Prep Kit -> BioSpyder900
- Chemistry -> Amplicon

(Section) Reads
- 50
	+ It's just there by itself

(Section) Settings
- adapter -> CTGTCTCTTATACACATCT

(Section) Data
- Sample_ID
	+ The ones defined by the users
- Description
	+ Empty by default
- I7_Index_ID
	+ 'R' index ID
- index
	+ actual index sequence
- I5_Index_ID
	+ THIS IS MEANT TO BE REVERSED, CHECK BETTER WHAT SEQUENCES YOURE GETTING 
- index2
	+ actual sequence 
- Sample_Project
	+ Empty by default











# To talk about
- Do we want to allow for the selection of more than one plate for the 24 and 48 wells plates too? (If yes we need the workaround for the indexes)
	+ No

- For the 24 well file structure - when the user has to insert the sample names, there is no enforcement for leaving the space for the pos/neg controls. Would this be the last two rows?
	+ ASK GARRET IF WE NEED TO LEAVE WELLS EMPTY

- Does the 48 well layout run on a 96 well plate? In the sample sheet there are 48 wells for the samples and 16 extra for the controls
	+ Ask Garret

- For the 96 well plate, do I need to enforce leaving the last two column empty?
	+ Ask Garret

- Is there a difference between the 800 and 900 probes? If yes should I prompt the user to pick one or the other?
	+ STOP USING 800 - CHECK WHICH FLAVOURS 

- HiSeq sample sheet loses user's metadata
	+ Is being discontinued, not too much of an issue

- MiniSeq TSV has uneven number of R and F indexes :( 



- How many pos/neg controls per sample sheet. Actively enforce this in pre-determined wells 






- If you order a 48 plate you can fill it all, how many for the 96? for pos neg 
- Do users only use flavour 800 or 900? 
- Why in flavour 900 there are E-F-G-H missing

















# Meeting with Garrett
- Customers can ask for up to 8 96-wells plates
	+ Lower priority is 7 96-well plates plus 1 48-well plate
- For 48 wells plate we give only 6 columns, the controls are managed in the rest of the table
- Only offer the 900 series flavors (e to l)
- In the 96-well plates
	+ 90 samples
	+ 6 controls
		* 4 positive
		* 2 negative
- For special orders, eg. 384-well (4 • 96-well plates)
	+ 3 plates have 6 controls each (18 total)
	+ 1 plate has no control

- Give user possibility to upload additional file with sample name and additional comments

*For future*
- Give possibility to create heatmap with mapping rate and gene count





