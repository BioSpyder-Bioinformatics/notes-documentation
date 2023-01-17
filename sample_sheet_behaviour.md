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
Breakdown:
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
	+ User-assigned index of reference (A-D800/A-T900)
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
































# To talk about
- Do we want to allow for the selection of more than one plate for the 24 and 48 wells plates too? (If yes we need the workaround for the indexes)

- For the 24 well file structure - when the user has to insert the sample names, there is no enforcement for leaving the space for the pos/neg controls. Would this be the last two rows?

- Does the 48 well layout run on a 96 well plate? In the sample sheet there are 48 wells for the samples and 16 extra for the controls

- For the 96 well plate, do I need to enforce leaving the last two rows empty?

- Is there a difference between the 800 and 900 probes? If yes should I prompt the user to pick one or the other?
