# Generic pipeline
- TempO-SeqR _makes_ Sample Sheet 
- Sequencer _makes_ Fastqs 
- TempO-SeqR _produces_ Alignment -> _which produces_ Counts Table 
- Counts table is _used for_ QC 
	+ Descriptive statistics 
	+ Barplot 
	+ Dendogram
	+ Scatterplot 
	+ PCA
- Counts table is used for *Differential expression analysis*




# Sample Sheet
The sample sheeet is a 'library' that associates _specific index sequences_ with _appropriate samples_.
The sample sheet generator is used to convert a _plate map of sample names_ into a _Sample Sheet for NGS instruments_

The NGS leverages the Sample Sheet in order to produce separate FastQ files, one for each sample (or group of index sequences)

Once obtained the FastQ files for each sample, these can be aligned to a reference genome to produce a count table

There are the options to
- Check if a manually produced sample sheet is compatible with the assay (enforces checks for the naming and indices)
- Generate a sample sheets
	+ Requires
		* Plate layout (24,48,96 wells), number of plates, project/experiment name. After downloading the CSV, fill up the cells with the name of the samples and reupload. On reupload it asks to select index, ask about this (i suppose it only determines the index sequence based on the index gave as it has to be unique)
	+ Produces
		* Sample sheets for MiSeq, NextSeq, MiniSeq, HiSeq
			- All sample sheets but the one of HiSeq have metadata including
				+ Workflow, assay, application, type of output and adapter (each has its specifics)
			- For all sample sheets:
				+ It takes the CSV and 'melts' the table on the column side (so that the cells are represented as A1, B1, C1, ..., A4, B4, C4, ...)
				+ Each column receives the same index 1 (eg A1, B1, C1)
				+ Each row receives the same index 2 (eg A1, A2, A3), this allows to have a unique combination across all the plate

!Have a look for the Illumina tool Salvo mentioned to double check sample sheets!!

!Edit straight on web page, with automated checks on input (but for later on)

Ability to upload different files formats eg
- Generic
	+ CSV
- Propietary 
	+ xlsx (excel)
- Non-propietary
	+ ods (excel-like)
	+ This is mainly if we decide to append a second sheet with metadata to the original sample sheet. In this case we need a workaround for CSV uploads (will be later on)
	
Given that the number of indexes will grow in the future, set this up to have 8 96-wells plate at the time instead of only 4



# Count table
Table showing the counts (signal levels) for each probe, *in* each sample. This is what is used for the analysis
This is the result of alignment of the produced FastQ to a reference table. It displays the amount of reads mapped to a particular probe for each sample.

When displaying the count table:
- Add filters (radio/checkbox buttons) to isolate expression of specific groups/genes
	+ Allow for inserting comma separated list of genes/samples of interest
	+ Filter table by
		* Average read count
		* Normalised read count eg count per million
	+ Ability to download the table with this specific filtering



# Components of interest for app

*Side menu (all pages)*
- Home
- Sample sheet generator
- Aligner (?)
- Quality control
- Differential expression
- Support
- Logout

*Landing Page*
- Some sort of animation/image while page loads
- Alignment history
- Latest updates (if any)
- User's info

*Sample Sheet Generator*
_Check manual to see what is required_
- This has to have a lot of forceful checks on form
- Maybe do that can be filled both on page or downloaded and reuploaded!



*Aligner*
Ignore for now, probably do an GUI to do locally (Maybe refer user to download link and guide for tool)



*Upload section for count table*
Required as users can no longer align online 
This is a double upload window asking for the count table and Sample sheet. Sample sheet is optional but recemmended for identifying and grouping samples as they were placed in the plate. (The sample sheet is the plate layout one, the one you re-upload to obtain the one to feed to the sequencer)



*QC Tools*
For all graphs you need to allow to select a subset of data to be plotted (Also doublecheck how to do in graph selection). 
For all graphs and tables you need to allow for download and selection download
_For graphs/statistics_ select samples with checkboxes or even better side by side list that allow you to swop elements between them


_Display Count Table_
(So that user can have a look at the table straight on the app. This requires a search function or something)

_Sample statistics_
Displays the total mapped reads, average reads per probe and sample type. Here user can specify positive and negative controls
Other parameters useful for the analysis include:
- P80
	+ Metric describing number of probes required to capture 80% of change in expression. This is useful as QC tool in case number of probes required vary too much between samples (could be given by sequencing depth - required adjustments in attenuation etc)
- Number of detected genes
- Number of probes with at least count 1
	+ This is mainly for single cell
- Statistics by group
	+ eg pos-neg-control-group x/y/z
	+ Useful for mediating depth difference between samples


_Barplot_
Shows number of reads per sample

_Scatterplot_
Shows degree of correlation between two samples and displays R-squared value. Selecting in graph shows the names of the selected genes as well as the level of expression in the two plotted samples. Allows to show the plot in linear, log2 or log10 scale

_Dendrogram_
Show clustering of samples based on similarity (low reads samples not displayed <1000)

_PCA_
Shows clustering of samples based on similarity. This plot is also offered in 3D (when 3 PC are selected)
Samples with same labels need to be displayed with the same color. (Labels are given in the column 'description' of the sample sheet. Maybe add option to label dynamically online)
Need to understand how to extract PC metadata to determine components


_Boxplot_
Boxplots repressenting the number of genes. Useful to be able to colour boxes based on group, plus add statistics on hover (maybe checkboxes for user to select which metrics they're interested to see on hover)


_Dendrogram_
Have a look for some fancy template. Having a circular dendrogram would be better








*Differential expression tools*

_MAPlot_
Visual representation of differentially expressed data. Plots baseMean (mean of normalised counts across all samples) vs log2FoldChange. Requires selection of one control sample and one experimental one. Samples with p<0.05 have different colour dots.
Can filter reults by p-value and/or BaseMean value.
Once again select dots to get name of the gene, baseMean, log2FoldChange and p-val



_DeSeq_
Page representing the raw data in a CSV table obtained by the DESeq analysis








# Naming restrictions

*FastQ*
- No spaces or special characters in the filename
- Spaces are replaced by underscores or hyphens, reach an agreement and stick with it (suggest underscores)
- File extension can only be .fastq, .fastq.gz, .fq


*Sample naming*
- No spaces or special characters in sample name
- Again, spaces replaced by _ or -
- Unique sample names
- No use of BSPC BSNC (biospyder pos/neg control)
	+ Any other??

























