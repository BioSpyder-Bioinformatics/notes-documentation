# Check probes efficiency
Step-by-step guide to prepare the probe data for the RAS-RAF project. 
Required
- Count table
- Mapped-unmapped
- Experiment plan

*WT -> wild type*

## Probes naming convention
```{gene_name}_s_c{base position}{WT base}-{mutant base}_{often another base}_{sometimes suffix}

KRAS_s_c34G-A_c35G-A_2D
```

## Samples naming convention

*MAF -> minor allele frequency (referring to mutant) It refers to the % of mutant probe spiking the sample*
*gDNA -> genomic DNA, in these samples there should be no probe signal*
*NC* -> negative control


# Analysis On the count table
Only a few probes are analysed for efficiency at the time, however all of them are sequenced. (Get which probes you need to focus on from the Experiment plan docx)

We want to produce a table of signal WT to mutant ratio for each of the conditions and SNPs. We need
- Ratio wt-mutant
- Ratio between replicates
- Find the best signal 
	+ Present all this in a powerpoint
- For the mapped-unmapped file 
	+ Need the % mapped barplot

Steps:
- In the first step we want to get the ratio between the WT and the mutant base for each of the conditions
	+ The ratio is obtained as mutant/WT•100 
		* For example, if in the sample sheet we find that one of the probes of interest is KRAS Exon 2 - 34G-A
			- Find the mutant probe for it based on the name, in this case KRAS_s_c34G-A_c35G_2D
			- Find the WT probe, very similar to the mutant one, only changes that does not give the mutation base, in this case KRAS_s_c34G_c35G_2D
	+ Once you have the ratio, transpose it in another sheet and plot an histogram of all the conditions' ratio
	+ If the conditions are repeated, do an average of the ratios
		* Prepare a table with the average of the ratios (if any), with a row for each SNP and a column for each condition, grouped by NC - gDNA - MAF%

# Analysis on count table
Get the % of mapped reads per sample in a barplot

group by condition boxplot showing difference between conditions



