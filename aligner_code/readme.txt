# temposeq_aligner
TempoSeq aligner is a tool to align your TempoSeq sequences and extract a feature count table. It allows for alignment with star, bwa and kallisto. It is dependent on the script makeGtf.py and on the packages star, kallisto, bwa and samtools.

### Usage
python temposeq_aligner.py [options]
*Options*:
- -i, --input-directory
    + Required input. Directory where the fastq files are. Requires either '.' for current directory or the full directory path.
- -a, --aligner
    + Required input. Aligner of preference. Options: bwa, star, kallisto.
- -g, --reference-genome
    + Required input. Reference genome of preference. Input either one of the options, or the complete path to the .fa file. Options: rat_w_1.0, human_s1500_1.2, human_w_2.0, human_w_2.1, mouse_s1500_1.2, mouse_w_1.0-
- -o, --output-name
    + Required input. Prefix for the output name for the count table.
- t, --threads
    + Not required. Number of threads used by the aligners. Default: 8.
- -u, --unzipped
    + Not required. Include if the input files are unzipped - Expected file extension _.fastq_. Default: False - Default behaviour: Expected file extension _.fastq.gz_.

*Quirks*
- Can only run on Alpha as of now as it is dependent on makeGtf, ask salvo what's best doing (maybe import it on)
- The shortcuts for the reference genomes can also only be used on alpha as they point to gio's home path.
- Find stuff to manually adjust to run in different servers with ctrl+f #!!!
- Need to run in environment with dependencies installed
