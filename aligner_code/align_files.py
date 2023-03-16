import sys
import os
import subprocess

# Entry folder
entry_folder = '.'
# Reference index path
reference_index = '/home/gioele/bwaIndices/TempO-Seq_Rat_Whole_Transcriptome_1.0.fa'
# Outname
output_name = 'new_name'
# Are the input fastq zipped?
zipped = True




# Some variables:
current_directory = os.getcwd()
file_list = os.listdir()

# Get the list of files in the directory (either .gz if zipped or .fastq if not zipped!)
expected_extension = 'gz' if zipped else 'fastq'
file_list = [file for file in file_list if file.split('.')[-1] == expected_extension]
print('This is the file list i found!')


# Get a list of files in directory (artificial at this point)
file_list =["first.fastq.gz", "second.fastq.gz", "third.fastq.gz"] 

print('The following files will be included in the analysis: ', ', '.join(file_list))




def align_star(filename):
    print('starting ', filename)
    # Activate conda environment
    #subprocess.Popen('conda activate aligner', shell=True)

    temp_name = filename.split('.')[0]
    # how is the temporary directory called
    temp_dir = f'temp_{temp_name}'
    # Make temporary directory 
    os.mkdir(f'{current_directory}/{temp_dir}/')
    # Copy file in directory
    subprocess.Popen(f'cp {current_directory}/{filename} {current_directory}/{temp_dir}', shell=True)
    # Move into temporary directory
    os.chdir(f'{current_directory}/{temp_dir}')
    
    # Make GTF reference (script is in the folder alignment_scripts in gio's home)
    subprocess.Popen(f'python /home/gioele/alignment_scripts/makeGtf.py {reference_index} > input.gtf', shell=True)

    # Make Hash Table
    subprocess.Popen(f'STAR --runMode genomeGenerate --runThreadN 7 --genomeDir ./ --genomeFastaFiles {reference_index} --genomeSAindexNbases 4 --sjdbGTFfile input.gtf --sjdbGTFfeatureExon exon')

    # Run STAR alignment
    subprocess.Popen(f'STAR --genomeDir ./ --readFilesIn {f"<(gunzip -c ./{filename})" if zipped else filename} --runThreadN 8 --outSAMtype BAM SortedByCoordinate --scoreDelOpen -10000 --scoreInsOpen -10000 --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMunmapped Within --outFileNamePrefix {temp_name}')

    # Produce count table
    subprocess.Popen(f'featureCounts -a input.gtf -o read_count.txt {temp_name}Aligned.sortedByCoord.out.bam')

    print('done ', filename)
    



# Depending on the aligner chosen

# For each file
for file in file_list:
    align_star(file)
    




# REMOVE FILE
#os.unlink(path)
# Remove directory
#os.rmdir(path)




