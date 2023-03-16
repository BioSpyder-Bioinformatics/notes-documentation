import os
import subprocess
import pandas as pd
#For sending email
import smtplib
from email.message import EmailMessage


# CHANGE THESE PARAMETERS!
# Reference index path
reference_index = '/home/gioele/bwaIndices/TempO-Seq_Rat_Whole_Transcriptome_1.0.fa'
# Outname
output_name = 'new_name'
# Are the input fastq zipped?
zipped = True
# How many threads?
threads = 8
# User email
email = 'gioelemook97@gmail.com'



# Capture environment variables:
# !!! if you move this getcwd, change the align() fn to exit the temporary directory every time as this will change otherwise and nest directories into each others 
current_directory = os.getcwd()
file_list = os.listdir()

# Get the list of files in the directory (either .gz if zipped or .fastq if not zipped!)
expected_extension = 'gz' if zipped else 'fastq'
# Get the file list to process discriminating based on extension
file_list = [file for file in file_list if file.split('.')[-1] == expected_extension]
print('Detecting files in the current working directory ', current_directory)

#file_list =["first.fastq.gz"]#, "second.fastq.gz", "third.fastq.gz"] 

print('The following files will be included in the analysis: ', ', '.join(file_list))

# Get a list of all the temporary files that get made
temp_dir_list = []

def align(filename):
    print('Starting ', filename)
    
    # Get temporary name to make a directory (filename without extension)
    temp_name = filename.split('.')[0]
    # how is the temporary directory called
    temp_dir = f'temp_{temp_name}'

    # append it to list of temporary directories
    temp_dir_list.append(temp_dir)
    # Make temporary directory 
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Check if already exists? Not sure as this is practically a check already
    #if not os.path.exists(f'{current_directory}/{temp_dir}/'):
    os.mkdir(f'{current_directory}/{temp_dir}/')
    # Copy fastq file in directory
    subprocess.Popen(f'cp {current_directory}/{filename} {current_directory}/{temp_dir}', shell=True).wait()
    # Move into temporary directory
    os.chdir(f'{current_directory}/{temp_dir}')
    
    # Make GTF reference (script is in the folder alignment_scripts in gio's home)
    subprocess.Popen(f'python /home/gioele/alignment_scripts/makeGtf.py {reference_index} > input.gtf', shell=True).wait()

    # Make Hash Table
    subprocess.Popen(f'STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {current_directory}/{temp_dir}/ --genomeFastaFiles {reference_index} --genomeSAindexNbases 4 --sjdbGTFfile input.gtf --sjdbGTFfeatureExon exon', shell=True).wait()

    # Run STAR alignment
    # This works on normal scripts but not on sh (bc of file flocking, just expanding the file)
    # _subprocess.Popen(f'STAR --genomeDir {current_directory}/{temp_dir}/ --readFilesIn {f"<(gunzip -c {current_directory}/{temp_dir}/{filename})" if zipped else filename} --runThreadN 8 --outSAMtype BAM SortedByCoordinate --scoreDelOpen -10000 --scoreInsOpen -10000 --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMunmapped Within --outFileNamePrefix {temp_name}' , shell=True)

    # Expand file with zcat if zipped
    if zipped:
        subprocess.Popen(f'zcat {current_directory}/{temp_dir}/{filename} > {current_directory}/{temp_dir}/{temp_name}.fastq', shell=True).wait()
        filename = f'{temp_name}.fastq'

    subprocess.Popen(f'STAR --genomeDir {current_directory}/{temp_dir}/ --readFilesIn {filename} --runThreadN {threads} --outSAMtype BAM SortedByCoordinate --scoreDelOpen -10000 --scoreInsOpen -10000 --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMunmapped Within --outFileNamePrefix {temp_name}' , shell=True).wait()
    

    print('Alignment complete for file ', filename)

    # Produce count table
    subprocess.Popen(f'featureCounts -a input.gtf -o read_count.txt {temp_name}Aligned.sortedByCoord.out.bam', shell=True).wait() 
    # Last argument indicates the bam file

    # Remove everything but read_count.txt
    # Get file list
    file_list = os.listdir()
    #Remove read_count.txt from file list
    file_list.remove('read_count.txt')

    # Remove all files but read_count.txt
    #subprocess.Popen(f'rm {" ".join(file_list)}').wait()
    subprocess.Popen(' '.join([f'rm {x};' for x in file_list]), shell=True).wait() # basically sends the command as rm file1; rm file2; ...
    print('Done ', filename)
    


def append_to_df(append_to, df):
    # IF append to df is empty
    if append_to.shape == (0,0):
        # append first row
        append_to['Geneid'] = df['Geneid']
    
    # Remove excess cols from df and rename it
    cols = df.columns
    # Get only cols of interest (first and last)
    cols = [cols[0], cols[-1]]
    # Get Df with cols of interest
    df = df[cols].copy()

    # Make mapping dict to rename cols (remove all that from last col name)
    rename_cols = {x:x.removesuffix('Aligned.sortedByCoord.out.bam') for x in cols}
    #remap col names
    df = df.rename(columns=rename_cols)
    # Merge dfs
    append_to = append_to.merge(df, how='outer', on='Geneid')

    return append_to


if __name__ == '__main__':
    # #trying 
    # output_name = 'new_analysis'


    # # For each file run alignment
    # for file in file_list:
    #     align(file)

    # # temp_dir_list = ['temp_first', 'temp_second', 'temp_third']


    # # Compile all read counts
    # print('Starting compiling read count table')

    # # For each temp folder, load df and append it to general one
    # # declare empty df to append to
    # append_to = pd.DataFrame()
    # # For each temporary directory, create df and append cols of interest
    # for dir in temp_dir_list:
    #     df = pd.read_csv(f'{current_directory}/{dir}/read_count.txt', delimiter='\t', skiprows=[0])
    #     append_to = append_to_df(append_to, df)

    # #Write created df to file
    # outstring = f'{current_directory}/{output_name}_count_table.csv'
    # append_to.to_csv(outstring, index=False)

    # print('Wrote count table to ', outstring)

    # # Remove all temp folders 
    # for dir in temp_dir_list:
    #     # Do it in 2 steps just in case
    #     subprocess.Popen(f'rm {current_directory}/{dir}/read_count.txt', shell=True).wait()
    #     subprocess.Popen(f'rmdir {current_directory}/{dir}/', shell=True).wait()
    #     print('Removed temporary directory ', dir)
        
    print('Processes completed, sending email')

    # Send email to the user!
    if email != None:
        # declare plain message
        msg = EmailMessage()
        # set message content
        msg.set_content("""Hello,
        You alignment on TempoSeq-Portal is completed.
        Regards,
        BioSpyder Team""")
        msg['Subject'] = 'Alignment completed'
        msg['From'] = 'random@email.com'
        msg['To'] = email

        # Send message
        s = smtplib.SMTP('localhost')
        s.send_message(msg)
        s.quit()
        
