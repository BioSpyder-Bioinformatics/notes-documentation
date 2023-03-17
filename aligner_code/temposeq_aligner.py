import os
import subprocess
import pandas as pd
#For sending email
import smtplib
import argparse


# Get a list of all the temporary files that get made

def align_star(filename, reference_index, threads, zipped, temp_dir_list, current_directory):
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
    return temp_dir_list
    


def align_bwa(filename, reference_index, threads, zipped, temp_dir_list, current_directory):
    print('Starting', filename)

    # Get temporary name to make a directory (filename without extension)
    temp_name = filename.split('.')[0]
    # how is the temporary directory called
    temp_dir = f'temp_{temp_name}'

    # append it to list of temporary directories
    temp_dir_list.append(temp_dir)
    # Make temporary directory 
    os.mkdir(f'{current_directory}/{temp_dir}/')
    # Copy fastq file in directory
    subprocess.Popen(f'cp {current_directory}/{filename} {current_directory}/{temp_dir}', shell=True).wait()
    # Move into temporary directory
    os.chdir(f'{current_directory}/{temp_dir}')
    
    # Make GTF reference (script is in the folder alignment_scripts in gio's home)
    subprocess.Popen(f'python /home/gioele/alignment_scripts/makeGtf.py {reference_index} > {current_directory}/{temp_dir}/input.gtf', shell=True).wait()

    #make bwa index
    subprocess.Popen(f'bwa index {reference_index}', shell=True).wait()

    # Expand file with zcat if zipped
    if zipped:
        subprocess.Popen(f'zcat {current_directory}/{temp_dir}/{filename} > {current_directory}/{temp_dir}/{temp_name}.fastq', shell=True).wait()
        filename = f'{temp_name}.fastq'

    # Align with BWA
    subprocess.Popen(f'bwa mem -v 1 -c 2 -L 100 -t {threads} {reference_index} {current_directory}/{temp_dir}/{filename} > {current_directory}/{temp_dir}/alignment.sam', shell=True).wait()

    # Convert sam to bam
    subprocess.Popen(f'samtools view -bS {current_directory}/{temp_dir}/alignment.sam > {current_directory}/{temp_dir}/alignment.bam', shell=True).wait()

    print('Alignment complete for file ', filename)

    # Run feature counts
    subprocess.Popen(f'featureCounts -a {current_directory}/{temp_dir}/input.gtf -o {current_directory}/{temp_dir}/read_count.txt {current_directory}/{temp_dir}/alignment.bam', shell=True).wait()

    # Remove everything but read_count.txt
    # Get file list
    file_list = os.listdir()
    #Remove read_count.txt from file list
    file_list.remove('read_count.txt')

    # Remove all files but read_count.txt
    #subprocess.Popen(f'rm {" ".join(file_list)}').wait()
    subprocess.Popen(' '.join([f'rm {x};' for x in file_list]), shell=True).wait() # basically sends the command as rm file1; rm file2; ...
    print('Done ', filename)
    return temp_dir_list



def align_kallisto(filename, reference_index, threads, zipped, temp_dir_list, current_directory):
    print('Starting', filename)

    # Get temporary name to make a directory (filename without extension)
    temp_name = filename.split('.')[0]
    # how is the temporary directory called
    temp_dir = f'temp_{temp_name}'

    # append it to list of temporary directories
    temp_dir_list.append(temp_dir)
    # Make temporary directory 
    os.mkdir(f'{current_directory}/{temp_dir}/')
    # Copy fastq file in directory
    subprocess.Popen(f'cp {current_directory}/{filename} {current_directory}/{temp_dir}', shell=True).wait()
    # Move into temporary directory
    os.chdir(f'{current_directory}/{temp_dir}')

    #Make kallisto index
    subprocess.Popen(f'kallisto index -i kallisto.index {reference_index}', shell=True).wait()

    # Expand file with zcat if zipped
    if zipped:
        subprocess.Popen(f'zcat {current_directory}/{temp_dir}/{filename} > {current_directory}/{temp_dir}/{temp_name}.fastq', shell=True).wait()
        filename = f'{temp_name}.fastq'

    # Run kallisto quantifier
    subprocess.Popen(f'kallisto quant --single -i kallisto.index -o {current_directory}/{temp_dir} -l 50 -s 1 -t {threads} --bias --single-overhang {current_directory}/{temp_dir}/{filename}', shell=True).wait()

    # Remove everything but read_count.txt
    # Get file list
    file_list = os.listdir()
    #Remove abundance.tsv from file list
    file_list.remove('abundance.tsv')

    # Remove all files but read_count.txt
    subprocess.Popen(' '.join([f'rm {x};' for x in file_list]), shell=True).wait() # basically sends the command as rm file1; rm file2; ...
    print('Done ', filename)
    return temp_dir_list








def append_to_df(append_to, df, filename = None):
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
    # If filename is not specified capture it from the DF
    if filename == None:
        # Make mapping dict to rename cols (remove all that from last col name)
        rename_cols = {x:x.removesuffix('Aligned.sortedByCoord.out.bam') for x in cols}
        #remap col names
        df = df.rename(columns=rename_cols)
    else:
        # This is specific to kallisto
        df.est_counts = df.est_counts.astype(int)
        df = df.rename(columns={'est_counts': filename})
        
    # Merge dfs
    append_to = append_to.merge(df, how='outer', on='Geneid')

    return append_to



# def send_email(recipient):
#     body = 'Subject: Alignment complete \n\n' + 'Hello,\n\nYour alignment has completed!' + '\n\nRegards,\nBioSpyder Team'
#     sender = 'biospyder_aligner@outlook.com'

#     # Initialise smtp objects
#     try:
#         smtpObj = smtplib.SMTP('smtp-mail.outlook.com', 587)  #smtp.office365.com
#     except Exception as e:
#         print(e)
#         smtpObj = smtplib.SMTP_SSL('smtp-mail.outlook.com', 465)

#     smtpObj.ehlo()
#     smtpObj.starttls()
#     smtpObj.login('tempo-portal@biospyder.com', '8nXiM7eTW-')
#     smtpObj.sendmail(sender, recipient, body)
#     smtpObj.quit()
#     print(f'Email sent to {recipient}')






def run_aligner(aligner, reference_genome, input_directory, output_name, input_zipped=None, email=None, threads=8):
    # User choices

    # Capture reference genome path 
    index_reference_dict = {
        'rat_w_1.0': '/home/gioele/bwaIndices/TempO-Seq_Rat_Whole_Transcriptome_1.0.fa',
        'human_s1500_1.2':'/home/gioele/bwaIndices/TempO-Seq_Human_S1500+_Surrogate_1.2.fa',
        'human_w_2.0': '/home/gioele/bwaIndices/TempO-Seq_Human_Whole_Transcriptome_2.0_ManifestV2.fa',
        'human_w_2.1': '/home/gioele/bwaIndices/TempO-Seq_Human_Whole_Transcriptome_2.1.fa',
        'mouse_s1500_1.2':'/home/gioele/bwaIndices/TempO-Seq_Mouse_S1500+_Surrogate_1.2.fa',
        'mouse_w_1.0': '/home/gioele/bwaIndices/TempO-Seq_Mouse_Whole_Transcriptome_1.0.fa'
    } 

    #This is derived from the reference genome
    if reference_genome in index_reference_dict.keys():
        reference_index = index_reference_dict[reference_genome]
    else: 
        reference_index = reference_genome
  
    # Are the input fastq zipped?
    zipped = input_zipped or True

    # Capture environment variables:
    # !!! if you move this getcwd, change the align() fn to exit the temporary directory every time as this will change otherwise and nest directories into each others 
    current_directory = input_directory or os.getcwd()
    file_list = os.listdir(current_directory)

    # Get the list of files in the directory (either .gz if zipped or .fastq if not zipped!)
    expected_extension = 'gz' if zipped else 'fastq'
    # Get the file list to process discriminating based on extension
    file_list = [file for file in file_list if file.split('.')[-1] == expected_extension]
    
    #file_list =["first.fastq.gz"]#, "second.fastq.gz", "third.fastq.gz"] 

    print('The following files will be included in the analysis: ', ', '.join(file_list))


    temp_dir_list = []
    # For each file run alignment
    if aligner == 'star':
        for file in file_list:
            temp_dir_list = align_star(file, reference_index, threads, zipped, temp_dir_list, current_directory)
    elif aligner == 'bwa':
        for file in file_list:
            temp_dir_list = align_bwa(file, reference_index, threads, zipped, temp_dir_list, current_directory)
    elif aligner == 'kallisto':
        for file in file_list:
            temp_dir_list = align_kallisto(file, reference_index, threads, zipped, temp_dir_list, current_directory)
    else:
        raise Exception('Aligner not recognised')

    # temp_dir_list = ['temp_first', 'temp_second', 'temp_third']


    # Compile all read counts
    print('Starting compiling read count table')
    # For each temp folder, load df and append it to general one
    # declare empty df to append to
    append_to = pd.DataFrame()
    # For each temporary directory, create df and append cols of interest
    for dir in temp_dir_list:
        # if the aligner is kallisto, no need to skip first row, but need to skip last col
        if aligner == 'kallisto':
            df = pd.read_csv(f'{current_directory}/{dir}/abundance.tsv', delimiter='\t')
            #remove last col + change first col name to Geneid
            df = df.iloc[:,:-1]
            df = df.rename(columns={'target_id': 'Geneid'})
            append_to = append_to_df(append_to, df, dir.removeprefix('temp_')) # Kallisto produces files differently from featureCount, so you need to adapt them to be the same
        else:
            df = pd.read_csv(f'{current_directory}/{dir}/read_count.txt', delimiter='\t', skiprows=[0])
            append_to = append_to_df(append_to, df)

    #Write created df to file
    outstring = f'{current_directory}/{output_name}_count_table.csv'
    append_to.to_csv(outstring, index=False)

    print('Wrote count table to ', outstring)

    # Remove all temp folders 
    for dir in temp_dir_list:
        # Do it in 2 steps just in case
        subprocess.Popen(f'rm {current_directory}/{dir}/{"read_count.txt" if aligner != "kallisto" else "abundance.tsv"}', shell=True).wait()
        subprocess.Popen(f'rmdir {current_directory}/{dir}/', shell=True).wait()
        print('Removed temporary directory ', dir)
        
    print('Processes completed')

    # # Send email to the user!
    # if email != None:
    #     send_email(email)
        

if __name__ == '__main__':
    # Reference index path
    # How many threads?
    # threads = 8
    # # User email
    # email = 'gioelemook97@gmail.com'
    # aligner = 'kallisto' #'star'
    # reference_genome = 'rat'
    # input_directory = '/home/gioele/rat_fasta/loop'
    # input_zipped = None
    # #This is derived from the input 
    # output_name = 'kallisto'
    #Call aligner

    parser = argparse.ArgumentParser(description="Tool to align your TempoSeq sequences and extract a feature count table. It allows to alignment with star, bwa and kallisto")
    parser.add_argument('-i', '--input-directory', required=True, help='Directory where the fastq files are. Input "." for current directory, or the full directory path')
    parser.add_argument('-a', '--aligner', required=True, help='Select the required aligner, options: star, bwa, kallisto')
    parser.add_argument('-g', '--reference-genome', required=True, help='Select the required reference genome, options: rat_w_1.0, human_s1500_1.2, human_w_2.0, human_w_2.1, mouse_s1500_1.2, mouse_w_1.0, or give a path to .fa file') # If you add here add to index_reference_dict
    parser.add_argument('-o', '--output-name', required=True, help='Prefix for output name')
    parser.add_argument('-t', '--threads', required=False, help='Number of thread used, default: 8')
    parser.add_argument('-u', '--unzipped', required=False, action='store_false', help='The input files are unzipped, default: False, default expexted extension .fastq.gz')
    args = vars(parser.parse_args())
    input_directory = args['input_directory']
    aligner = args['aligner']
    reference_genome = args['reference_genome']
    output_name = args['output_name']

    threads = args['threads'] or 8
    input_zipped = args['unzipped']

    if input_directory == '.':
        input_directory = os.getcwd()
    
    print(f'Options: input directory {input_directory}, aligner {aligner}, reference genome {reference_genome}, output name {output_name}, threads {threads}, input zipped {input_zipped}')
    run_aligner(aligner, reference_genome, input_directory, output_name, input_zipped=input_zipped, threads=threads)


