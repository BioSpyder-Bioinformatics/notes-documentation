import os
import re

# Report
# Total number of files
# Total files with unique names
# How many folders created
# If all wells present eg A1 - H12


# Dry run, eg only get report
dry = True


# Standard BIOS, has to be the same
bios_number = None
correct_file_extension = 'fastq.gz'
# Current directory
dir = os.getcwd() + '/to_correct' # Correct this


# Declare report Dictionary
report_dict=dict(
    #Plates names
    plates = [],
    #Dictionary with list for each plate
    wells_listing = {
        #ex
        #  plate1 = ['A1', 'B5', 'C1']
    },
    unformatted_filenames = [],
    controls = [],
    out_of_ordinary = []
)


# Get a list of all folders and check if
# There are already directories, notify user if so
files = []
for (dirpath, dirnames, filenames) in os.walk(dir):
    files.extend(filenames)
    if dirnames:
        print('You are acting on a folder that has already subdirectories, do you want to continue?')
        choice = input('Y/n?')
        if choice == 'n' or choice == 'no':
            print('Stopping execution')
            raise SystemExit(0)
    print('These are your files: ', len(files))
    break



# There are BIOS which are different! In which case stop the iteration immediately
# This made as separate step to avoid interrupting mid copying
bios_experiments = []
for file in files:
    bios = re.findall(r'BIOS(\d{4})-.*.fastq.gz', file)
    if bios:
        bios_experiments.append(bios[0])

# Make set of Bios_experiments, if length is more than 1 interrupts iteration and lists which experiments are present
if len(set(bios_experiments)) > 1:
    print(f'There is more than one experiment in this folder, please double check! You have files for experiment BIOS:', [x for x in set(bios_experiments)])
    raise SystemExit()

bios_number = bios_experiments[0]



# Iterate through each file:
# - check that experiment is equal to the bios number (i know i checked already but just in case)
# - check that there is the folder for the plate number, if not create it
#      - Add it to report dictionary
# - check that the well does not appear in the file already
#      - Add in well listing with plate as key
# - copy file to folder
# - change file name removing bios
for file in files:
    try:
        # For each file try to extract the info from the filename, if it fails, add it to report and file to 'unformatted' directory
        experiment, plate, well = re.findall(r'BIOS(\d{4})-([A-Z]{2}\d{8})_([A-Z]\d{2}).fastq.gz', file)[0]
        # Check info were extracted correctly
        if experiment and plate and well:
            # If experiment is not the same one as defined before interrupt execution
            # Really doubt this will ever happen as they check pretty much the same regular expression, but just in case
            if experiment != bios_number:
                print(f'There is more than one experiment in this folder, please double check!')
                raise SystemExit()

            # Add file to report
            # Add plate
            report_dict['plates'].append(plate)
            # Add well
            if plate not in report_dict['wells_listing']:
                report_dict['wells_listing'][plate] = []
            report_dict['wells_listing'][plate].append(well)



            # Check that there is a folder with the plate number, if not create it
            for (dirpath, dirnames, filenames) in os.walk(dir):
                # If the plate folder does not exist make it
                if plate not in dirnames:
                    if not dry:
                        os.system(f'mkdir {dir}/{plate}')
                else:
                    #if it exists already just add the file and copy it to change the name
                    if not dry:
                        #Copy file
                        os.system(f'cp {file} {plate}/')
                        # Change filename
                        os.system(f'cp {plate}/{file} {plate}/{plate}_{well}.fastq.gz')
                break

            


        else:
            print('Something wrong with this file!')
            # Add filename to out_of_the_ordinary
            report_dict['out_of_the_ordinary'].append(file)
            # Add it to unformatted
    except:
        print(f'{file} did not match expected pattern! ')
        # Check if it contains bspc/bsnc and add it to the control subfolder, add it to report (ignorecase)
        is_control = re.findall(r'\w.*(bspc|bsnc).*.fastq.gz', file, re.I)
        if is_control:
            report_dict['controls'].append(file)
        # If not even a control, it's probably unformatted, so add it to the unformatted dictionary
        else:
            report_dict['unformatted_filenames'].append(file)




# Make the report plates a set
report_dict['plates'] = set(report_dict['plates'])

print(report_dict)


# Make folder with plate number name
# Copy file in that folder + remove BIOS#### in front








# Check that the BIOS is never different, if it is stop with exception!!






