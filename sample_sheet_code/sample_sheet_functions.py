########################## Functions sections ##########################
from dictionaries import *
import pandas as pd
##################### Make template

# 24 well plate and
# 48 well plate
# Function to make 24 and 48 well plate sample sheet
def make_24_48_well_sample_sheet(filename, twentyfour):
    # Check if filename already exists!
    if os.path.isfile(filename):
        raise Exception('This file already exists') # Make this so the user can decide to overwrite
    
    #If true it's for the 24, if false it's for the 48
    if twentyfour:
        write = different_24_dict
    else:
        write = different_48_dict
    
    # Convert dictionary to pd data frame
    table_to_write = pd.DataFrame.from_dict(write)
    
    with open(filename, 'w') as file:
        table_to_write.to_csv(file, mode='w', index=True, header=True)

    return 'done'


# 96 well plate(s)
import os


# Function requires filename and the number of plates to be concatenated (1<x<8)
def make_96_well_sample_sheet(filename, n_plates):
    # Load table in df format
    table_96_well = pd.DataFrame.from_dict(different_96_dict)

    # n_plates has to be 1 < x < 8
    if n_plates < 1 or n_plates > 8:
        raise Exception('Cannot process that quantity')
    
    # Check if filename already exists!
    if os.path.isfile(filename):
        raise Exception('This file already exists') # Make this so the user can decide to overwrite

    # Open file once and append tables to it. The file is opened in 'w' so to over-write anything else (anyway you get blocked by the exception above if the file already exists)
    with open(filename, 'w') as file:
        for _ in range(n_plates):
            table_96_well.to_csv(file, mode='a', index=True, header=True)
    
    # In dash we return the file instance??
    return 'done'



####################################









##################### Import filled sample sheets

# Function recognising type of plate and if 96 well plate split in many DFs
#_ requires recognising the type of plate

# Check shape of DF


#


# Function to get 96-plates (even if many • doesn't enforce-check if it's doing it right)
def get_tables(sample_sheet):
    # Number of rows in sample sheet
    df_len = len(sample_sheet)
    #dfs to be
    dfs = []
    # Slice the sample sheet so that the recurring 10th line gets skipped (the empty one)
    for x in range(8, df_len+1, 10):
        low_end = 0 if x == 8 else x - 8
        dfs.append(sample_sheet.iloc[low_end:x, :])
    return dfs

# Function for enforcing naming restrictions on samples


####################################








##################### Prepare melted table

#







