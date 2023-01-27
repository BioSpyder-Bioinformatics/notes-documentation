########################## Functions sections ##########################
import os
from dictionaries import *
import pandas as pd
from datetime import datetime
import numpy as np
import re
##################### Make template

# 24 well plate and
# 48 well plate
# Function to make 24 and 48 well plate sample sheet
def make_24_48_well_sample_sheet(filename, twentyfour, is_csv=True):
    
    #If true it's for the 24, if false it's for the 48
    if twentyfour:
        write = different_24_dict
    else:
        write = different_48_dict
    
    # Convert dictionary to pd data frame
    table_to_write = pd.DataFrame.from_dict(write)
    
    string = table_to_write.to_csv()

    return string

    # with open(filename, 'w') as file:
    #     table_to_write.to_csv(file, mode='w', index=True, header=True)

    # return 'done'


# 96 well plate(s)


# Function requires filename and the number of plates to be concatenated (1<x<8)
def make_96_well_sample_sheet(filename, n_plates):
    # Load table in df format
    table_96_well = pd.DataFrame.from_dict(different_96_dict)

    # n_plates has to be 1 < x < 8
    if n_plates < 1 or n_plates > 8:
        raise Exception('Cannot process that quantity')

    #Return string for ease
    string = ''
    for _ in range(n_plates):
        string += table_96_well.to_csv() + ',,,,,,,,,,,,\n'

    return string
    
    #Old method, best returning the string
    # # Open file once and append tables to it. The file is opened in 'w' so to over-write anything else (anyway you get blocked by the exception above if the file already exists)
    # with open(filename, 'w') as file:
    #     for _ in range(n_plates):
    #         table_96_well.to_csv(file, mode='a', index=True, header=True)
    
    # # In dash we return the file instance??
    # return 'done'



####################################









##################### Import filled sample sheets

# Function recognising type of plate and if 96 well plate split in many DFs
#_ requires recognising the type of plate

# Check shape of DF


#


# Function to get 96-plates (even if many • doesn't enforce-check if it's doing it right)
# Need to read in the sample sheet before
# 
# sheet = pd.read_csv('file_name.csv', index_col=0)
# OORRRRR
#sheet = pd.read_excel('filename.csv', index_col=0)

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








##################### Prepare melted table (without header, one DF at the time)
def map_plate(to_map, flavour, next_seq=False, mini_seq=False):
    # Mapped_plate is a constant, flavour changes for 24 and 48 to 24/48-Well
    mapped_plate = pd.DataFrame.from_dict(plate_wells_mapping)

    df = pd.DataFrame({'Sample_Name':to_map.values.ravel(), 'Sample_Well': mapped_plate.values.ravel()})
    df['Sample_Plate'] = flavour 

    #OTHER CHECKS?? IF BLANK REMOVE? !!!!!!!!!!!!!!!!!!!!!!


    # Add indexes from dictionaries
    indexes_flavours = flavours_list[flavour]
    df['I5_Index_ID'] = df['Sample_Well'].apply(lambda x: indexes_flavours[x][0])
    df['I7_Index_ID'] = df['Sample_Well'].apply(lambda x: indexes_flavours[x][1])

    #Map sequences from index name
    # df['I5_Index_ID'] = df['Sample_Well'].map(indexes_flavours)

    df['index'] = df['I7_Index_ID'].map(r_primers_to_sequence)

    #If next seq map to other dict
    if next_seq:
        df['index2'] = df['I5_Index_ID'].map(f_primers_to_reverse_complement)
    else:
        df['index2'] = df['I5_Index_ID'].map(f_primers_to_sequence)

    # Add Sample_Project and Description rows (empty rows)
    df['Sample_Project'], df['Description']  = np.nan, np.nan
    
    # Add Sample_ID (eg row index, but with index 1 instead of 0)
    df['Sample_ID'] = df.index +1

    if mini_seq:
        df = df[['Sample_ID', 'Description', 'I7_Index_ID', 'index', 'I5_Index_ID', 'index2', 'Sample_Project']]
    else:
        df = df[['Sample_ID', 'Sample_Name', 'Sample_Plate', 'Sample_Well', 'I7_Index_ID', 'index', 'I5_Index_ID', 'index2', 'Sample_Project', 'Description']]

    return df

#map_plate(to_map, flavour) # Add kwargs next_seq/mini_seq if required


def loop_and_stack(tables, flavours, **kwargs):
    dfs = []
    for table, flavour in zip(tables, flavours):
        dfs.append(map_plate(table, flavour, **kwargs))
    return dfs




def concatenate_and_reindex(tables):
    merged = pd.concat([table for table in tables]).reset_index(drop=True)
    merged['Sample_ID'] = merged.index + 1
    return merged



def check_table(table):
    cells = list(table.values.ravel())
    blanks = 0
    pos_controls = 0
    neg_controls = 0
    # Regular expression extracts all cells that start with one or more underscores and are not positive or negative controls
    to_change = []
    for cell in list(cells):
        # Check for regex
        match = re.findall(r'(_+(?!bsnc|bspc)\w+)', cell)
        if match:
            cells.remove(cell)
            to_change.append(match[0])

        #Check if blank
        if cell == 'blank':
            cells.remove(cell)
            blanks +=1 

        #Check if pos control
        if cell == '_bspc':
            cells.remove(cell)
            pos_controls +=1 

        #Check if neg control
        if cell == '_bsnc':
            cells.remove(cell)
            neg_controls +=1 

    total_samples=len(cells)
    unique_sample_names=len(set(cells))
    something_wrong = len(to_change) > 0 or total_samples != unique_sample_names
    
    specs = dict(blanks = blanks, positive_controls = pos_controls, negative_controls=neg_controls, to_change = to_change, total_samples=total_samples, unique_sample_names=unique_sample_names, something_wrong=something_wrong)
    return specs












######## For new solution

def map_sample_to_well(table, flavour):
    mapped_plate = pd.DataFrame.from_dict(plate_wells_mapping)

    df = pd.DataFrame({'Sample_Name':table.values.ravel(), 'Sample_Well': mapped_plate.values.ravel()})
    df['Sample_Plate'] = flavour

    return df



def stack_data_frames(dfs):
    return pd.concat([df for df in dfs]).reset_index(drop=True)




def compare_rows(complete, user):
    traslated_dict = user.T.to_dict()
    traslated_arr_of_dict = [traslated_dict[x] for x in range(len(traslated_dict))]

    # Get from complete DF well name the plate flavor and well position (_ is placeholder)
    _, plate, well = complete.strip('_').split('_')

    #Make case 
    match = False
    name = None

    for row in traslated_arr_of_dict:
        if row['Sample_Well'] == well and row['Sample_Plate'] == plate:
            match = True
            name = row['Sample_Name']
       
    return name if match else complete




def migrate_sample_name(complete, user):
    match = False
    name = None
    
    complete['Sample_Name'] = complete['Sample_Name'].apply(lambda x: compare_rows(x, user))
    #pd.DataFrame.apply(compare_rows(complete.iterrows(), user))
    return complete



def filter_df(df):
    df2 = df[~df['Sample_Name'].str.startswith('__')]
    return df2







