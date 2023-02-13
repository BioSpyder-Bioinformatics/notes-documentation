import pandas as pd
from io import StringIO
from base64 import b64decode
import re
from dictionaries import plate_wells_mapping, complete_sample_sheet_string

####### HELPER FUNCTIONS

# Get metadata on reupload
def get_metadata(lines):
    metadata = {}
    for line in lines:
        cells = line.split(',')
        key = cells[0].replace(' ', '_').lower()
        metadata[key] = cells[1]

    return metadata


# Get tables from DF
def get_tables(sample_sheet):
    """Returns plates in format [{plate_selection:plate(df)}, ...]"""
    # Number of rows in sample sheet
    df_len = len(sample_sheet)
    #dfs to be
    dfs = {}
    # Slice the sample sheet so that the recurring 10th line gets skipped (the empty one)
    for x in range(9, df_len+1, 10):
        low_end = 0 if x == 9 else x - 9
        # Grab the table based on indexes
        new_table = sample_sheet.iloc[low_end:x, :]
        
        # Make first row equal to header
        new_header = new_table.iloc[0] #grab the first row for the header
        new_table = new_table[1:] #take the data less the header row
        new_table.columns = new_header #set the header row as the df header

        plate_flavour = new_header[0]

        # Make first column equal to index
        # Get first column
        new_index = list(new_table.pop(plate_flavour))
        # Make that col the new index
        new_table.index = new_index

        #dfs.append({plate_flavour:new_table})
        dfs[plate_flavour] = new_table

    return dfs



# Checks if the sample names inputed by user are in line with the suggested naming convention
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

    # Check what type of plate it is and return it in the dict!
    plate_type = ''
    
    specs = dict(blanks = blanks, positive_controls = pos_controls, negative_controls=neg_controls, to_change = to_change, total_samples=total_samples, unique_sample_names=unique_sample_names, something_wrong=something_wrong)
    return specs



# Map sample to well
# Creates mapping of user sample names to the designated wells
def map_sample_to_well(table, flavour):
    mapped_plate = pd.DataFrame.from_dict(plate_wells_mapping)

    df = pd.DataFrame({'Sample_Name':table.values.ravel(), 'Sample_Well': mapped_plate.values.ravel()})
    df['Sample_Plate'] = flavour

    return df


# Stack data frames
# Simply concatenates DFs and resets indexes
def stack_data_frames(dfs):
    return pd.concat([df for df in dfs]).reset_index(drop=True)


# Migrate sample name + helper
# Helper function supporting migrate_sample_name(), compares well position
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

# Migrate sample name
# Map user's sample names to the large pre-built sample sheet
def migrate_sample_name(complete, user):
    match = False
    name = None
    
    complete['Sample_Name'] = complete['Sample_Name'].apply(lambda x: compare_rows(x, user))
    #pd.DataFrame.apply(compare_rows(complete.iterrows(), user))
    return complete



# Remove rows that start with double underscore from complete sample sheet DF
def filter_df(df):
    df2 = df[~df['Sample_Name'].str.startswith('__')]
    df2 = df2.reset_index(drop=True)
    df2['Sample_ID'] = df2.index +1
    return df2




# Get body of sample sheet (without 24/48 handling)
#Wrapper get body of sample sheet
def get_sample_sheet_body(tables, complete):
    # Get complete sample sheet from string
    complete_sample_sheet = pd.read_csv(StringIO(complete_sample_sheet_string), index_col=0)

    # Map sample names based on well position
    mapped_sheets = [map_sample_to_well(table, flavour.capitalize()) for flavour, table in tables.items()]

    # Concatenate sample sheets
    user_sample_sheet = stack_data_frames(mapped_sheets)

    # Now iterate through the large manifest/sample sheet data frame and substitute sample names with user's sample names if well and plate position matches
    complete_sample_sheet_to_return = migrate_sample_name(complete_sample_sheet, user_sample_sheet)

    # If it's requested the complete one return it, otherwise filter it and return it
    if complete:
        return complete_sample_sheet_to_return
    else:
        return filter_df(complete_sample_sheet_to_return)







##########################
######## WRAPPER FUNCTIONS
