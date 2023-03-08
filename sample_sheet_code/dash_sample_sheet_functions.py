import pandas as pd
from io import StringIO
from base64 import b64decode
import re
from dictionaries import *

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
def check_table(key, table):
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

    # Check if there is something wrong (eg >1 well to change, or sample names are not unique)
    something_wrong = len(to_change) > 0 or total_samples != unique_sample_names

    # Check what type of plate it is based on key
    plate_type = key.split('_')[-1]

    # If the plate is 24/48 wells, double check if something is actually wrong (they have a lot of blanks)
    if plate_type in ['24', '48']:
        if plate_type == '24':
            something_wrong = len(to_change) > 0 or total_samples != unique_sample_names or blanks != 72
        if plate_type == '48':
            something_wrong = len(to_change) > 0 or total_samples != unique_sample_names or blanks != 48

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

    # Get from complete DF well name the plate flavor and well position
    p, flavour, well = complete.strip('_').split('_')
    plate = f'{p}_{flavour}'

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
    # match = False
    # name = None
    
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
    mapped_sheets = [map_sample_to_well(table, flavour) for flavour, table in tables.items()]

    # Concatenate sample sheets
    user_sample_sheet = stack_data_frames(mapped_sheets)

    # Now iterate through the large manifest/sample sheet data frame and substitute sample names with user's sample names if well and plate position matches
    complete_sample_sheet_to_return = migrate_sample_name(complete_sample_sheet, user_sample_sheet)

    # If it's requested the complete one return it, otherwise filter it and return it
    if complete:
        return complete_sample_sheet_to_return
    else:
        return filter_df(complete_sample_sheet_to_return)



# Helper function to manage the adding of the 24/48 wells plates to the main sample sheet
def add_tf_fe(sample_sheet, tuple, reversed):
    flavour, table = tuple

    # Map 96 well plate to sample sheet 
    mapped_plate = map_sample_to_well(table, flavour)
    
    # Remove blanks from mapping
    clean_plate = mapped_plate.drop(mapped_plate[mapped_plate['Sample_Name'] == 'blank'].index).reset_index(drop=True)

    # Map plate to indexes and get all the generic info needed
    dictionary = flavours_list[flavour]
    
    to_add = ["I7_Index_ID","index","I5_Index_ID","index2","index2_reverse"]
    for key in to_add:
        clean_plate[key] = clean_plate['Sample_Well'].map(dictionary[key])

    # If reversed drop the column index2 and rename index2_reverse to index2
    # Else drop index2_reverse

    if reversed:
        clean_plate = clean_plate.drop(columns=['index2'])
        clean_plate = clean_plate.rename(columns={'index2_reverse': 'index2'})

    else:
        clean_plate.drop(columns=['index2_reverse'], inplace=True)

    # Declare missing columns empty eg Sample_ID, Sample_Project and Description
    clean_plate['Sample_ID'], clean_plate['Sample_Project'], clean_plate['Description'], = '','','no_description'

    # Order columns as per sample sheet
    clean_plate = clean_plate[["Sample_ID","Sample_Name","Sample_Plate","Sample_Well","I7_Index_ID","index","I5_Index_ID","index2","Sample_Project","Description"]]

    # Merge DFs and reset indexes
    merged = pd.concat([sample_sheet, clean_plate])
    merged.reset_index(drop=True, inplace=True)
    merged["Sample_ID"] = merged.index + 1

    return merged


# Function to check the presence of duplicates among all tables
def check_no_duplicates_or_unwanted_chars(tables):
    tables_cells = [list(table.values.ravel()) for key, table in tables.items()]
    all_cells = []
    for cells in tables_cells:
        all_cells.extend(iter(cells))

    # No unwanted cells to start with
    unwanted = False
    # If there are no duplicates, check that there are no unwanted characters
    for cell in list(all_cells):
        # Check for regex (checks for non-word characters)
        match = re.findall(r'[^a-zA-Z0-9_-]', cell)
        if match:
            # If unwanted is still false, make it a dict
            if not unwanted:
                unwanted = {'unwanted': [cell]}
            # If already is a dict, append the unwanted cell
            else:
                unwanted['unwanted'].append(cell)
    # If any unwanted cell, return it
    if unwanted:
        return unwanted

    # If the length of the list is not the same as the set there are duplicates
    duplicates = len(all_cells) != len(set(all_cells))

    # If there are duplicates, remove single-occurrence cells and return duplicates
    if duplicates:
        for x in set(all_cells):
            all_cells.remove(x)

        # Key value pair the repeats with the occurrences
        occurrences = {x: all_cells.count(x)+1 for x in set(all_cells)}

        occurrences.pop('blank', 'none')

        # Return duplicate cells
        return occurrences

    # If there are unwanted cells return it, else return false
    return False




##########################
######## WRAPPER FUNCTIONS


# Wrapper to produce the template based on selection, and metadata, returns filename and sample sheet to be written 
def wrapper_make_template(selection, project_name, experiment_name, additional_comments, file_format='csv'):
    """
    Wrapper to produce the template based on selection, and metadata, returns filename and sample sheet to be written 
    """

    sample_sheet = ''
    filename = ''

    # Handle empty filename
    project_name = project_name if project_name != '' else 'unnamed_project'
        

    # Iterate through the selection, if a key is true, add it to the sample sheet.
    # If the key is tf (24) or fe (48), handle it accordingle

    for flavour, value in selection.items():
        # If the selection is false, continue
        if not value:
            continue
        # If the key is 24/48 act accordingly
        if flavour == 'Plate_24':
            import_df = pd.DataFrame.from_dict(different_24_dict)
            string = import_df.to_csv(index=False)
            
        elif flavour == 'Plate_48':
            import_df = pd.DataFrame.from_dict(different_48_dict)
            string = import_df.to_csv(index=False)
        
        else:
            # Standard 96 well plate, add flavour to it
            import_df = pd.DataFrame.from_dict(different_96_dict)
            new_cols = [x if not x.startswith('Unnamed') else flavour for x in import_df.columns]
            import_df.columns = new_cols

            string = import_df.to_csv(index=False)

        #Append string to sample sheet collective string + new line
        sample_sheet += string + ',,,,,,,,,,,,\n'


    # Add metadata to sample sheet
    metadata = f"""Project Name,{project_name},,,,,,,,,,,\nExperiment Name,{experiment_name},,,,,,,,,,,\nComments,{additional_comments},,,,,,,,,,,\n"""
    sample_sheet = metadata + sample_sheet

    project_name = project_name.replace(' ', '_')
    filename = f'template_{project_name}.{file_format}'

    return filename, sample_sheet




# Obtain information from the upload file and extract tables, reports and metadata
def wrapper_handle_upload(content, filename):
    """Obtain information from the upload file and extract tables, reports and metadata. This function also handles the decoding of the file content."""    
    
    # Handle file
    # content_type, content_string = content.split(',')
    # decoded_content = b64decode(content_string)
    # file_content = StringIO(decoded_content.decode('utf-8')).readlines()

    # Only for debug
    file_content = content


    # Slice metadata eg first 3 lines
    metadata = file_content[:3]
    # convert metadata
    metadata_dict = get_metadata(metadata)

    
    # Get rest of file, read it in with pandas and get the tables
    file_string = ''.join(file_content[3:])
    # Read file in with pandas
    df = pd.read_csv(StringIO(file_string), header=None)
    # Helper function to extract the tables and return it in dictionary where flavour:df
    tables = get_tables(df)

    # Perform checks on tables
    checks_tables = {key:check_table(key, value) for key, value in tables.items()}

    # Get general check report
    report = check_no_duplicates_or_unwanted_chars(tables)

    return tables, checks_tables, metadata_dict, report






# Wrapper create complete sample sheet!!
def wrapper_create_complete_sample_sheet(tables, complete, machinery, project, experiment, comments):
    # Declare list of machineries based on required I5 index
    possible_machinery_selection_forward_strand = ['mini_seq_standard', 'mi_seq', 'hi_seq_2500', 'hi_seq_2000','nova_seq_6000_v_1.5']
    possible_machinery_selection_reverse_strand = ['iSeq', 'next_seq_system', 'hi_seq_x', 'hi_seq_4000', 'hi_seq_3000', 'mini_seq', 'nova_seq_6000_v_1.0']

    #Determine if the reverse complement is required!
    if machinery in possible_machinery_selection_forward_strand:
        reverse_complement = False
    elif machinery in possible_machinery_selection_reverse_strand:
        reverse_complement = True
    else:
        raise Exception('Machinery not recognised!')
    

    # Handle the presence of 24/48 plates
    tf_or_fe = None
    # Retain 24/48 well layouts, if any, and produce the sample sheet for the rest before
    if 'Plate_24' in tables.keys():
        tf_or_fe = ('Plate_24', tables['Plate_24']) 
        del tables['Plate_24']
    if 'Plate_48' in tables.keys():
        tf_or_fe = ('Plate_48', tables['Plate_48']) 
        del tables['Plate_48']
    
    # Make header from standard one
    header = make_generic_header(project, experiment, comments, machinery)

    # Get body of sample sheet
    sample_sheet_body = get_sample_sheet_body(tables, complete)

    # Make filename string starting from project name
    # If the filename has spaces in it, replace it with underscores
    project_name = project if ' ' not in project else project.replace(' ', '_')
    filename = f'{project_name}_sample_sheet{"_all_primers" if complete else ""}.csv'
    
    # IF the machinery requires REVERSE I5, substitute it in the sample sheet!
    # (map index names to correct sequence from the dictionary)
    if reverse_complement:
        sample_sheet_body['index2'] = sample_sheet_body['I5_Index_ID'].map(f_primers_to_reverse_complement)

    # Append 24/48 plate to it, if any
    if tf_or_fe:
        sample_sheet_body = add_tf_fe(sample_sheet_body, tf_or_fe, reverse_complement)

    # At this point in the old versions, you'd reduce the number of cols for mini_seq, given that we're standardising we won't (won't break the pipeline)
    
    #Prepare text to return
    text = header + sample_sheet_body.to_csv(index=False, header=False)

    return filename, text
