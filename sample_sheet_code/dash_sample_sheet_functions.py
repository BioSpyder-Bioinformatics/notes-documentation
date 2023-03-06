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
        print('reversed')
        clean_plate = clean_plate.drop(columns=['index2'])
        clean_plate = clean_plate.rename(columns={'index2_reverse': 'index2'})
        print(clean_plate.keys())

    else:
        clean_plate.drop(columns=['index2_reverse'], inplace=True)

    # Declare missing columns empty eg Sample_ID, Sample_Project and Description
    clean_plate['Sample_ID'], clean_plate['Sample_Project'], clean_plate['Description'], = '','',''

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
        # Check for regex
        match = re.findall(r'[^a-zA-Z0-9_-]', cell)#######################################################################################################################################################################################################################################################################################################
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

    # Iterate through the selection, if a key is true, add it to the sample sheet.
    # If the key is tf (24) or fe (48), handle it accordingle

    for flavour, value in selection.items():
        # If the selection is false, continue
        if not value:
            continue
        # If the key is 24/48 act accordingly
        if flavour == 'tf':
            import_df = pd.DataFrame.from_dict(different_24_dict)
            string = import_df.to_csv(index=False)
            
        elif flavour == 'fe':
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
    """Obtain information from the upload file and extract tables, reports and metadata"""

    #Checks before the function starts for filename and content
    if filename is not None and filename.split('.')[1] != 'csv':
        #return genericLabelWrapper("Please only upload 'csv' files"), None
        return "Please only upload 'csv' files", None, None
    if content is None:
        # return genericLabelWrapper('The content of the file you uploaded is not valid'), None
        return 'The content of the file you uploaded is not valid', None, None
        
    
    # Handle file
    # UNCOMMENT THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    checks_tables = {key:check_table(value) for key, value in tables.items()}

    # Get general check report
    report = check_no_duplicates_or_unwanted_chars(tables)

    return tables, checks_tables, metadata_dict, report






# Get the complete sample sheet
def wrapper_create_complete_sample_sheet(tables, complete, machinery, project, experiment, comments):
    # Handle the presence of 24/48 plates
    tf_or_fe = None
    # Need to flag if the reverse complement is required
    reverse_complement = False

    # Retain 24/48 well layouts and produce the sample sheet for the rest
    if 'tf' in tables.keys():
        tf_or_fe = ('tf', tables['tf']) 
        del tables['tf']
    if 'fe' in tables.keys():
        tf_or_fe = ('fe', tables['fe']) 
        del tables['fe']

    # Get body of sample sheet
    sample_sheet_body = get_sample_sheet_body(tables, complete)


    # Make header and adapt file based on machine chosen
    header = None
    filename=None
    # 
    if machinery == 'hi_seq':
        # The header is basically empty
        header = make_hiseq_header()
        filename = f'{project}_SR_HS{"_all_primers" if complete else ""}.csv'
        # File does not change
    #
    elif machinery == 'mini_seq':
        # Header only has experiment 
        header = make_miniseq_header(experiment)
        # Refactor the columns to reduce them to the required number
        sheet = sheet[['Sample_ID', 'Description', 'I7_Index_ID', 'index', 'I5_Index_ID', 'index2', 'Sample_Project']]
        filename = f'{project}_MiniSeq{"_all_primers" if complete else ""}.csv'
    #
    elif machinery == 'mi_seq':
        # Nothing special
        header = make_miseq_header(project, experiment, comments)
        filename = f'{project}_SR_MS{"_all_primers" if complete else ""}.csv'
    #
    elif machinery == 'next_seq':
        header = make_nextseq_header(project, experiment, comments)
        # Next seq requires reverse I5, map the index names to the correct sequence from the dictionary
        sample_sheet_body['index2'] = sample_sheet_body['I5_Index_ID'].map(f_primers_to_reverse_complement)
        reverse_complement = True
        filename = f'{project}_PE_NS{"_all_primers" if complete else ""}.csv'
    else:
        raise Exception('I do not recognise this machinery instrument!')

    # Append 24/48 plate to it
    if tf_or_fe:
        sample_sheet_body = add_tf_fe(sample_sheet_body, tf_or_fe, reverse_complement)

    text = header + sample_sheet_body.to_csv(index=False, header=False)

    return filename, text





