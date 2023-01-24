########################## Functions sections ##########################
from dictionaries import *
import pandas as pd
##################### Make template

# 24 well plate
# upload everithn


# 48 well plate



# 96 well plate
import os

# Import dictionary only once so to don't have to re-process at each call
standard_96_well_table = pd.DataFrame.from_dict(standard_96_dict)
# Function requires filename and the number of plates to be concatenated (1<x<8)
def make_96_well_sample_sheet(filename, n_plates):
    # n_plates has to be 1 < x < 8
    if n_plates < 1 or n_plates > 8:
        raise Exception('Cannot process that quantity')
    
    # Check if filename already exists!
    if os.path.isfile(filename):
        raise Exception('This file already exists') # Make this so the user can decide to overwrite

    # Open file once and append tables to it. The file is opened in 'w' so to over-write anything else (anyway you get blocked by the exception above if the file already exists)
    with open(filename, 'w') as file:
        for _ in range(n_plates):
            standard_96_well_table.to_csv(file, mode='a', index=True, header=True)
    
    # In dash we return the file instance??




####################################









##################### Import filled sample sheets

# Function recognising type of plate and if 96 well plate split in many DFs
#_ requires recognising the type of plate





# Function for enforcing naming restrictions on samples


####################################








##################### Prepare melted table

#







