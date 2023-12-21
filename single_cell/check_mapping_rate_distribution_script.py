import pandas as pd
import seaborn as sns
import os

# Variables
threshold = '10000'
indexes = ['S01_S1', 'S02_S2', 'S03_S3', 'S04_S4', 'S05_S5', 'S06_S6', 'S07_S7', 'S08_S8', 'S09_S9', 'S10_S10', 'S11_S11', 'S12_S12', 'S13_S13', 'S14_S14', 'S15_S15', 'S16_S16', 'S17_S17', 'S18_S18', 'S19_S19', 'S20_S20', 'S21_S21', 'S22_S22', 'S23_S23', 'S24_S24', 'S25_S25', 'S26_S26']

folders = [f'{threshold}_{x}' for x in indexes]

outdir = 'distribution_plots'


for folder in folders:
    # Load in barcode quantification and make it a dict
    barcodes = pd.read_csv(f'{folder}/barcode_quantification.csv')

    # Make the barcode the index so to have it as key, drop that col
    barcodes.index = barcodes['Complete_barcode']

    barcodes = barcodes[['Number_of_reads']]

    # Barcodes as dict
    barcodes_dict = barcodes.to_dict()['Number_of_reads']

    del barcodes # rm initial df

    # Load in the reads count table
    reads_df = pd.read_csv(f'{folder}/human_read_count.csv')

    # Make Geneid the index before transposing, drop that col
    ids = list(reads_df['Geneid'])

    reads_df.index = ids

    reads_df = reads_df.drop('Geneid', axis=1)

    # Transpose DF and make a column with all counts summed up
    reads_df = reads_df.T


    # Make a column with the sum of counts
    reads_df['Sum'] = reads_df.sum(axis=1)

    # Remove all the other columns
    reads_df = reads_df['Sum']

    # Get a dictionary out of the df
    reads_dict = reads_df.to_dict()

    del reads_df


    # Get a dictionary of how much is the % of mapping rate
    mapping_dict = {}

    for sequence, total_count in reads_dict.items():
        # Put in the dict the % of each cell's mapping
        mapping_dict[folder+sequence] = round(total_count/barcodes_dict[sequence]*100, 2)

    all_values = mapping_dict.values()


    fig=sns.displot(all_values)


    if not os.path.exists(outdir):
        os.mkdir(outdir)

    fig.savefig(f'{outdir}/{folder}_distribution.png')














