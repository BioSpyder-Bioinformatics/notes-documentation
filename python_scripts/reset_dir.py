import os

rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
cols = range(1, 6) #up to 13


_name = 'BIOS2444-TC00002695_P25.fastq.gz'

dir = os.getcwd()


# Remove content of directory

os.system(f'rm -r {dir}/to_correct/*')




for row in rows:
    for col in cols:
        # Missing file
        # if row == 'C' and col == 3:
        #   b continue
        name = f'BIOS2444-TC00002695_{row}{col:02d}.fastq.gz'
        os.system(f'echo "sample" > {dir}/to_correct/{name}')



spike = True

if spike:
    files = [
            'BIOS2444TC00002695_P05.fastq.gz', #Wrong format
            'BIOS2444-TC00002692_A03.fastq.gz', #Different plate
            #'BIOS2434-TC00002695_P25.fastq.gz', #Different experiment
            'BIOS2444-HLiver01_bspc.fastq.gz', # Control
            'receipt.txt', # Wrong extension 
            ]
    for name in files:
        os.system(f'echo "sample" > {dir}/to_correct/{name}')
        print('done')
