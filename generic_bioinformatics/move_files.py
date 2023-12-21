import subprocess 
import os
import sys

file_list = sys.argv[1]

with open(file_list, 'r') as f:
	files = f.readlines()

files = [x.strip() for x in files]

# get list of actual files and map them
actual_files = os.listdir()

mapper = {}


for file in files:
	actual_file = [x for x in actual_files if x.startswith(file)][0]
	mapper[file] = actual_file
	print(f"found {file} in {actual_file}")
	subprocess.check_call(f'mv {actual_file} deliver', shell=True)
	print('moved')



