"""
script to create the sample sheet file from the content of the target parent directory;
the current version of the script looks into subfolders,
for a different structure of folders and files, the script needs to be revised
"""

"""
RUN AS:
python3 make_sample_sheet.py --target_dir=<path_to_files> --output_dir=<config_folder> --label=<sequencing_tag>
"""

## IMPORT LIBRARIES
import re
import csv
import argparse
import pandas as pd
from pathlib import Path
from collections import OrderedDict

# Create the parser
parser = argparse.ArgumentParser(description='Make the sample sheet file for the resequencing-mem pipeline')

# Add arguments
parser.add_argument('-t', '--target_dir', type=str, required=True, 
                    help='name of the target folder from where the R1 and R2 fastq files must be read')
parser.add_argument('-s', '--output_dir', type=str, required=True, 
                    help='directory where the sample sheet is to be stored (created if needed)')
parser.add_argument('-l', '--label', type=str, required=True, 
                    help='sequencing label (e.g. IGA, BGI, etc.)')

# Parse the argument
args = parser.parse_args()

# Print to check arguments values
print('Target folder is:', args.target_dir)
print('Output folder is:', args.output_dir)
print('Sequencing run label is:', args.label)

## PARAMETERS
path_to_files = args.target_dir
outdir = args.output_dir
label = args.label

#%% read files from folder
p = Path(path_to_files)

## create list of sample folders (one subfolder per sample)
sample_folders = []
for x in p.iterdir():
    if x.is_dir():
        sample_folders.append(x.name)

## read sample file names
sample_dict = []
for smp in range(len(sample_folders)):
    sample = sample_folders[smp]
    p = Path(path_to_files).joinpath(sample)
    temp_dict = {'sample':sample}
    for x in p.iterdir():
        print('sample: ',sample,'file: ',x.name)
        if '_1.fq.gz' in x.name:
            temp_dict['fastq_1'] = x
        if '_2.fq.gz' in x.name:
            temp_dict['fastq_2'] = x

    sample_dict.append(temp_dict)
    #sample_dict.append(OrderedDict(sorted(temp_dict.items())))

print('n. of records in sample_dict:',len(sample_dict))

fname = Path(outdir).joinpath(label + '_sample_sheet.csv')
print('writing file names to {}'.format(fname))
#with open(fname, 'w', newline='') as output_file:
#    dict_writer = csv.DictWriter(output_file, keys)
#    dict_writer.writeheader()
#    dict_writer.writerows(sample_dict)

## converting to Pandas dataframe and writing out
df = pd.DataFrame(sample_dict)
cols = ['sample', 'fastq_1', 'fastq_2']
df = df[cols]
df.to_csv(fname, index=False)

print("DONE!")

