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
import os
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

#%%
## PARAMETERS
path_to_files = args.target_dir
outdir = args.output_dir
label = args.label

#path_to_files = "/home/filippo/Documents/freeclimb/VariantCalling/mock_data"
#outdir = "/home/filippo/Documents/freeclimb/VariantCalling/Config"
#label = "BGI" ## identifier of the sequencing batch (e.g. IGA, BGI etc.)
r1_identifier = '_1.fq.gz'
r2_identifier = '_2.fq.gz'

#%% read files from folder
p = Path(path_to_files)

## create list of sample folders (one subfolder per sample)
sample_folders = []
for x in p.iterdir():
    if x.is_dir():
        sample_folders.append(x.name)

## read sample file names
sample_files = []
for smp in range(len(sample_folders)):
    sample = sample_folders[smp]
    p = Path(path_to_files).joinpath(sample)
    for x in p.iterdir():
        print('sample: ',sample,'file: ',x.name)
        temp_dict = {'sample':sample}
        if r1_identifier in x.name:
            temp_dict['fastq_1'] = x
        if r2_identifier in x.name:
            temp_dict['fastq_2'] = x
        #temp_dict['file'] = x
        sample_files.append(temp_dict)
    
print('n. of samples from sample_folders:',len(sample_folders))
print('n. of records in sample_dict:',len(sample_files))

#%% convert list of dict to pandas dataframe
## add r1 and r2 file columns, remove NAs, take root file name
df = pd.DataFrame(sample_files)
df_long = pd.melt(df, id_vars=['sample'], value_vars=['fastq_1','fastq_2'], 
                  var_name = 'paired-end_file', value_name='file_path', ignore_index=True)
df_long = df_long.dropna()
df_long['root_name'] = [re.sub("\\..*$","",os.path.basename(x)) for x in df_long['file_path']]

#%% now convert from long to wide format
df = pd.pivot_table(df_long, index = ['root_name','sample'], columns = ['paired-end_file'], 
                         values = ['file_path'], aggfunc=lambda x: x) #Reshape from long to wide
df.index.name = None

#%% Normalize column names
df.reset_index(inplace=True)
df.columns = [x[1] if x[1] != '' else x[0] for x in df.columns]

#%% save out file
fname = Path(outdir).joinpath(label + '_sample_sheet.csv')
print('writing file names to {}'.format(fname))
#with open(fname, 'w', newline='') as output_file:
#    dict_writer = csv.DictWriter(output_file, keys)
#    dict_writer.writeheader()
#    dict_writer.writerows(sample_dict)

## converting to Pandas dataframe and writing out
#df = pd.DataFrame(sample_dict)
cols = ['sample', 'fastq_1', 'fastq_2']
df = df[cols]
df.to_csv(fname, index=False)

print("DONE!")

