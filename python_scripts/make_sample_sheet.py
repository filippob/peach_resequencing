
import re
import csv
import pandas as pd
from pathlib import Path
from collections import OrderedDict

path_to_files = '/home/ngs/freeclimb_resequencing/BGI-FreeClimb'
outdir = '/home/freeclimb/Config'

#%%--
p = Path(path_to_files)

sample_folders = []
for x in p.iterdir():
    if x.is_dir():
        sample_folders.append(x.name)

sample_dict = []
for smp in range(len(sample_folders)):
    sample = sample_folders[smp]
    p = Path(path_to_files).joinpath(sample)
    temp_dict = {'sample':sample}
    for x in p.iterdir():
        print(x.name)
        if '_1.fq.gz' in x.name:
            temp_dict['fastq_1'] = x
        if '_2.fq.gz' in x.name:
            temp_dict['fastq_2'] = x

    sample_dict.append(temp_dict)
    #sample_dict.append(OrderedDict(sorted(temp_dict.items())))

print('n. of records in sample_dict:',len(sample_dict))

fname = Path(outdir).joinpath('prova.csv')
#with open(fname, 'w', newline='') as output_file:
#    dict_writer = csv.DictWriter(output_file, keys)
#    dict_writer.writeheader()
#    dict_writer.writerows(sample_dict)


df = pd.DataFrame(sample_dict)
cols = ['sample', 'fastq_1', 'fastq_2']
df = df[cols]
df.to_csv(fname, index=False)

