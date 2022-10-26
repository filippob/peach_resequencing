
import re
import csv
from pathlib import Path

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

print('n. of records in sample_dict:',len(sample_dict))
print(sample_dict[0])


keys = sample_dict[0].keys()

fname = Path(outdir).joinpath('prova.csv')
with open(fname, 'w', newline='') as output_file:
    dict_writer = csv.DictWriter(output_file, keys)
    dict_writer.writeheader()
    dict_writer.writerows(sample_dict)

