from pathlib import Path

path_to_files = '/home/filippo/Documents/freeclimb'

#%%--
p = Path(path_to_files)

for x in p.iterdir():
    print(x)

