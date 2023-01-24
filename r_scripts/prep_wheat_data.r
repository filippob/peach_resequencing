
## script to prepare example data

library("BGLR")
library("data.table")

base_folder = '~/Documents/freeclimb/g_x_e'
outdir = "data"


## get sample data
data(wheat)
Y = wheat.Y # grain yield evaluated in 4 different environments (599 samples, 4 columns)
X = wheat.X ## SNP genotypes
A = wheat.A
 
print("Writing out example data")
dir.create(file.path(base_folder, outdir), recursive = TRUE, showWarnings = FALSE)

writeLines(" - writing out phenotype data")
fname = file.path(base_folder, outdir, "phenotypes.csv")
fwrite(x = Y, file = fname, sep = ",")

writeLines(" - writing out marker data")
fname = file.path(base_folder, outdir, "markers.csv")
fwrite(x = X, file = fname, sep = ",")

writeLines(" - writing out pedigree relationship data")
fname = file.path(base_folder, outdir, "Amat.csv")
fwrite(x = A, file = fname, sep = ",")

print("DONE!")