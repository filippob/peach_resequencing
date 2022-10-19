
## load libraries
library("tidyverse")
library("data.table")

### SET UP ###
basefolder = "~/Documents/freeclimb/VariantCalling/"
sample_file = "Config/M07 04 Samples-Spreadsheet_v5 - UNIMI-DISAA.xlsx"
prefix = "1_ID2101_" ## beginning of file names in /home/ngs/freeclimb_resequencing/IGA_2022-01/delivery_20220110/raw_sequences/
r1_suffix = "_S1_L001_R1_001.fastq.gz"
r2_suffix = "_S1_L001_R2_001.fastq.gz"
server_path = "/home/ngs/freeclimb_resequencing/IGA_2022-01/delivery_20220110/raw_sequences"
outputdir = "/home/freeclimb/Config" ## on our server

### READ DATA ##
writeLines(" - reading sample sheet file")
fname = file.path(basefolder, sample_file)
samples = readxl::read_xlsx(fname, sheet = 1)
samples = samples |> 
  mutate(sample_name = paste(`IGA smp ID`, `Customer label`, `Plate #`, `Plate pos`, sep = "-"),
         sample = paste(`IGA smp ID`, `Customer label`, sep="-"))

## CREATE R1 and R2 file names
writeLines(" - creating absolute paths to R1 and R2 fastq files")
r1names = paste(prefix, samples$sample_name, r1_suffix, sep="")
r2names = paste(prefix, samples$sample_name, r2_suffix, sep="")

## create sample sheet
df <- data.frame("sample" = samples$sample,
           "fastq_1" = file.path(server_path, r1names),
           "fastq_2" = file.path(server_path, r2names))


## write out sample sheet
writeLines(" - writing out the sample sheet file")
fname = file.path(outputdir, "IGA_sample_sheet.csv")
print(paste("writing to ", fname))
fwrite(x = df, file = fname, sep = ",")

print("DONE!")