
# INPUT CONFIGURATION MANAGEMENT ------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 1){
  #loading the parameters
  source(args[1])
} else {
  #this is the default configuration, used for development and debug
  writeLines('Using default config')
  
  #this dataframe should be always present in config files, and declared
  #as follows
  config = NULL
  config = rbind(config, data.frame(
    base_folder = '~/Documents/freeclimb',
    # phenotypes = 'GxE/merge_ID_phen_corrected/merged_ID_phen_corrected_edited.csv',
    phenotypes = 'g_x_e/data/merged_ID_phen_corrected_edited.csv',
    genotypes = 'VariantCalling/Analysis/merged_vcf/extract_imputed.raw',
    trait = "MD",
    year = 2018,
    outdir = 'GxE/data',
    prefix = "MD_2018",
    force_overwrite = FALSE
  ))
  
}


# SETUP -------------------------------------------------------------------
library("tidyverse")
library("data.table")


# READ THE DATA -------------------------------------------------------------------
writeLines(" - reading the phenotypic data")
fname = file.path(config$base_folder, config$phenotypes)
pheno = fread(fname)

vec = c("sample","Reference","loc","year","platform",config$trait)
pheno = pheno |> filter(year == config$year) |> select(all_of(vec))
# filter(pheno, is.na(!!as.symbol(config$trait)))
pheno = filter(pheno, !is.na(.data[[config$trait]]))

print(paste("From the phenotype file", fname, " ", nrow(pheno), "non-missing records have been read for", config$trait, "in year", config$year))

locations = unique(pheno$loc)
print(paste("Phenotypic data from the following locations:", paste(locations, collapse = ",")))
n_locations = length(locations)
stopifnot(n_locations > 1)

writeLines(" - reading the genotypic data")
fname = file.path(config$base_folder, config$genotypes)
geno = fread(fname, header = FALSE, skip = 1)
geno <- rename(geno, sample=V2)

# PREPROCESSING -------------------------------------------------------------------
writeLines(" - preprocessing")
print("select columns and check for duplicate IDs")
vec = c("sample","loc",config$trait)
pheno <- pheno |> select(all_of(vec))
names(pheno) = c("sample","loc","trait")
duplicates = group_by(pheno, sample) |> summarise(N = n()) |> filter(N>n_locations) |> pull(sample)

print(paste("N. of duplicate IDs removed", length(duplicates)))
pheno <- filter(pheno, !sample %in% duplicates)

writeLines(" - matching phenotype and genotype data")
samples <- unique(pheno$sample)
vv <- geno$sample %in% samples
geno = filter(geno, sample %in% samples)

print(paste("N. of samples matched between phenotypes and genotypes", sum(vv)))

## reshape phenotypic data and reorder to have same order as genotype data
pheno = pheno |> spread(key = loc, value = trait)
idx <- match(geno$sample, pheno$sample)
pheno <- pheno[idx,]

miss_rate = pheno |> select(-sample) |> is.na() |> colSums()/nrow(pheno)
print("missing rate per location:")
print(miss_rate)

## remove missing columns
writeLines(" - removing columns (locations) with missing rate > 50% - if any-.")
print(paste("removing the following location(s):",names(which(miss_rate>0.5))))
vec = miss_rate < 0.5
vec = c(TRUE,as.vector(vec))
pheno = as.data.frame(pheno)
pheno = pheno[,vec]

print(head(pheno))

stopifnot(sum(pheno$sample == geno$sample) == nrow(geno))

dirname = file.path(config$base_folder,config$outdir,config$prefix)
dir.create(dirname, showWarnings = FALSE)

writeLines(" - writing out files phenotypes.csv and markers.csv")
fname = file.path(dirname,"phenotypes.csv")
fwrite(select(pheno, -sample), file = fname)

fname = file.path(dirname,"markers.csv")
geno = geno[,-c(1:6)]
fwrite(geno, file = fname, col.names = TRUE, row.names = FALSE)

print("DONE!!")
