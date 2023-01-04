
## R script to plot results from admixture (from get_admixture_results.sh)

library("ggpubr")
library("ggplot2")
library("tidyverse")
library("data.table")

prj_folder = "/home/filippo/Documents/freeclimb/VariantCalling"
exp_folder = "Analysis/admixture"
outdir = "results"
cvf = file.path(exp_folder, "CVerr")
niterf = file.path(exp_folder, "niter")

CVerr = fread(cvf)
niter = fread(niterf)

CVerr$nK = as.integer(gsub("[^0-9]","",as.character(CVerr$V3)))
CVerr <- CVerr[order(CVerr$nK),]
CVerr <- rename(CVerr, cv_error = V4) |> select(cv_error,nK)

A <- cbind(CVerr,niter$V1)
A <- rename(A, niter = V2)
# 
# plot(A$nK,A$cv_error,type="l",main="Predictive accuracy",xlab="#n K", ylab="Cross-validation error")
# plot(A$nK,A$niter,type="l",main="Reaching convergence",xlab="#n K", ylab="#n of iterations")

p1 <- ggplot(A, aes(x = nK, y = cv_error)) + geom_line() + ggtitle("Predictive accuracy")
p2 <- ggplot(A, aes(x = nK, y = niter)) + geom_line() + ggtitle("Reaching convergence")

g <- ggarrange(p1,p2, ncol = 1)
fname = file.path(prj_folder, outdir, "cv_admixture.png")
ggsave(filename = fname, plot = g, device = "png", width = 6, height = 8)

####################
## R: plot results!!
####################

igaf = "Config/IGA_sample_sheet.csv"
bgif = "Config/BGI_sample_sheet_fix.csv"
reseqf = "Config/Reseq_sample_sheet.csv"
ncbif = "Config/ncbi_sample_sheet.csv"

### METADATA ###
iga = fread(igaf)
iga$dataset <- "IGA"
iga <- select(iga, -c(fastq_1,fastq_2))
bgi = fread(bgif)
bgi$dataset <- "BGI"
bgi <- select(bgi, -c(fastq_1,fastq_2))
bgi <- unique(bgi)
reseq = fread(reseqf)
reseq$dataset <- "Reseq"
reseq <- select(reseq, -c(fastq_1,fastq_2))
reseq <- unique(reseq)
ncbi = fread(ncbif)
ncbi$dataset <- "NCBI"
ncbi <- select(ncbi, -c(fastq_1,fastq_2))
ncbi <- unique(ncbi)

metadata = bind_rows(iga,bgi,reseq,ncbi)

### outputfiles
fileName <- "imputed"
fileType <- "Q"

QFam <-read.csv(file.path(exp_folder, paste(fileName,".fam",sep="")),sep=" ",header=FALSE)
QFam <- select(QFam, V1) |> rename(sample = V1)
QFam <- QFam |> inner_join(metadata, by = "sample")

#### sets the number of populations to through
#lowest K
Start_K <- 5
#highest K
End_K <- 15

library("RColorBrewer")
n <- Start_K
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))


all_data <- tibble(sample=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

for (k in Start_K:End_K) {
  data <- read_delim(file.path(exp_folder, paste(fileName,k,fileType,sep=".")),
                     col_names = paste0("Q",seq(1:k)),
                     delim=" ")
  data$sample <- QFam$sample
  data$dataset <- QFam$dataset
  data$k <- k
  
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k, -dataset) -> data
  all_data <- rbind(all_data,data)
}

all_data$Q <- factor(all_data$Q, levels = paste("Q", seq(1:15), sep=""))
all_data %>%
  filter(k == Start_K) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2)) +
  scale_fill_manual(values = sample(col_vector, Start_K), aesthetics = "fill") + 
  guides(fill=guide_legend(title="Q"))


p <- all_data %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2),
        axis.text.y = element_text(size = 5),
        strip.text = element_text(size = 7)) +
  scale_fill_manual(values = sample(col_vector, k), aesthetics = "fill") +
  guides(fill=guide_legend(title="Q")) + 
  facet_wrap(~k,ncol=1)

fname = file.path(prj_folder, outdir, "admixture_plot.png")
ggsave(filename = fname, plot = p, device = "png", width = 8.5, height = 10.5)
