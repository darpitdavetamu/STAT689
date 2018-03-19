# Extracting COAD data from TCGA using TCGA2Stat


# Set working directory to the same location as the R script. Adjust as necessary by adding commented directory
# setwd("~/Dropbox/My_Articles/ForStudents/ClassificationMultimodal/Datasets/TCGA-COAD")

# Package to extract TCGA data
library(TCGA2STAT)
# normalized count RNASeq 2 data, not quite counts
rnaseq.coad <- getTCGA(disease = "COAD", data.type = "RNASeq2", type = "RSEM")

# 20501 by 328
# normalized read counts miRNA data
mirna.coad <- getTCGA(disease = "COAD", data.type = "miRNASeq", type = "rpmmm")
# 705 by 221

#Log-transform rnaseq and mirna
rna_log <- log2(rnaseq.coad$dat + 1)
mirna_log <- log2(mirna.coad$dat + 1)

#The fourth series of digits (e.g. 01A or 11A) indicates the tissue type. 01A = primary colon adenocarcinoma; 11A = normal tissue. Only use samples from cancer tissue
############################################################
# This correspond to normal tissue, 41 such samples
# colnames(rna_log)[grep("-11A-", colnames(rna_log))]
#subset_normal <- (1:ncol(rna_log))[grep("-11A-", colnames(rna_log))]
subset_cancer <- (1:ncol(rna_log))[grep("-01A-", colnames(rna_log))]
#subject_ID_n <- substr(colnames(rna_log)[subset_normal],1,12)
subject_ID<- substr(colnames(rna_log)[subset_cancer],1,12)
rna_log <- rna_log[, subset_cancer]
subset_cancer_miRNA <- (1:ncol(mirna_log))[grep("-01A-", colnames(mirna_log))]
mirna_log <- mirna_log[, subset_cancer_miRNA]
subject_mirna <- substr(colnames(mirna_log),1,12)


# Filter the data based on zero or close to zero counts
zero_rna <- apply(rna_log, 1, function(x) sum(x<0.1))
# zero_mirna <- apply(mirna_log, 1, function(x) sum(x<0.01))
subset_rna <- (1:nrow(rna_log))[zero_rna < round(ncol(rna_log)*1/3)]
# subset_mirna <- (1:nrow(mirna_log))[zero_mirna < round(ncol(mirna_log)*1/3)]
rna_log <- rna_log[subset_rna, ]
# mirna_log <- mirna_log[subset_mirna, ]

#Further filter rna data based on sd
sd_rna <- apply(rna_log, 1, sd)
subset_rna <- (1:nrow(rna_log))[sd_rna > 0.5]
length(subset_rna)
rna_log <- rna_log[subset_rna,]

# Connect by platforms and subjects
#####################################
# Intersect mRNA and miRNA
subject_RNAs <- intersect(subject_ID, subject_mirna) 
length(subject_RNAs) # 218 which is very nice, essentially all for which I had mirna

# Union mRNA and miRNA
subject_union2 <- union(subject_ID, subject_mirna)
length(subject_union2)

#################
# Subtype information from consensus subtypes, reference paper Guinney et al., 2015
labels <- read.table("~/Dropbox/STAT 689 Project/Data/Log_trash/cms_labels_public_all.txt", header = T)
tcga <- labels[labels$dataset=="tcga", ]
labelpresent <- match(subject_union2, tcga$sample)
subtype <- tcga$CMS_network[labelpresent] #CMS subtypes based on core consensus samples
subtype[subtype == "UNK"] = NA # Treat unclassified as NA
subtype = as.factor(as.character(subtype))

# How many samples are with complete set of data (group + 2 views) versus incomplete

# Missing group, both views are present
sum(is.na(subtype[(subject_union2 %in% subject_mirna)]))

# Missing group, only RNAseq is present
sum(is.na(subtype[!(subject_union2 %in% subject_mirna)]))

# Available group, only RNAseq is present
sum(!is.na(subtype[!(subject_union2 %in% subject_mirna)]))

# Available group and both views
sum(!is.na(subtype[(subject_union2 %in% subject_mirna)]))

# For the 37 + 51 for which no subtype information, were some of them later classified into one of CMS subtypes? Yes, see below, this can be used for training
indextest = subject_union2[is.na(subtype)]
label_test = match(indextest, tcga$sample)
subtype_test <- tcga[label_test, 5]
sum(subtype_test %in% c("CMS1", "CMS2", "CMS3", "CMS4"))
sum(subtype_test == "NOLBL", na.rm = T)

# Want to do semi-supervised learning and matching of platforms
# Create a list matched by subject/group/subtype - for Yunfeng (?)
X2 = matrix(NA, ncol(rna_log), nrow(mirna_log))
X2[subject_union2 %in% subject_mirna,] = t(mirna_log)

COAD <- list(sample = subject_union2, subtype = subtype, subtypeRF = tcga[labelpresent, 5], RNA = t(rna_log), miRNA = X2)

save(COAD, file = "~/Documents/darpitdave/TCGA-COAD/Data/STAT689_COAD.Rda")




