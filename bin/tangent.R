### Author: Wolu Chukwu, wchukwu@broadinstitute.org, Shahab Sarmashghi, ssarmash@broadinstitute.org
### Date last updated: April 18, 2025
### License: GNU GPL2, Copyright (C) 2024 Dana-Farber Cancer Institute

#################################################################
##                      READ IN LIBRARIES                      ##
#################################################################
library(useful)
library(tidyverse)
library(data.table)
library(pracma)
library(optparse)
library(rCGH)
library(dplyr)

##################################################################
##                DEFINE INPUT OPTIONS AND FLAGS                ##
##################################################################
option_list <- list(
  make_option(c("-n", "--pon"), type = "character", default = NA,
              help = "Path to panel of normals"),
  make_option(c("-t", "--tumor_cts"), type = "character", default = NA,
              help = "Path to tumor depth of coverage file"),
  make_option(c("-m", "--scale-mthd"), type = "character", default = "median",
              help = "Scaling method for signal intensities. Can be either mean or median. Defaults to median"),
  make_option(c("-r", "--ref_cli"), type = "character", default = NA,
              help = "Path to clinical data for the reference panel of normals"),
  make_option(c("-c", "--clin"), type = "character", action = "store_true",
              help = "Path to clinical data for any normals provided to supplement the reference panel"),
  make_option(c("-s", "--sample_sex"), type = "character", default = NA,
              help = "Sample's sex"),
  make_option(c("-d", "--n_dim"), type = "double", default = 10,
              help = "Number of dimension to use for SVD"),
  make_option(c("-f", "--fix_outliers"), type = "logical", default = FALSE,
              help = "Boolean indicating whether to adjust suspected outlier signals using median of neighboring loci"),
  make_option(c("-w", "--trailingN"), type = "double", default = 5,
              help = "Integer indicating the window of point around a locus to determine if it is an outlier. Defaults to 5")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# check that the depth of coverage file is provided
if (!is.na(opt$pon) & !is.na(opt$tumor_cts)) {
  pon <- opt$pon
  tumor_cts <- opt$tumor_cts
}else {
  stop("User must provide panel of normals, and the tumor counts file to be normalized. See script usage (--help)")
}

pon <- as.data.frame(readRDS(pon))


##################################################################
##      READ IN RELEVANT INPUTS/ GENERATE DOWNSTREAM INPUT      ##
##################################################################
#read in and extract tumor counts information formated in a similar manner as the pon file
sample.name <- sub('_counts.txt','',basename(tumor_cts))

extractCount <- function(df) {
  df.ext <- df %>%
    filter(CONTIG %in% paste0('chr', c(1:22, 'X', 'Y'))) %>%
    unite(col=position, c('START', 'END'), sep='-') %>%
    unite(col=loci, c('CONTIG', 'position'), sep=':') %>%
    column_to_rownames('loci')
  return(df.ext)
}

doc.t <- read_delim(tumor_cts, comment='@') %>%
  extractCount
colnames(doc.t) <- sample.name
#doc.t needs to be subject to the same scaling paradigms and order as the normal doc.
doc.t <- doc.t[rownames(pon), , drop=F]
  
  
if (opt$fix_outliers){
  
  #annotate probes
  hg38 <- hg38 %>%
    mutate(chrom=paste0('chr',chrom)) %>%
    mutate(chrom=case_when(chrom=='chr23' ~ 'chrX', chrom=='chr24' ~ 'chrY', TRUE ~ as.character(chrom)))
  
  probes <- rownames(pon) %>%
    data.table %>%
    setNames('locus') %>%
    separate(locus,into=c('Chr','loci'),sep = ':')%>%
    separate(loci,into=c('Start','End'),sep = '-') %>%
    mutate(across(all_of(c('Start','End')),as.numeric))%>%
    left_join(hg38, by=c('Chr'='chrom')) %>%
    mutate(arm=case_when(End < centromerStart ~ 'p',
                         Start > centromerEnd ~ 'q')) %>%
    mutate(locus=paste(Chr,paste(Start,End,sep = '-'),sep = ':'))
  
  trailingN <- opt$trailingN
  
  remove.outliers <- function(list) {
    list.new <- list %>%
      as.data.frame() %>%
      setNames('signal') %>%
      bind_cols(probes %>% select(locus, Chr, arm)) %>%
      group_by(Chr, arm)
    
    for (i in 1:trailingN) {
      column.name.last <- paste0('last', i)
      column.name.next <- paste0('next', i)
      list.new <- list.new %>%
        mutate(!!column.name.last := lag(signal, n=i, default=dplyr::first(signal))) %>%
        mutate(!!column.name.next := lead(signal, n=i, default=dplyr::last(signal))) %>%
        select(locus, Chr, arm, signal, starts_with('last'), starts_with('next'))
    }
    
    list.outlier.detected <- list.new %>%
      ungroup() %>%
      mutate(last_N_median=select(., starts_with('last')) %>% as.matrix %>% matrixStats::rowMedians()) %>%
      mutate(next_N_median=select(., starts_with('next')) %>% as.matrix %>% matrixStats::rowMedians()) %>%
      mutate(outlier=case_when(signal!=0 &
                                 abs(signal - last_N_median) > 0.3 &
                                 abs(signal - next_N_median) > 0.3 &
                                 abs(log(signal) - log(last_N_median)) > 4 &
                                 abs(log(signal) - log(next_N_median)) > 4 ~ TRUE,
                               TRUE ~ FALSE))
    
    list.replaced <- list.outlier.detected %>%
      mutate(signal.new=case_when(outlier==TRUE ~ select(., c('last1', 'signal', 'next1')) %>% as.matrix %>% matrixStats::rowMedians(), TRUE ~ signal)) %>%
      mutate(index=1:n())
    
    list.outlier.detected.return <- list.outlier.detected %>%
      column_to_rownames('locus') %>%
      select(outlier)
    
    list.replaced.return <- list.replaced %>%
      column_to_rownames('locus') %>%
      select(signal.new)
    
    return(list(new.signal=list.replaced.return, outlier.detected=list.outlier.detected.return))
  }
  
  doc.outlier.removed.list <- lapply(doc.t, remove.outliers)
  
  doc.new.signal.list <- lapply(doc.outlier.removed.list, '[[', 'new.signal')
  doc.t <- doc.new.signal.list %>%
    bind_cols() %>%
    setNames(colnames(doc.t))
}

samp.mean <- mean(doc.t[!grepl('chrX|chrY', rownames(doc.t)),])%>%
  as.data.frame() %>%
  setNames('mean')

## Options (Based on "early_gistic_prep.m") #ADOPTED FROM KEI'S CODE; NOT SURE HOW THIS THRESHOLD WAS DECIDED
floor.val.frac <- 0.001 # Threshold value specified as double within [0,1].
#       CN values will be floored at floor.val.frac x (data mean)

floor.vals <- samp.mean %>%
  mutate(floor.vals=mean * floor.val.frac)
doc.0rep <- doc.t
for (i in 1:length(doc.t)) {
  print(i)
  min.val <- floor.vals[i, 'floor.vals']
  vals <- doc.t[,i] %>% unlist()
  doc.0rep[vals < min.val, i] <- min.val
}

## 4.2. Log2 transformation
doc.log2 <- log2(doc.0rep)
scale_method <- opt$`scale-mthd`

if(is.null(scale_method)){
  stop("Error: Scaling method must be provided. Can be either mean or median")
}else{
  scale.method <- scale_method
}

if (scale.method=='mean') {
  print("Scaling method = 'mean'")
  mean_val <- mean(doc.log2[!grepl('chrX|chrY', rownames(doc.log2)),])
  doc.scaled <- t(t(doc.log2) - mean_val)
} else{
  print("Scaling method = 'median'")
  median_val <- median(doc.log2[!grepl('chrX|chrY', rownames(doc.log2)),])
  doc.scaled <- t(t(doc.log2) - median_val)
}



#################################################################
##                    TANGENT NORMALIZATION                    ##
#################################################################

#clinical data for the reference doc file and the current cohort
ref_clinical <- fread(opt$ref_cli)
if(opt$clin=="TRUE"){ #no clinical data provided
  clinical_data <- ref_clinical[,c('sample_id','gender')]
  clinical_data <- unique(clinical_data,by='sample_id')
}else{
  cohort_clinical <- fread(opt$clin)
  clinical_data <- rbind(ref_clinical[,c('sample_id','gender')],cohort_clinical[,c('sample_id','gender')])
  clinical_data <- unique(clinical_data,by='sample_id')
}

samples <- colnames(pon)

female.samples <- clinical_data %>%
  filter(sample_id %in% samples) %>%
  filter(gender=='female') %>%
  pull(sample_id) %>%
  unique()
male.samples <- clinical_data %>%
  filter(sample_id %in% samples) %>%
  filter(gender=='male') %>%
  pull(sample_id) %>%
  unique()

female.x.mean <- pon[grepl('chrX', rownames(pon)),] %>%
  select(all_of(female.samples)) %>%
  as.matrix() %>%
  mean()
male.x.mean <- pon[grepl('chrX', rownames(pon)),] %>%
  select(all_of(male.samples)) %>%
  as.matrix() %>%
  mean()
female.x.sd <- pon[grepl('chrX', rownames(pon)),] %>%
  select(all_of(female.samples)) %>%
  as.matrix() %>%
  sd()
male.x.sd <- pon[grepl('chrX', rownames(pon)),] %>%
  select(all_of(male.samples)) %>%
  as.matrix() %>%
  sd()

pon.xy.ztransformed <- pon[grepl('chrX|chrY', rownames(pon)), ] %>%
  rownames_to_column('locus') %>%
  separate(col=locus, into=c('chr', 'pos'), sep=':') %>%
  pivot_longer(names_to='sample_id', values_to='signal', cols=-c('chr', 'pos')) %>%
  mutate(signal=case_when(sample_id %in% male.samples & chr=='chrX' ~ ((signal - male.x.mean)/male.x.sd) * female.x.sd + female.x.mean,
                          TRUE ~ signal)) %>%
  pivot_wider(names_from='sample_id', values_from='signal') %>%
  unite(col=locus, c('chr', 'pos'), sep=':') %>%
  column_to_rownames('locus')

pon.ztransformed <- pon[!grepl('chrX|chrY', rownames(pon)), ] %>%
  bind_rows(pon.xy.ztransformed)

#becasue samples is derived from the colnames of the panel of normals
male.normals <- clinical_data %>%
  filter(sample_id %in% samples & gender=='male') %>%
  pull(sample_id) %>%
  unique()



sample_sex <- opt$sample_sex
doc.t <- doc.scaled
if (sample_sex=='male'){
  N <- pon[, male.normals] %>%
    as.matrix() #N refers to male normals
  
  X <- doc.t %>% 
    as.matrix()
  N.means <- apply(N, 1, mean)
  N0 <- N - N.means
  T0 <- X - N.means
  
  Npi <- pracma::pinv(N0)
  weights <- Npi %*% T0
  proj <- N0 %*% weights
  Tn <- T0 - proj #Tn: TUMOR NORMALIZED
  
  ## Re-scaling after Tangent (by median)
  print('Re-scaling...')
  Tn.medians <- Tn[!grepl('chrX|chrY', rownames(Tn)),] %>%
    median
  
  Tn.rescaled <- t(t(Tn)- Tn.medians)
  
  ## Adjust chrY so that it is relative to CN=2
  Tn.rescaled.auto <- Tn.rescaled[!grepl('chrX|chrY', rownames(Tn.rescaled)),]
  Tn.rescaled.xy <- Tn.rescaled[grepl('chrX|chrY', rownames(Tn.rescaled)),]
  Tn.rescaled.xy.adjusted <- Tn.rescaled.xy - 1
  Tn.sexmatched.male <- Tn.rescaled.auto %>%
    as.data.frame() %>%
    bind_rows(Tn.rescaled.xy.adjusted %>% as.data.frame()) %>%
    as.matrix()
  
}else{ #generates an empty data table with the same size as the expected chrY portion of the result
  Tn.sexmatched.male <- setNames(data.frame(matrix(ncol = 1, nrow = nrow(subset(pon,grepl('chrY',rownames(pon)))))), colnames(doc.t))
  rownames(Tn.sexmatched.male) <- rownames(subset(pon,grepl('chrY',rownames(pon))))
}

N.autox <- pon.ztransformed[!grepl('chrY', rownames(pon.ztransformed)), ] %>%
  as.matrix() #here N.autox is the subset of the z.transformed data that has just autosomes and X signal from the normals
T.autox <- doc.t[!grepl('chrY', rownames(doc.t)),] %>%
  as.matrix() #T.autox is the same as above but for the tumors


## SVD
N.autox.svd <- svd(N.autox)
colnames(N.autox.svd$u) <- colnames(N.autox)
rownames(N.autox.svd$u) <- rownames(N.autox)
rownames(N.autox.svd$v) <- colnames(N.autox)


dim.i <- min(ncol(pon),opt$n_dim) #using n_dim dimensions; can't use more dimensions that there are dimensions in the normal subspace

N.lower.i <- N.autox.svd$u[, 1:dim.i] %*% diag(N.autox.svd$d[1:dim.i]) %*% t(N.autox.svd$v[,1:dim.i])
## Get the origin in the normal subspace
N.means <- apply(N.lower.i, 1, mean)
N.lower.i0 <- N.lower.i - N.means
T.autox0 <- T.autox - N.means

Npi.autox <- pracma::pinv(N.lower.i0)
weights.autox <- Npi.autox %*% T.autox0
proj.autox <- N.lower.i0 %*% weights.autox
Tn.autox <- T.autox0 - proj.autox #Tn: TUMOR NORMALIZED

## Re-scaling after Tangent (by median)
print('Re-scaling...')
Tn.auto.medians <- Tn.autox[!grepl('chrX', rownames(Tn.autox)),] %>%
  median

Tn.normalized <- t(t(Tn.autox)- Tn.auto.medians)
Tn.i <- Tn.normalized %>%
  as.data.frame()

#get the final normalized file by binding the rows of the autox normalized output with that of the sex-matched male
Tn.sexmatched.male <- Tn.sexmatched.male %>%
  as.data.frame()
colnames(Tn.sexmatched.male) <- colnames(Tn.i)
Tn.combined <- Tn.i %>%
  bind_rows(subset(Tn.sexmatched.male, grepl('chrY', rownames(Tn.sexmatched.male)))) %>%
  replace(is.na(.), -10)

colnames(Tn.combined) <- colnames(Tn.i)


tangent_path <- 'posttangent_log2.rds'
saveRDS(Tn.combined, tangent_path)
