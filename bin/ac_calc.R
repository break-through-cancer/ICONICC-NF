#!/usr/bin/env Rscript
library(data.table)
library(reshape2)
library(scales)
library(latex2exp)
library(tools)
library(ggplot2)
library(rhdf5)


args <- commandArgs(trailingOnly = T)

case <- args[1] #the participant id in case we want to use this to recognize the outputs
tumor_counts <- args[2]


#if the file is provided as a hdf5 or as a txt, treat 'em differently
if (file_ext(tumor_counts)=='txt'){
  read_depth <- fread(tumor_counts, skip = "CONTIG")
} else{
  if(file_ext(tumor_counts)=='hdf5'){ #assuming the naming conventions are all the same for the hdf5 counts files
    counts <- h5read(file = tumor_counts,name = "counts/values")
    intervals <- h5read(file = tumor_counts,name = "intervals/transposed_index_start_end")
    contigs <- h5read(file = tumor_counts,name = "intervals/indexed_contig_names")
    intervals <- data.table(intervals)
    intervals[,V1:=contigs[V1+1]]
    read_depth <- cbind(intervals,counts)
    colnames(read_depth) <- c('CONTIG','START','END','COUNT')
  }
  else{
    stop("Please provide the read counts file as either a .txt or .hdf5 file")
  }
}

autosomes <- sprintf('chr%s',c(1:22))
canonical_contigs <- c(autosomes, sprintf('chr%s',c('X', 'Y')))

read_depth <- read_depth[CONTIG %in% canonical_contigs]

#Checking if the bins are of equal sizes
test <- read_depth[,.(LENGTH=END-START+1)]
bins <- unique(test$LENGTH)[which.max(tabulate(match(test$LENGTH,unique(test$LENGTH))))]#the mode of bin lengths
#we assume consistent bin intervals were used
bin_size <- test[LENGTH==bins,.N]/nrow(test)
paste0("This ratio should be close to 1 if the bins in the tumor counts file are mostly of equal sizes: ",bin_size)

scale_count <- function(x){
  return ((sum(x)*(desired_binsize*1000))/(length(x)*bins))
}
re_locate <- function(x){
  return(cbind('START'=min(x),'END'=max(x)))
}
aggregate_into_bins <- function(chr,data) { #do aggregation by chr and lapply
  #first find gaps
  #data <- as.data.table(data)
  rd_sub <- data[CONTIG==chr,] #subset by the chr of interest
  gap_length <- floor((matrix(rd_sub$START[2:nrow(rd_sub)])-matrix(rd_sub$END[1:(nrow(rd_sub)-1)]))/bins)
  rd_sub[,"gap_length":=c(0,gap_length)]#if the length is 1 at a locus, it tell us that right before
  #that locus, there is a 'bins' sized gap so the current locus is discontinuous from the previous one
  # Calculate the bin number for each row
  rd_sub$bin <- floor(rd_sub$START / (desired_binsize*1000))
  for (i in which(rd_sub$gap_length!=0)){#push things into a new bin if there's a preceeding gap or not
    rd_sub$bin[i:nrow(rd_sub)] <- rd_sub$bin[i:nrow(rd_sub)] + rd_sub$gap_length[i]
  }
  #rd_sub[, "N" := (.N), bin]
  # Group by chromosome and bin, then sum the reads
  aggregate_count <- aggregate(COUNT ~ bin, rd_sub, FUN = scale_count)
  aggregate_locus <- aggregate(cbind(START,END) ~ bin, rd_sub, FUN = re_locate)
  aggregate_locus <- cbind('CONTIG'=chr,'bin'=aggregate_locus$bin,
                           'START'=as.numeric(aggregate_locus$START[,1]),'END'=as.numeric(aggregate_locus$END[,2]))
  aggregated_data <- merge(aggregate_count,aggregate_locus,by='bin')
  aggregated_data <- data.table(aggregated_data)
  aggregated_data[,bin:=NULL]
  return(aggregated_data)
}

aggregate_bins <- args[3] #should be true or false indicating if the user wishes to aggregate bins or not
if(aggregate_bins){
  desired_binsize <- args[4] #the desired bin_size must be provided right after
  desired_binsize <- as.integer(desired_binsize) #convert the entry to an integer
  not_int <- is.na(desired_binsize) #check if the result is not an integer
  if(not_int){
    stop("Please provide the desired bin size in kb to aggregate your counts data. E.g., input 50 for a desired bin size of 50kb")
  }else{
    cat('Warning: Please ensure that the desired bin size in provided in kb multiples (e.g., input 50 for 50kb bins NOT 50000) otherwise results will be inaccurate.\n
        We recommend aggregating read depth data over 25kb bins if your intervals are less than 25kb.\n')
    #how many contiguous bins at a time need to be taken together to make the aggregate
    bin_no <- (desired_binsize*1000)/bins #number of bins that should be grouped
    aggregate_readdepth <- as.data.table(rbindlist(lapply(canonical_contigs,aggregate_into_bins,read_depth)))
    #check if the resulting aggregate bin sizes are about the same size
    test2 <- aggregate_readdepth[,.(LENGTH=as.numeric(END)-(as.numeric(START)+1))]
    bins2 <- unique(test2$LENGTH)[which.max(tabulate(match(test2$LENGTH,unique(test2$LENGTH))))]#the mode of bin lengths
    #we assume consistent bin intervals were used
    bin_size2 <- test2[LENGTH==bins2,.N]/nrow(test2)
    cat(paste0("This ratio should be close to 1 if the bins in the resulting tumor counts aggregate file are mostly of equal sizes: ",bin_size2),'\n')
    
    normal_sample_id <- args[5] #changing this to access the control_sample_id/normal sample id
    tumor_sample_id <- args[6] #this also accesses the case_sample_id
    all_input_files <- args[ - (1:6)] #takes from the 7th argument onward, 
  }
  
}else{
  cat('Warning: We recommend aggregating read depth data over 25kb bins if your intervals are less than 25kb.\n')
  aggregate_readdepth <- read_depth
  normal_sample_id <- args[4] #changing this to access the control_sample_id/normal sample id
  tumor_sample_id <- args[5] #this also accesses the case_sample_id
  all_input_files <- args[ - (1:5)] #takes from the 6th argument onward
  
}


#Normalizing the total read depth
mrd <- aggregate_readdepth[CONTIG%in%autosomes,.(med=median(COUNT, na.rm = T))][[1]]
aggregate_readdepth[,NORMALIZED_COUNT:=COUNT/mrd]

#Allelic fractions
read_data <- function(chr, directory, prefix, suffix, columns, sep="auto"){
  data_file <- paste0(directory,'/', prefix, sprintf(suffix, chr))
  data_table <- fread(data_file, col.names = columns, sep = sep)
  return(data_table)
}



normal_prefix <- normal_sample_id #these new lines are no longer dependent on a different filename to tell what the tumor and normal prefixes should be.
tumor_prefix <- tumor_sample_id

#the rest of the files are entered into a collective variable argument
#all_input_files <- args[ - (1:4)] #takes from the 5th argument onward
directories <- unique(lapply(all_input_files, FUN = dirname)) #extracts the unique directory names from all the input file


phased_GT_directory <- directories[1]

allelic_depths_directory <- directories[2]

hetinfo_files_directory <- directories[3]
het_prefix <- paste0(case,'_')

cat('Reading phased_GT files...\n')
l <- lapply(canonical_contigs[1:23], read_data, phased_GT_directory, normal_prefix, ".phasedGT.%s.1.txt", c("GT1", "GT2"), "|")
phased_GT <- rbindlist(l, idcol = "CHR")
cat('done \n')

cat('Reading AD files... \n')
l <- lapply(canonical_contigs[1:23], read_data, allelic_depths_directory, tumor_prefix, ".AD.%s.txt", c("REF", "ALT"))#1:23 because we don't have info on '24'[Y]
allelic_depth <- rbindlist(l, idcol = "CHR")
cat('done \n')

cat('Reading hetinfo files... \n')
l <- lapply(canonical_contigs[1:23], read_data, hetinfo_files_directory, het_prefix, "%s_HetInfo.txt", c("CHROM","POS", "REF", "ALT", "DP"))
Het_info <- rbindlist(l)
cat('done \n')

variants <- cbind(Het_info[,.(CHROM,POS)], allelic_depth[,.(REF, ALT)], phased_GT[,.(GT1, GT2)])

all_het <- variants[GT1!=GT2, .N]/nrow(variants)
paste0("This ratio should be close to 1 if all the variants are heterozygous: ", all_het)

variants[, c("AD1", "AD2") := list(REF*(1-GT1)+ALT*GT1, REF*(1-GT2)+ALT*GT2)]

#calculating allelic read depth
data_list <- list()
for (chr in canonical_contigs[1:23]){
  rdc <- aggregate_readdepth[CONTIG==chr]
  vc <- variants[CHROM==chr]
  bin_breaks <- sort(c(rdc$START,rdc$END))
  #variants_grouped <- cut(vc$POS, bin_breaks, labels=F)
  vc[, "bin" := cut(POS, bin_breaks, labels=F)]
  vc[, c("AF1", "AF2", "snp_count") := list(mean(AD1/(AD1+AD2)), mean(AD2/(AD1+AD2)), .N), bin]
  af <- na.omit(unique(vc[,.(bin, AF1, AF2, snp_count)], by="bin"))
  rdc[, "bin" := 2*seq(.N)-1]
  md <- af[rdc, on=.(bin)][,-c("bin")] #merge it in such a way that bins that have read counts but no variants are still retained with NA for this stage
  #md[, c("C1", "C2") := list(NORMALIZED_COUNT*AF1, NORMALIZED_COUNT*AF2)]
  data_list <- append(data_list, list(md))
}
allelic_counts <- rbindlist(data_list)

#generate a random assortment of 0's and 1's for how ever many NA's there are in the 'N' column
#then fill those into the AF1 column where there are NAs and make the AF2 column abs(1-AF1)
set.seed(1234)
rand_assign <- rbinom(nrow(allelic_counts[is.na(snp_count),]),1,0.5)
na_indices <- which(is.na(allelic_counts$snp_count))
for (i in 1:length(na_indices)){
  allelic_counts$AF1[na_indices[i]] <- rand_assign[i]
  allelic_counts$AF2[na_indices[i]] <- abs(1-allelic_counts$AF1[na_indices[i]])
}

allelic_counts[,c("C1", "C2") := list(NORMALIZED_COUNT*AF1, NORMALIZED_COUNT*AF2)]
allelic_counts[,snp_count:=ifelse(is.na(snp_count),0,snp_count)]

read_depth_out <- paste0(case,"_read_depth.txt")
allelic_counts_out <- paste0(case, "_allelic_counts.txt")
plots_file <- paste0(case,'_plots.pdf')


cat(plots_file,'\n')
#add plotting code
pdf(file = plots_file)
gg <- ggplot(aes(C1, C2), data=allelic_counts)+
  theme_bw()+
  geom_point(alpha=0.05, size=0.1)+
  scale_x_continuous(name=TeX("$CN_1$"), limits = c(0,2))+
  scale_y_continuous(name=TeX("$CN_2$"), limits = c(0,2))+ggtitle("All Allelic Copy Numbers")
gg
g <- ggplot(aes(C1, C2), data=allelic_counts[snp_count>5])+
  theme_bw()+
  geom_point(alpha=0.05, size=0.1)+
  scale_x_continuous(name=TeX("$CN_1$"), limits = c(0,2))+
  scale_y_continuous(name=TeX("$CN_2$"), limits = c(0,2)) +ggtitle("Allelic Copy Numbers with >5 variants")
g
dev.off()


write.table(aggregate_readdepth, read_depth_out, row.names = F) #now outputs the aggregated read_depth info
write.table(allelic_counts, allelic_counts_out, row.names = F)