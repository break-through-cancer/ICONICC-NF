#!/usr/bin/env Rscript
#################################################################
##                      READ IN LIBRARIES                      ##
#################################################################
#tangent wdl
library(useful)
library(tidyverse)
library(pracma)
library(optparse)
library(data.table)
library(reshape2)
library(scales)
library(latex2exp)
library(tools)
library(ggplot2)
library(rhdf5)
library("gridExtra")
library(DNAcopy)
library(rCGH)
library(Hmisc)
library(gUtils)
library(GenomicRanges)
library(parallel)

##################################################################
##                DEFINE INPUT OPTIONS AND FLAGS                ##
##################################################################
option_list <- list(
  make_option(c("-f", "--al_cts"), type = "character", default = NA,
              help = "Path to allelic counts file"),
  make_option(c("-t", "--tangent_cts"), type = "character", default = NA,
              help = "Path to tangent normalized read counts"),
  make_option(c("-p", "--participant_id"), type = "character", default = NA,
              help = "Participant identifier for the sample"),
  make_option(c("-r", "--smth_reg"), type = "integer", default = 10,
              help = "Smoothing region for smooth.CNA method"),
  make_option(c("-o", "--outlier_SD"), type = "integer", default = 4,
              help = "Number of SDs from nearest point to call outlier during initial segment generation"),
  make_option(c("-s", "--smooth_SD"), type = "integer", default = 2,
              help = "Number of SDs from median of smoothing region to position smoothed point"),
  make_option(c("-a", "--alpha"), type = "double", default = 0.01,
              help = "significance level for test to accept new segments"),
  make_option(c("-w", "--min_width"), type = "integer", default = 2,
              help = "Minimum number of markers to include in a segment"),
  make_option(c("-x", "--undo_split"), type = "character", default = "none",
              help = "method for undoing change points. Options include c('none','prune','sdundo')"),
  make_option(c("-n", "--undo_prune"), type = "double", default = 0.05,
              help = "Proportional increase in sum of squares when undo_split='prune'"),
  make_option(c("-d", "--undo_SD"), type = "integer", default = 2,
              help = "SD between means to retain split when undo_split='sdundo'"),
  make_option(c("-g", "--sample_sex"), type = "character", default = NA,
              help = "Sex of sample. Takes values 'male' or 'female'"),
  make_option(c("-z", "--run_recovery"), type = "logical", default = TRUE,
              help = "Boolean indicating whether to recover outliers from Tangent Segments (defaults to true)"),
  make_option(c("-c", "--recovery_thresh"), type = "double", default = 20,
              help = "Amplitude/Standard Error threshold to annotate and recover outliers from Tangent Segments"),
  make_option(c("-i", "--use_hetinfo"), type = "logical", default = FALSE,
              help = "Boolean indicating whether to check if an outlier region overlaps bins without heterozygous snps (defaults to false)"),
  make_option(c("-m", "--miss_thresh"), type = "double", default = 1.0,
              help = "Fraction of missing heterozygous snps in a segment to confirm outlier status. This value must be set if use_hetinfo is TRUE (defaults to 1)")

)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


##################################################################
##                          WDL FORMAT                          ##
##################################################################
##INPUTS
#1. Allelic counts
#2. Tangent Normalized Count
#3. Participant/sample id

allelic_counts <- fread(opt$al_cts)
tangent_norm <- readRDS(opt$tangent_cts)
colnames(tangent_norm) <- "log2_tangent"
case_id <- opt$participant_id #participant_id
sample_sex <- opt$sample_sex

#optional arguments for segmentation proper
smooth_reg=opt$smth_reg
outlier_SD_sc=opt$outlier_SD
smooth_SD_sc=opt$smooth_SD
alp=opt$alpha
min_wid=opt$min_width
undo_spl=opt$undo_split
undo_pr=opt$undo_prune
undo_SD=opt$undo_SD

#arguments for segment recovery
run_recovery <- opt$run_recovery
outlier_thresh <- opt$recovery_thresh
use_hetinfo <- opt$use_hetinfo
miss_thresh <- opt$miss_thresh



##################################################################
##                        FUNCTION BLOCK                        ##
##################################################################

#Function to generate weighted segmentation profiles
generate_segments_wghts <- function(mat,df,wghts,smooth.reg=10,outlier.SD.sc=4,smooth.SD.sc=2,
                                    alp=0.01,min.wid=2,undo.spl="none",undo.pr=0.05,undo.SD=3){ #generate the matrix outside,feed source df, case
  CNA.mat <- CNA(mat, as.factor(df$CONTIG),
                 as.numeric(df$START),data.type = 'logratio', sampleid = case_id)
  smoothed.CNA <- smooth.CNA(CNA.mat,smooth.region =smooth.reg,outlier.SD.scale=outlier.SD.sc,smooth.SD.scale=smooth.SD.sc)
  segment.CNA <- segment(smoothed.CNA,weights = wghts,
                         alpha=alp,min.width=min.wid,undo.splits=undo.spl,undo.prune=undo.pr,undo.SD=undo.SD)
  segments <- as.data.table(segment.CNA[["output"]])
  segments <- segments[order(factor(chrom,levels=canonical_contigs),as.numeric(loc.start))]
  segments[,SegmentID:=1:nrow(segments)]
  segments$loc.start <- as.numeric(segments$loc.start);segments$loc.end <- as.numeric(segments$loc.end)
  return(segments)
}

#Function to assign loci to defined segments from segmentation
assign_group <- function(itter,dt,segments_CNA,seg.source){ #a function designed to assign a 'unique' group to a signal in the tumor space based on the ar_segments so we can aggregate on these groups
  #cat(itter,'\n')
  row_c <- dt[itter,]
  chr <- row_c$CONTIG
  segments_info <- subset(segments_CNA,chrom==chr)
  seg_group <- segments_info[which((as.numeric(row_c$START)>=loc.start)&(as.numeric(row_c$END)<=loc.end)),]$SegmentID
  if (identical(integer(0),seg_group)){
    seg_group <- NA
  }
  return(row_c[,paste0(seg.source,'_segment'):=seg_group])
}

#Generating genome-wide plots
canonical_contigs <- c(paste0("chr",1:22),"chrX","chrY")

hg38 <- rCGH::hg38 %>%
  mutate(chrom=case_when(chrom==23 ~ 'X', chrom==24 ~ 'Y', TRUE ~ as.character(chrom)))%>%
  mutate(chrom=paste0('chr', chrom)) %>%
  as.data.table
hg38[,frac_genome:=as.numeric(length)/sum(as.numeric(hg38$length))]
max_perc <- max(hg38$frac_genome);min_perc <- min(hg38$frac_genome);range_perc<-max_perc-min_perc
new_max_perc <- max_perc-(0.25*range_perc);new_min_perc<-min_perc+(0.25*range_perc);new_range_perc<-new_max_perc-new_min_perc
hg38[,new_frac_genome:=((frac_genome*new_range_perc)-(min_perc*new_max_perc)+(max_perc*new_min_perc))/range_perc,by=chrom]
hg38[,new_length:=new_frac_genome*sum(as.numeric(hg38$length))] #new length of each chromosome
hg38[,new_centromereStart:=(new_length*centromerStart)/length] #new length of each chromosome
hg38[,new_centromereEnd:=(new_length*centromerEnd)/length] #new length of each chromosome
hg38[,cum_end := cumsum(new_length)]

#data for genome plot
strip_data <- as.data.table(cbind(as.vector(c(1,as.numeric(hg38$cum_end[1:22])+1)),as.numeric(hg38$cum_end)))
strip_data[,contigs := as.vector(canonical_contigs)]
colnames(strip_data) <- c('xmin','xmax','CONTIG')
strip_data[,c('ymin','ymax'):=list(-2,2)] #changed from 0 to 3 --> -2 to 2 since we're now working in log 2 only. Some values may still fall outside this range
strip_data[,fill:= rep(c('a','b'),length.out=nrow(strip_data))]
strip_data$xmin<-as.numeric(strip_data$xmin);strip_data$xmax<-as.numeric(strip_data$xmax)
strip_data$ymin<-as.numeric(strip_data$ymin);strip_data$ymax<-as.numeric(strip_data$ymax)

hg38[,cum_START:=strip_data$xmin]; hg38[,cum_END:=strip_data$xmax]
#lose chrY in strip data
strip_data <- strip_data[-24,]

get_segdf <- function(itter,seg_file){
  row_l = seg_file[itter,]
  new_tot_length = as.numeric(as.character(hg38[chrom==row_l$chrom,]$new_length))
  old_tot_length = as.numeric(as.character(hg38[chrom==row_l$chrom,]$length))
  new.loc.start <- as.numeric(as.character(row_l$loc.start))*(new_tot_length/old_tot_length)
  x <- as.numeric(hg38[chrom==row_l$chrom,]$cum_START)+new.loc.start
  new.loc.end <- as.numeric(as.character(row_l$loc.end))*(new_tot_length/old_tot_length)
  x.end <- as.numeric(hg38[chrom==row_l$chrom,]$cum_START)+new.loc.end
  #print(new.loc.end)
  y <- row_l$seg.mean
  y.end <- row_l$seg.mean
  #print(row_l$new.loc.end)
  return(data.table(cbind(x,y,x.end,y.end)))
}

#Function to generate segments without weighting (relevant for tangent)
generate_segments <- function(mat,df,smooth.reg=10,outlier.SD.sc=4,smooth.SD.sc=2,
                              alp=0.01,min.wid=2,undo.spl="none",undo.pr=0.05,undo.SD=3){ #generate the matrix outside,feed source df, case
  CNA.mat <- CNA(mat, as.factor(df$CONTIG),
                 as.numeric(df$START),data.type = 'logratio', sampleid = case_id)
  smoothed.CNA <- smooth.CNA(CNA.mat,smooth.region =smooth.reg,outlier.SD.scale=outlier.SD.sc,smooth.SD.scale=smooth.SD.sc)
  segment.CNA <- segment(smoothed.CNA,
                         alpha=alp,min.width=min.wid,undo.splits=undo.spl,undo.prune=undo.pr,undo.SD=undo.SD)
  segments <- as.data.table(segment.CNA[["output"]])
  segments <- segments[order(factor(chrom,levels=canonical_contigs),as.numeric(loc.start))]
  segments[,SegmentID:=1:nrow(segments)]
  segments$loc.start <- as.numeric(segments$loc.start);segments$loc.end <- as.numeric(segments$loc.end)
  return(segments)
}

#Funtion to annotate outliers among Tangent segments
outlier_segments <- function(itter,segments.ar,full_dt,outlier_thresh){#check the standard error of the cluster of points that contribute to that segment
  #if the ratio between the distance between the mean of the segments and the standard error of the points that constitute the segment is
  #above the recovery threshold value, then call it an outlier
  segment.i <- segments.ar[itter,]
  segment.ar.i <- dt2gr(segment.i)
  ar.tan.overlap.i <- gr2dt(segments.tan.gr %&% segment.ar.i) #tangent segments that overlap the ar segment
  ar.tan.overlap.i[,perc_overlap := ifelse(start>=segment.i$loc.start&end<=segment.i$loc.end,1,((max(end,segment.i$loc.end)-min(start,segment.i$loc.start))-(abs(start-segment.i$loc.start)+abs(end-segment.i$loc.end)))/(end-start)),by='SegmentID']
  ar.tan.overlap.i <- ar.tan.overlap.i[perc_overlap>0.5]

  if (nrow(ar.tan.overlap.i>=1)){ #if at least one overlap exists
    ar.mean <- segments.ar[itter,]$seg.mean
    ar.tan.overlap.i[,seg.sem:=sd(full_dt[tangent_segment==SegmentID,]$log2_tangent,na.rm=T)/sqrt(nrow(full_dt[tangent_segment==SegmentID,])),by='SegmentID'] #using the standard error of the mean rather than standard deviation
    ar.tan.overlap.i[,SDR_ratio:=(abs(seg.mean-ar.mean))/seg.sem]
    ar.tan.overlap.i[,paste0(outlier_thresh,'SDR'):=ifelse(SDR_ratio<=outlier_thresh,'normal','outlier')]
    return(ar.tan.overlap.i)
  }else{
    return(NULL)
  }
}



#################################################################
##                 PROCESSING AND SEGMENTATION                 ##
#################################################################

allelic_counts_df <- allelic_counts %>%
  filter(CONTIG %in% paste0('chr', c(1:22, 'X', 'Y'))) %>%
  unite(col=loci, c('CONTIG', 'START'), sep=':')

tangent_norm$loci <- gsub('-.*','',rownames(tangent_norm))
tangent_norm$CONTIG <- gsub(':.*','',rownames(tangent_norm))
tangent_norm <- tangent_norm %>%
  filter(CONTIG %in% paste0('chr', c(1:22, 'X')))

full_dt <- full_join(allelic_counts_df,tangent_norm,by='loci')
full_dt$CONTIG <- gsub(':.*','',full_dt$loci)
full_dt <- as.data.table(full_dt)
full_dt[,START:=gsub("^.*:", "", loci)]
full_dt[,c('AFMAX','AFMIN'):=list(pmax(AF1,AF2),pmin(AF1,AF2))]
full_dt[,AF_ratio:=AFMIN/AFMAX]
full_dt$cum_start=cumsum(as.numeric(full_dt$START))
full_dt$cum_end=cumsum(as.numeric(full_dt$END))


set.seed(1234)

#Generate Allelic ratio segments
ar.mat <- matrix(full_dt[snp_count>5]$AF_ratio)
segments.ar <- generate_segments_wghts(ar.mat,full_dt[snp_count>5],wghts = as.numeric(full_dt[snp_count>5]$snp_count),
                                       smooth.reg=smooth_reg,outlier.SD.sc=outlier_SD_sc,
                                       smooth.SD.sc=smooth_SD_sc,alp=alp,min.wid=min_wid,undo.spl=undo_spl,
                                       undo.pr=undo_pr,undo.SD=undo_SD)
segments.ar[,loc.end:=loc.end+24999] #first adjust the end locations of the segments
rm(allelic_counts,allelic_counts_df) #cleanup

full_dt <- full_dt[order(factor(CONTIG,levels = canonical_contigs),as.numeric(START))]
full_dt <- rbindlist(lapply(1:nrow(full_dt),assign_group,full_dt,segments.ar,seg.source='AR'))
cts.mean.ar <- aggregate(full_dt$log2_tangent,by=list(full_dt$AR_segment),FUN = function(x){return (mean(x,na.rm=T))})
colnames(cts.mean.ar) <- c('SegmentID','seg.mean')
segments.ar <- merge(segments.ar[,!'seg.mean'],cts.mean.ar,by='SegmentID')

#modify the AR segments to be adjacent
segments.ar <- segments.ar %>%
  group_by(chrom) %>%
  mutate(loc.end=if_else(row_number()==n(),hg38$length[match(chrom,hg38$chrom)],loc.end)) %>%
  ungroup()

#First: adjust ends (forward-looking)
segments.ar.adj <- segments.ar %>%
  group_by(chrom) %>%
  mutate(
    next_start = lead(loc.start),
    next_seg_mean = lead(seg.mean),
    gap = next_start - loc.end,
    midpoint = floor((gap - 1) / 2),
    new_end = case_when(
      !is.na(gap) & gap > 1 & abs(seg.mean) > abs(coalesce(next_seg_mean, 0))~ next_start - 1,
      !is.na(gap) & gap > 1 & abs(seg.mean) == abs(coalesce(next_seg_mean, 0)) ~ loc.end + midpoint,
      TRUE ~ loc.end
    )
  ) %>%
  mutate(loc.end = new_end) %>%
  select(-new_end, -next_start, -next_seg_mean, -gap, -midpoint) %>%
  ungroup()

# Second: adjust starts (backward-looking)
segments.ar.adj <- segments.ar.adj %>%
  group_by(chrom) %>%
  mutate(
    prev_end = lag(loc.end),
    prev_seg_mean = lag(seg.mean),
    gap = loc.start - prev_end,
    new_start = case_when(
      # Adjust only if there's a gap and current segment has stronger signal
      !is.na(gap) & gap > 1 & abs(seg.mean) > abs(coalesce(prev_seg_mean, 0)) ~ prev_end + 1,

      # Equal signal -> meet in the middle
      !is.na(gap) & gap > 1 & abs(seg.mean) == abs(coalesce(prev_seg_mean, 0)) ~ prev_end + 1 + floor((gap - 1) / 2),

      # Otherwise, leave loc.start unchanged
      TRUE ~ loc.start
    )

  ) %>%
  mutate(loc.start = new_start) %>%
  select(-prev_end, -prev_seg_mean, -gap, -new_start) %>%
  ungroup() %>%
  as.data.table()

segments.ar <- segments.ar.adj
ar.plotdf <- rbindlist(lapply(1:nrow(segments.ar),get_segdf,segments.ar))

ar_old <- ggplot(data=full_dt, aes(cum_start, log2_tangent)) +
  theme_classic()+
  ylim(-2,2) +ggtitle(paste0('AR segments\n',nrow(ar.plotdf),' Segments'))+
  geom_rect(data=strip_data, aes(xmin=xmin,xmax=xmax, ymin=ymin, ymax=ymax,fill=fill),
            inherit.aes=F, alpha = 0.1,show.legend = F)+
  geom_segment(data = ar.plotdf,
               mapping = aes(x=x, y=y, xend=x.end, yend=y.end),
               inherit.aes = FALSE,lwd=1,col='red')+
  theme(plot.title = element_text(hjust = 0.5)) +xlab("Start Position") + ylab("Segment Mean")

#generate Tangent segments
tan.mat <- matrix(full_dt$log2_tangent)
segments.tan <- generate_segments(tan.mat,full_dt,smooth.reg=smooth_reg,outlier.SD.sc=outlier_SD_sc,
                                  smooth.SD.sc=smooth_SD_sc,alp=alp,min.wid=min_wid,undo.spl=undo_spl,
                                  undo.pr=undo_pr,undo.SD=undo_SD)

#correct X segmentation if the sample is male
is_male <- ifelse(sample_sex=='male',TRUE,FALSE)
if(is_male){
  segments.tan.stringent <- generate_segments(tan.mat,full_dt,smooth.reg=smooth_reg,outlier.SD.sc=outlier_SD_sc,
                                              smooth.SD.sc=smooth_SD_sc,alp=1e-3,min.wid=min_wid,undo.spl=undo_spl,
                                              undo.pr=undo_pr,undo.SD=undo_SD)

  segments.tan <- rbind(segments.tan[!grepl('chrX',chrom),],segments.tan.stringent[chrom=='chrX'])
}
segments.tan[,loc.end:=loc.end+24999] #adjust the end locations of the segments
segments.tan[,SegmentID:=1:nrow(segments.tan)]

rm(ar.mat,tan.mat,tangent_norm) #cleanup

if (!run_recovery){ #end all analyses here
  #save the corrected version of segments.ar as the final segmentation output
  segfile <- paste0(case_id,".seg.txt")
  if(is_male){
    seg_file <- rbind(segments.ar[!grepl('chrX',chrom),],segments.tan.stringent[chrom=='chrX',])
    seg_file[,SegmentID:=1:nrow(seg_file)]

    full_dt <- full_dt[order(factor(CONTIG,levels = canonical_contigs),as.numeric(START))]
    full_dt <- rbindlist(lapply(1:nrow(full_dt),assign_group,full_dt,seg_file,seg.source='AR'))
    cts.mean.ar <- aggregate(full_dt$log2_tangent,by=list(full_dt$AR_segment),FUN = function(x){return (mean(x,na.rm=T))})
    colnames(cts.mean.ar) <- c('SegmentID','seg.mean')
    seg_file <- merge(seg_file[,!'seg.mean'],cts.mean.ar,by='SegmentID')

    ar.plotdf <- rbindlist(lapply(1:nrow(seg_file),get_segdf,seg_file))

    ar_old <- ggplot(data=full_dt, aes(cum_start, log2_tangent)) +
      theme_classic()+
      ylim(-2,2) +ggtitle(paste0('AR segments\n',nrow(ar.plotdf),' Segments'))+
      geom_rect(data=strip_data, aes(xmin=xmin,xmax=xmax, ymin=ymin, ymax=ymax,fill=fill),
                inherit.aes=F, alpha = 0.1,show.legend = F)+
      geom_segment(data = ar.plotdf,
                   mapping = aes(x=x, y=y, xend=x.end, yend=y.end),
                   inherit.aes = FALSE,lwd=1,col='red')+
      theme(plot.title = element_text(hjust = 0.5)) +xlab("Start Position") + ylab("Segment Mean")

    seg_file <- seg_file[,SegmentID:=NULL]
  }else{
    seg_file <- segments.ar[,SegmentID:=NULL]
  }
  seg_file <- seg_file[order(factor(chrom,levels = canonical_contigs),as.numeric(loc.start))]
  write.table(seg_file, segfile, row.names = F,quote = F,col.names = T,sep = '\t')

  requiredPlots <- paste0(case_id,'_segmentPlots.pdf')
  pdf(file=requiredPlots)
  print(ar_old)
  dev.off()

  allelic_counts_file <- paste0(case_id,"_processed_counts.txt")
  accts_file <- full_dt[order(factor(CONTIG,levels = canonical_contigs),as.numeric(START))]
  write.table(accts_file, allelic_counts_file, row.names = F,quote = F,col.names = T,sep = '\t')
}


if(run_recovery){
  segments.tan <- segments.tan %>%
    as.data.frame() %>%
    left_join(hg38 %>% select(chrom, centromerStart, centromerEnd), by='chrom') %>%
    mutate(arm=case_when(loc.end < centromerStart ~ 'p',
                         loc.start > centromerEnd ~ 'q'))
  segments.tan <- as.data.table(segments.tan)
  segments.tan$centromerStart <- as.numeric(segments.tan$centromerStart);segments.tan$centromerEnd <- as.numeric(segments.tan$centromerEnd)

  #check if most of the arm (>30%) overlaps with the centromere, if so, call that a centromeric segment
  segments.tan[,arm_adj := ifelse(is.na(arm),ifelse((((max(centromerEnd,loc.end)-min(centromerStart,loc.start))-(abs(loc.start-centromerStart)+abs(loc.end-centromerEnd)))/(loc.end-loc.start))<=0.3,ifelse(loc.start<centromerStart,'p','q'),arm),arm),by='SegmentID']
  segments.tan[,chr_arm:= ifelse(!is.na(arm_adj),paste0(chrom,':',arm_adj),NA)]

  full_dt <- full_dt[order(factor(CONTIG,levels = canonical_contigs),as.numeric(START))]
  full_dt <- rbindlist(lapply(1:nrow(full_dt),assign_group,full_dt,segments.tan,seg.source='tangent'))

  segments.tan.gr <- dt2gr(segments.tan)

  if(is_male){
    segments.ar.ref <- segments.ar[!grepl('chrX',chrom),] #if the sample is male, don't consider chrX for recovery
  }else{
    segments.ar.ref <- segments.ar
  }

  segments.ar.ref <- na.omit(segments.ar.ref,cols="seg.mean") #lose any entries for which the seg.mean is NA; these were AR segments for which they were no tangent probes in the entire bin
  outlier.segments.annot <- rbindlist(lapply(1:nrow(segments.ar.ref),outlier_segments,segments.ar.ref,full_dt,outlier_thresh))

  if(is_male){
    outlier.segments.annot<-rbind(outlier.segments.annot,gr2dt(segments.tan.gr)[seqnames=='chrX',],fill=T) #if the sample is male, bring back the chrX segments
    outlier.segments.annot[which(grepl('chrX',outlier.segments.annot$seqnames)), paste0(outlier_thresh,"SDR") := 'normal'] #changed to reflect new outlier thresholds
  }

  #add the centromeric segments back and flag them as NA
  outlier.segments.annot <- outlier.segments.annot[order(factor(seqnames,levels=canonical_contigs),start)]
  colnames(outlier.segments.annot)[1:3] <- c('chrom','loc.start','loc.end')

  #get final set of outlier segments
  if(nrow(outlier.segments.annot)!=nrow(segments.tan)){
    cat('The following tangent segments did not overlap any AR segments and were not considered in recovery \n')
    print(setdiff(segments.tan$SegmentID,outlier.segments.annot$SegmentID))
    #if there are segments that aren't considered, add them back and just set them as 'normal"?
    rem_segments <- setdiff(segments.tan$SegmentID,outlier.segments.annot$SegmentID)
    segments.tan.rem <- subset(gr2dt(segments.tan.gr),gr2dt(segments.tan.gr)$SegmentID %in% rem_segments)
    colnames(segments.tan.rem)[1:3] <- c('chrom','loc.start','loc.end')
    outlier.segments.annot <- rbind(outlier.segments.annot,segments.tan.rem,fill=T)
    outlier.segments.annot[which(SegmentID%in% rem_segments), paste0(outlier_thresh,"SDR") := 'normal']
    outlier.segments.annot <- outlier.segments.annot[order(factor(chrom,levels=canonical_contigs),loc.start)]
  }

  #generate tangent plots
  tan_plotdf <- rbindlist(lapply(1:nrow(segments.tan),get_segdf,segments.tan))
  outlier.colname <- paste0(outlier_thresh,'SDR')
  tan_plotdf[,outlier.class:=outlier.segments.annot[,..outlier.colname]]
  tan_plotdf[,chr_arm:= outlier.segments.annot[,chr_arm]]
  tan_plotdf[,outlier.class:=ifelse(is.na(chr_arm),'centromeric',outlier.class)]
  outlier.class <- tan_plotdf$outlier.class

  tangent_plot <- ggplot(data=full_dt, aes(cum_start, log2_tangent)) +
    theme_classic()+
    ylim(-2,2) +
    geom_rect(data=strip_data, aes(xmin=xmin,xmax=xmax, ymin=ymin, ymax=ymax,fill=fill),
              inherit.aes=F, alpha = 0.1,show.legend = F)+
    geom_segment(data = tan_plotdf,
                 mapping = aes(x=x, y=y, xend=x.end, yend=y.end,col=outlier.class),
                 inherit.aes = FALSE,lwd=1)+scale_color_manual(values = c('outlier'='blue','normal'='red','centromeric'='green'))+
    theme(plot.title = element_text(hjust = 0.5)) +xlab("Start Position") + ylab("Segment Mean") +
    ggtitle(paste0('Tangent Segments: ',nrow(tan_plotdf),' Segments\n Outliers at ', outlier_thresh, ' SDR(s) from AR neighborhood mean annotated'))+
    guides(color = guide_legend(title = "Outlier Class")) +theme(legend.position = "none")


  #prepare outliers to recover
  outliers <- outlier.segments.annot[get(paste0(outlier_thresh,'SDR'))=='outlier'] #work with outliers at the set outlier threshold

  segments.ar[,seg_nature:='original']
  segments.ar$loc.start <- as.numeric(segments.ar$loc.start);segments.ar$loc.end <- as.numeric(segments.ar$loc.end)

  outliers$loc.start <- as.numeric(outliers$loc.start);outliers$loc.end <- as.numeric(outliers$loc.end)
  if(nrow(outliers)>0){
    if (use_hetinfo){ #if overlapping outliers with missing locations
      cat('Using the default missingness threshold of 1. \nIf you would like a lower threshold for the fraction of bins within a segment that lack heterozygous snps before that segment can be tagged for recovery, set this value using the --miss_thresh flag \n')
      missed_loci <- full_dt[snp_count<5]
      missed_GR <- makeGRangesFromDataFrame(missed_loci, keep.extra.columns = T,
                                            seqnames.field="CONTIG",
                                            start.field="START",
                                            end.field="END",
                                            ignore.strand = T)
      outliers.mod <- data.table(matrix(nrow = 1,ncol = ncol(outliers)+1));colnames(outliers.mod)=c(colnames(outliers),'is_seg_missing')
      for (i in 1:nrow(outliers)){
        row_l <- outliers[i,]
        row_l_gr <- GRanges(row_l$chrom, IRanges(as.numeric(as.character(row_l$loc.start)), as.numeric(as.character(row_l$loc.end))))
        tangent_segid <- row_l$SegmentID
        perc_missing <- nrow(full_dt[(tangent_segment==tangent_segid) & !is.na(log2_tangent) & snp_count<5,])/nrow(full_dt[tangent_segment==tangent_segid & !is.na(log2_tangent),])
        is_seg_missing <- ifelse(perc_missing==miss_thresh,1,0)
        row_l <- cbind(row_l,is_seg_missing)
        outliers.mod <- rbind(outliers.mod,row_l,use.names=F)
      }
      outliers <- outliers.mod[is_seg_missing==1,]
    }
    #check again if outliers is empty or not
    if (nrow(outliers)>0){
      for (i in 1:nrow(outliers)){
        segment_i <- outliers[i,]
        segment_i[,seg_nature:='recovered']
        ref_GR <- makeGRangesFromDataFrame(segments.ar, keep.extra.columns = T,
                                           seqnames.field="chrom",
                                           start.field="loc.start",
                                           end.field="loc.end",
                                           ignore.strand = T)
        ref_sub <- segments.ar[ chrom == segment_i$chrom]
        ref_sub_ranges <- ref_GR[seqnames(ref_GR) == as.character(segment_i$chrom)]
        seg_strt <- GRanges(as.character(segment_i$chrom),IRanges(start = as.numeric(segment_i$loc.start),width = 1))
        check_str_overlap <- ref_sub_ranges %&% seg_strt
        #if it finds a segment that overlaps, store start overlap segment id
        if(length(check_str_overlap)>=1){
          seg_strt_id <- check_str_overlap@elementMetadata@listData$SegmentID
          mod_segment_strt <- ref_sub[which(SegmentID==seg_strt_id),] #the original segment that houses the start location of tangent seg
          if(mod_segment_strt$loc.start==segment_i$loc.start){ #if the start location exactly overlaps that of an existing segment, look at end
            #modify the segment such that it's start location is 1bp from the end of segment_i to its original end
            mod_segment_strt[,loc.start:=as.numeric(segment_i$loc.end)+1]
          }else{ #if not overlapping, the other option is that the tangent segment start falls within an existing segment or right on the edge of one i.e., starts after
            mod_segment_strt[,loc.end:=as.numeric(segment_i$loc.start)-1] #then modify it such that it ends right before the start of tangent
          }
          mod_segment_strt$seg.mean <- mean(full_dt[CONTIG==mod_segment_strt$chrom & as.numeric(START)>=mod_segment_strt$loc.start & as.numeric(END)<=mod_segment_strt$loc.end,]$log2_tangent,na.rm=T)
          mod_segment_strt[,num.mark:=((loc.end-loc.start)%/%25000)+1] #approximate/assuming all bins from loc.end -> loc.start were taken
        }else{
          seg_strt_id <- integer(0)
          mod_segment_strt <- data.table(matrix(ncol = ncol(ref_sub),nrow=0));colnames(mod_segment_strt) <- colnames(ref_sub)
        }
        #do the same for the end location of the tangent segment
        seg_end <- GRanges(as.character(segment_i$chrom),IRanges(start = as.numeric(segment_i$loc.end),width = 1))
        check_end_overlap <- ref_sub_ranges %&% seg_end
        #if it finds a segment that overlaps, store end overlap segment id
        if(length(check_end_overlap)>=1){
          seg_end_id <- check_end_overlap@elementMetadata@listData$SegmentID
          mod_segment_end <- ref_sub[which(SegmentID==seg_end_id),]
          if(mod_segment_end$loc.end==segment_i$loc.end){ #if the ends overlap exactly
            #modify the segment such that it's end location is 1bp inward from the start of segment_i
            mod_segment_end[,loc.end:=as.numeric(segment_i$loc.start)-1]
          }else{
            mod_segment_end[,loc.start:=as.numeric(segment_i$loc.end)+1]
          }
          mod_segment_end$seg.mean <- mean(full_dt[CONTIG==mod_segment_end$chrom & as.numeric(START)>=mod_segment_end$loc.start & as.numeric(END)<=mod_segment_end$loc.end,]$log2_tangent,na.rm=T)
          mod_segment_end[,num.mark:=((loc.end-loc.start)%/%25000)+1] #approximate/assuming all bins from loc.end -> loc.start were taken
        }else{
          seg_end_id <- integer(0)
          mod_segment_end <- data.table(matrix(ncol = ncol(ref_sub),nrow=0));colnames(mod_segment_end) <- colnames(ref_sub)
        }
        #take the range of the recovered&original segments and remove all segments from the
        #original segment profile whose start and end locations fall within this new range
        left_bound <- min(mod_segment_strt$loc.start,mod_segment_end$loc.start,segment_i$loc.start)
        right_bound <- max(mod_segment_strt$loc.end,mod_segment_end$loc.end,segment_i$loc.end)
        remove_segs <- as.vector(segments.ar[chrom==segment_i$chrom&(loc.start>=left_bound & loc.end<=right_bound),]$SegmentID)
        #also take the segments where the breakpoints overlapped into account
        bkpt_overlaps <- ifelse(isempty(seg_strt_id)|isempty(seg_end_id),
                                c(seg_strt_id,seg_end_id),
                                c(seg_strt_id:seg_end_id)) #retain just the non-empty value or generate the range of values
        remove_segs <- unique(sort(c(remove_segs,bkpt_overlaps)))
        if(!(isempty(remove_segs))){
          segments.ar <-segments.ar[-remove_segs,]
        }
        segments.ar <- rbind(segments.ar,mod_segment_strt,mod_segment_end,segment_i[,c('ID','chrom','loc.start','loc.end','num.mark','seg.mean','SegmentID','seg_nature')])
        segments.ar <- segments.ar[loc.end>loc.start,] #only retain segments that have reasonable start and end locations
        segments.ar <- unique(segments.ar,by=c('chrom','loc.start','loc.end'))
        segments.ar <- segments.ar[order(factor(chrom,levels = canonical_contigs),as.numeric(loc.start))]
        segments.ar[,SegmentID:=1:nrow(segments.ar)]
      }
    }
  }else{
    segments.ar <- segments.ar
  }


  if(is_male){
    X_segs <- segments.tan[chrom=='chrX'] %>%
      select(any_of(colnames(segments.ar)))
    X_segs[,seg_nature:='recovered'] #because for the male, we recover every segment from Tangent
    segments.ar <- rbind(segments.ar[-which(chrom=='chrX'),],X_segs)
  }

  #save the corrected version of segments.ar as the final segmentation output
  segments.ar <- segments.ar[,c('ID','chrom','loc.start','loc.end','num.mark','seg_nature','seg.mean')]
  segfile <- paste0(case_id,".seg.txt")
  seg_file <- segments.ar[,SegmentID:=NULL]
  seg_file <- seg_file[order(factor(chrom,levels = canonical_contigs),as.numeric(loc.start))]
  write.table(seg_file, segfile, row.names = F,quote = F,col.names = T,sep = '\t')

  new.ar.plotdf <- rbindlist(lapply(1:nrow(segments.ar),get_segdf,segments.ar))
  new.ar.plotdf <- cbind(new.ar.plotdf,segments.ar[,'seg_nature'])

  #generate final ICONICC plot
  ar_new <- ggplot(data=full_dt, aes(cum_start, log2_tangent)) +
    theme_classic()+
    ylim(-2,2) +ggtitle(paste0('AR segments corrected\n',nrow(new.ar.plotdf),' Segments'))+
    geom_rect(data=strip_data, aes(xmin=xmin,xmax=xmax, ymin=ymin, ymax=ymax,fill=fill),
              inherit.aes=F, alpha = 0.1,show.legend = F)+
    geom_segment(data = new.ar.plotdf,
                 mapping = aes(x=x, y=y, xend=x.end, yend=y.end,col=seg_nature),
                 inherit.aes = FALSE,lwd=1)+scale_color_manual(values = c('recovered'='blue','original'='red'))+
    theme(plot.title = element_text(hjust = 0.5)) +xlab("Start Position") + ylab("Segment Mean")+
    guides(color = guide_legend(title = "Segment Identity")) +theme(legend.position = "none")

  requiredPlots <- paste0(case_id,'_segmentPlots.pdf')
  pdf(file=requiredPlots,width = 16,height = 12)
  par(mar=c(2,2,0.5,0.5))
  grid.arrange(ar_old,tangent_plot,ar_new,
               ncol = 1,nrow=3,widths=15,heights=c(4,4,4),
               left = 0, right = 0, padding = unit(0, "line"))
  dev.off()

  allelic_counts_file <- paste0(case_id,"_processed_counts.txt")
  accts_file <- full_dt[order(factor(CONTIG,levels = canonical_contigs),as.numeric(START))]
  write.table(accts_file, allelic_counts_file, row.names = F,quote = F,col.names = T,sep = '\t')

  #store the number of segments generated due to each method. [To be excluded in the final version]
  seg_no <- matrix(data=c(nrow(tan_plotdf),nrow(ar.plotdf),nrow(new.ar.plotdf),nrow(segments.ar[seg_nature=='recovered',])),ncol = 4, nrow = 1,
                   dimnames = list(case_id,c('Tangent','Allelic Ratio','Allelic Ratio post-recovery','Recovered Segments')))
  segments_meth <- paste0(case_id,"_segment_cts.txt")
  write.table(seg_no, segments_meth, row.names = T,quote = F,col.names = T,sep = '\t')
  print(seg_no)
}
