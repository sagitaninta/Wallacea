# How to plot results from angsd realSFS
# http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/angsd-activity-sfs-fst-pbs/
# https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md
# Mandatory library and auxiliary functions
nnorm <- function(x) x/sum(x)
library(tidyverse)

#######################################################3
########### FUNCTIONS FOR SINGLE SFS FILES #############
#######################################################3

read_sfs_freq<-function(f){
  x<-read.delim(f,header=F,sep = " ",strip.white =T)
  x<-x[,colSums(!is.na(x))==nrow(x)]
  x<-x[-1]
  colnames(x)<-1:length(x)
  x<-x/sum(x)
  x<-x[,-ncol(x)]
  x$region<-as.factor(rownames(x))
  # change into a data frame
  x_df<-gather(x, key="sites",
               value = "freq",
               -region)
}

read_sfs_count<-function(f){
  x<-read.delim(f,header=F,sep = " ",strip.white =T)
  x<-x[,colSums(!is.na(x))==nrow(x)]
  x<-x[-1]
  colnames(x)<-1:length(x)
  return(x)
}

#######################################################3
########### FUNCTIONS FOR UNFOLDED SFS #################
#######################################################3

# Getting the raw site count values in one dataframe
siteCounts_df<-function(l,taxa,island){
  x<-lapply(l,scan)
  x_names<-unlist(regmatches(l,gregexpr("[A-Za-z]*_pos",l)))
  for(i in 1:length(x)){
    names(x)[[i]] <- (x_names[i])
  }
  df<-do.call(rbind,x)
  #df<-df[,-1] # removing number of sites with 0 derived
  #df<-df[,-ncol(df)] # count of max allele
  df<-t(rbind(df,polymorphicSites=seq(1,ncol(df))))
  df<-as.data.frame(df) %>% gather(key="mutation",
                                   value="count",x_names)
  df$taxa<-rep(taxa,nrow(df))
  df$island<-rep(island,nrow(df))
  return(df)
}

# Getting all sites frequency
allAF_df<-function(l, taxa, island){
  x<-lapply(l,scan)
  x_names<-unlist((regmatches(l,gregexpr("[A-Za-z]*_pos",l))))
  for(i in 1:length(x)){
    names(x)[[i]] <- (x_names[i])
  }
  df<-do.call(rbind,x)
  # we take counts of sites with derived allele frequency 0 to 2N
  colnames(df)<-0:(ncol(df)-1)
  # count the proportion of those sites
  df<-t(apply(df,1,nnorm))
  df<-t(rbind(df,derivedAllele=1:(ncol(df)-1)))
  df<-as.data.frame(df) %>% gather(key = "mutation",
                                   value = "alleleFrequency",
                                   x_names)
  df$taxa<-rep(taxa,nrow(df))
  df$island<-rep(island,nrow(df))
  return(df)
}

# Getting derived or non-ancestral sites
derivedAF_df<-function(l, taxa, island){
  x<-lapply(l,scan)
  x_names<-unlist((regmatches(l,gregexpr("[A-Za-z]*_pos",l)))) # regex needs to be changed according to file name
  for(i in 1:length(x)){
    names(x)[[i]] <- x_names[i]
  }
  df<-do.call(rbind,x)
  df<-df[,-1]
  # we take counts of sites with derived allele frequency 1 to 2N
  colnames(df)<-1:ncol(df)
  # count the proportion of those sites
  df<-t(apply(df,1,nnorm))
  df<-t(rbind(df,derivedAllele=1:ncol(df)))
  df<-as.data.frame(df) %>% gather(key = "mutation",
                                   value = "derivedAlleleFrequency",
                                   x_names)
  df$taxa<-rep(taxa,nrow(df))
  df$island<-rep(island,nrow(df))
  return(df)
}

# Getting polymorphic sites

# Functions and objects to get SFS values from 1D SFS, 
# Only LoF and synonymous sites --------------------------------
# Getting derived sites
derivedSFS_df<-function(x,y){
  data <- rbind(
    LOF=scan(x)[-1],
    SYN=scan(y)[-1]
  )
  colnames(data) <- 1:ncol(data)
  data<-t(apply(data, 1, nnorm)) # note that 2N sites is used here
  df<-data.frame(LOF=t(data)[,1],
                 SYN=t(data)[,2],
                 non_ancestral=seq(1,ncol(data)))
  df<-df %>% gather(key = "mutation",
                    value = "frequency",
                    LOF,SYN)
  return(df)
}

# Getting polymorphic sites for two categories only
poly_df<-function(x,y){
  data <- rbind(
    LOF=scan(x)[-1],
    SYN=scan(y)[-1]
  )
  colnames(data) <- 1:ncol(data)
  # note that compared to previous function, 
  # this includes ony 2N-1 sites
  data<-t(apply(data[,-ncol(data)],1,nnorm)) 
  df<-data.frame(LOF=t(data)[,1],
                 SYN=t(data)[,2],
                 polymorphicSites=seq(1,ncol(data)))
  df<-df %>% gather(key = "mutation",
                    value = "frequency",
                    LOF,SYN)
  return(df)
}

# More categories for polymorphic sites -------------------- 
# (the category names must ended with suffix _pos)
poly_df2<-function(l,taxa,island){
  x<-lapply(l,scan)
  x_names<-unlist(regmatches(l,gregexpr("[A-Za-z]*_pos",l)))
  for(i in 1:length(x)){
    names(x)[[i]] <- (x_names[i])
  }
  df<-do.call(rbind,x)
  df<-df[,-1]
  colnames(df)<-1:ncol(df)
  # removing last allele category
  df<-t(apply(df,1,nnorm)) 
  df<-t(apply(df[,-ncol(df)],1,nnorm)) 
  df<-t(rbind(df,polymorphicSites=seq(1,ncol(df))))
  df<-as.data.frame(df) %>% gather(key="mutation",
                                   value="frequency",x_names)
  df$taxa<-rep(taxa,nrow(df))
  df$island<-rep(island,nrow(df))
  return(df)
}

######################################################
################## FOR FOLDED SFS ###################3
#####################################################3

# for folded sfs, 23 diploid individuals of buton anoa will have 46 chr
# thus, 47 categories, but 0 and 47 got collapsed, 1 and 46, 2 and 45 etc
# so 23 categories
# read here: https://github.com/ANGSD/angsd/issues/321
poly_folded_df<-function(x,y){
  data<-rbind(
    LOF=scan(x)[-1],
    SYN=scan(y)[-1]
  )
  data<-rbind(
    data[1,][which(data[1,]>0)],
    data[2,][which(data[2,]>0)]
  )
  colnames(data) <- 1:ncol(data)
  data<-t(apply(data[,-ncol(data)],1,nnorm))
  df<-data.frame(LOF=t(data)[,1],
                 SYN=t(data)[,2],
                 polymorphicSites=seq(1,ncol(data)))
  df<-df %>% gather(key = "mutation",
                    value = "frequency",
                    LOF,SYN)
  return(df)
}  

# Expanding the poly_folded function to accommodate the incorporation of more sites --------------------------
poly_folded_df2<-function(l,taxa,island){
  x<-lapply(l,scan)
  x_names<-unlist(regmatches(l,gregexpr("[A-Za-z]*_pos",l)))
  for(i in 1:length(x)){
    names(x)[[i]] <- x_names[i]
  }
  df<-do.call(rbind,x)
  df<-df[,-1]
  df<-df[,which(colSums(df)>0)]
  colnames(df)<-1:ncol(df)
  df<-t(apply(df[,-ncol(df)],1,nnorm))
  df<-t(rbind(df,polymorphicSites=seq(1,ncol(df))))
  df<-as.data.frame(df) %>% gather(key="mutation",
                                   value="frequency",x_names)
  df$taxa<-rep(taxa,nrow(df))
  df$island<-rep(island,nrow(df))
  return(df)
}  

##################################################3
########### FUNCTIONS TO DOWNSAMPLE ###############
##################################################3

#down sample to 5 individuals (10 chromosome) and exclude fixed derived
downsampleSFS <- function(x,chr){ #x 1:2n , chr < 2n
  n<-length(x) # first you store the original number of chr of the pop
  mat <- sapply(1:chr,function(i) choose(1:n,i)*choose(n-(1:n),chr-i)/choose(n,chr))
  nnorm( as.vector(t(mat) %*% x)[-chr] )
}

downsampleSFSderived <- function(x,chr){ #x 1:2n , chr < 2n
  n<-length(x) # first you store the original number of chr of the pop
  mat <- sapply(1:chr,function(i) choose(1:n,i)*choose(n-(1:n),chr-i)/choose(n,chr))
  nnorm( as.vector(t(mat) %*% x) )
}

# downsample polymorphic sites from folded SFS
poly_folded_downsampled_df2<-function(l,taxa,island,n_chr){
  x<-lapply(l,scan)
  x_names<-unlist(regmatches(l,gregexpr("[A-Za-z]*_pos",l)))
  for(i in 1:length(x)){
    names(x)[[i]] <- (x_names[i])
  }
  df<-do.call(rbind,x)
  df<-df[,-1]
  df<-df[,which(colSums(df)>0)]
  colnames(df)<-1:ncol(df)
  df<-t(apply(df[,-ncol(df)],1,nnorm))
  df<-t(rbind(df,polymorphicSites=seq(1,ncol(df))))
  df<-as.data.frame(df) %>% gather(key="mutation",
                                   value="frequency",x_names)
  mat<-df[,1:3] %>% spread(key=polymorphicSites,value=frequency)
  mat_ds<-t(apply(mat[,-1],1,downsampleSFS,chr=n_chr))
  rownames(mat_ds)<-mat[,1]
  mat_ds<-t(rbind(mat_ds,polymorphicSites=seq(1,ncol(mat_ds))))
  mat_df<-as.data.frame(mat_ds) %>% gather(key="mutation",
                                           value="frequency",x_names)
  mat_df$taxa<-rep(taxa,nrow(mat_df))
  mat_df$island<-rep(island,nrow(mat_df))
  return(mat_df)
}

# Testing function
# for folded sfs, 8 diploid individuals of togian babirusa will have 16 chr
# thus, 17 categories
derivedSFS_downsampled_df<-function(l,taxa,island,n_chr){
  x<-lapply(l,scan)
  x_names<-unlist(regmatches(l,gregexpr("[A-Za-z]*_pos",l)))
  for(i in 1:length(x)){
    names(x)[[i]] <- (x_names[i])
  }
  df<-do.call(rbind,x)
  df<-df[,-1]
  df<-t(apply(df,1,nnorm))
  dfDown<-t(apply(df,1,downsampleSFSderived,chr=n_chr))
  dfDown<-t(rbind(dfDown,polymorphicSites=seq(1,ncol(dfDown))))
  dfDown<-as.data.frame(dfDown) %>% gather(key="mutation",
                                           value="frequency",x_names)
  dfDown$taxa<-rep(taxa,nrow(dfDown))
  dfDown$island<-rep(island,nrow(dfDown))
  return(dfDown)
}

derivedFreqSFS_downsampled_df<-function(l,island,n_chr){
  x<-lapply(l,scan)
  x_names<-unlist(regmatches(l,gregexpr("[A-Za-z]*_pos",l)))
  for(i in 1:length(x)){
    names(x)[[i]] <- (x_names[i])
  }
  df<-do.call(rbind,x)
  df<-df/sum(df) # Ricard's idea
  df<-df[,-1]
  dfDown<-t(apply(df,1,downsampleSFS,chr=n_chr))
  dfDown2<-t(rbind(dfDown,polymorphicSites=seq(1,ncol(dfDown))))
  dfDown3<-as.data.frame(dfDown2) %>% gather(key="mutation",
                                           value="frequency",x_names)
  dfDown3$island<-rep(island,nrow(dfDown3))
  return(dfDown3)
}
