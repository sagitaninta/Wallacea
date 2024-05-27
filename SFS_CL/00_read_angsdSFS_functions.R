# Library

library(tidyverse)

# Function to read ANGSD SFS

#####################################3
########## UNFOLDED SFS #############
####################################3

# Angsd unfolded SFS give us 2N+1 output from a population of diploid size N
# First column is invariant sites
# Last column is fixed sites

# reading all SFS counts means 
# reading observed count
read_sfs_count<-function(file){
  v<-scan(file)
  return(v)
}

# reading derived allele count
read_sfs_count_derived<-function(file){
  v<-scan(file)[-1] # excluding invariant sites
  return(v)
}

# reading polymorphic sites only
read_sfs_count_polymorphic<-function(file){
  v<-scan(file)[-1]
  v<-v[-(length(v))] # excluding fixed sites
  return(v)
}

# reading all SFS frequency usually means only considering derived counts
# function to calculate frequency
nnorm <- function(x) x/sum(x)

# reading derived counts
read_sfs_freq_derived<-function(file){
  v<-scan(file)[-1]
  v<-nnorm(v) # densiy of SNPs
  return(v)
}

read_sfs_freq_polymorphic<-function(file){
  v<-scan(file)[-1]
  v<-nnorm(v)
  v<-v[-(length(v))] #excluding fixed sites
  return(v)
}

# When used on a list of files of SFS, these functions can be used with lapply
# Example
readList_sfsDerivedFreq<-function(l){
  x<-lapply(l,read_sfs_freq_derived)
  x_names<-unlist(regmatches(l,gregexpr("[A-Za-z]*_pos",l)))
  for(i in 1:length(x)){
    names(x)[[i]] <- (x_names[i])
  }
  df<-do.call(rbind,x)
  return(df)
}

readList_sfsDerivedCount<-function(l){
  x<-lapply(l,read_sfs_count_derived)
  x_names<-unlist(regmatches(l,gregexpr("[A-Za-z]*_pos",l)))
  for(i in 1:length(x)){
    names(x)[[i]] <- (x_names[i])
  }
  df<-do.call(rbind,x)
  return(df)
}
