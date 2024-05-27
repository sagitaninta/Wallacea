source("01_read_unfoldedSFS_list.R")

# Composite likelihood bootstrapping
# calculate composite likelihood
# https://www.stat.tamu.edu/~suhasini/teaching613/ED.pdf
calcCL <- function(expected, observed) {
  sum(observed * log(expected))
}

nnorm <- function(x) x/sum(x)

# function to get zeros in proper place while bootstrapping
getSFS<-function(counts, labels) {
  df=as.data.frame(table(counts))
  sfs=rep(0, length(labels))
  sfs[as.numeric(as.character(df$counts))] <- as.numeric(as.character(df$Freq))
  return(sfs)
}

# bootstrap-then-composite likelihood 
cls_bootstrap<-function(obs,mle,n,ds,reps){
  bs<-c()
  for (i in 1:reps){
    x<-sample(seq(1,n),prob=round(obs), replace=T, size=ds)
    x<-getSFS(x,c(1:n))
    clr<-calcCL(mle,x)
    bs<-rbind(bs,clr)
  }
  return(bs)
}

# Testing downsampling to "high" sites
downsampledSitesCLR_df<-function(obs_count,exp_freq,chr,island_name){
  obs<-as.data.frame(apply(obs_count,1,cls_bootstrap,mle=exp_freq,n=chr,ds=1000,reps=1000))
  colnames(obs)[1:4]<-c("high","low","moderate","modifier")
  df<-gather(obs,key = impact, value=clr, colnames(obs)[1:4])
  df$island<-island_name
  return(df)
}

anoaSE_exp<-read_sfs_freq_derived("input/anoa_SE_anoa_largeHIGH_AA_noOut_AAisREF_SNPs_AF_noFixed_LOW_pos_unfolded.sfs")
anoaBT_exp<-read_sfs_freq_derived("input/anoa_BT_anoa_largeHIGH_AA_noOut_AAisREF_SNPs_AF_noFixed_LOW_pos_unfolded.sfs")
anoaNO_exp<-read_sfs_freq_derived("input/anoa_NO_anoa_largeHIGH_AA_noOut_AAisREF_SNPs_AF_noFixed_LOW_pos_unfolded.sfs")

anoaSE_clrDF<-downsampledSitesCLR_df(anoaSE_count,anoaSE_exp,30,"SE")
anoaBT_clrDF<-downsampledSitesCLR_df(anoaBT_count,anoaBT_exp,46,"BT")
anoaNO_clrDF<-downsampledSitesCLR_df(anoaNO_count,anoaNO_exp,50,"NO")

anoa_clrDF<-rbind(
               anoaSE_clrDF,
               anoaBT_clrDF,
               anoaNO_clrDF)
anoa_clrDF<-anoa_clrDF %>% mutate(taxa="anoa",
                                  impact = factor(impact,levels=c("low","modifier","moderate","high")),
                                  island = factor(island, levels=c("NO","SE","BT")))

# babirus
bbrsSE_exp<-read_sfs_freq_derived("input/bbrs_SE_babi_largeHIGH_AA_noOut_AAisREF_SNPs_vep_CDS_cowOrtho_noFixed_LOW_pos_unfolded.sfs")
bbrsTO_exp<-read_sfs_freq_derived("input/bbrs_TO_babi_largeHIGH_AA_noOut_AAisREF_SNPs_vep_CDS_cowOrtho_noFixed_LOW_pos_unfolded.sfs")
bbrsNO_exp<-read_sfs_freq_derived("input/bbrs_NO_babi_largeHIGH_AA_noOut_AAisREF_SNPs_vep_CDS_cowOrtho_noFixed_LOW_pos_unfolded.sfs")

bbrsSE_clrDF<-downsampledSitesCLR_df(bbrsSE_count,bbrsSE_exp,28,"SE")
bbrsTO_clrDF<-downsampledSitesCLR_df(bbrsTO_count,bbrsTO_exp,16,"TO")
bbrsNO_clrDF<-downsampledSitesCLR_df(bbrsNO_count,bbrsNO_exp,48,"NO")

head(bbrsSE_clrDF)

bbrs_clrDF<-rbind(
              bbrsSE_clrDF,
              bbrsTO_clrDF,
              bbrsNO_clrDF)
bbrs_clrDF<-bbrs_clrDF %>% mutate(taxa="babirusa",
                                  impact = factor(impact,levels=c("low","modifier","moderate","high")),
                                  island = factor(island, levels=c("NO","SE","TO")))

clr_df<-rbind(anoa_clrDF,bbrs_clrDF)