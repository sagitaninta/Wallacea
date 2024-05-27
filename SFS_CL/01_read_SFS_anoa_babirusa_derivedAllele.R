source("00_angsdSFSplottingFunctions.R")

bbrsNOuCO<-list.files(path="input", pattern =glob2rx("bbrs_NO*CDS*unfolded.sfs"), full.names = T, all.files = T)
bbrsSEuCO<-list.files(path="input", pattern =glob2rx("bbrs_SE*CDS*unfolded.sfs"), full.names = T, all.files = T)
bbrsTOuCO<-list.files(path="input", pattern =glob2rx("bbrs_TO*CDS*unfolded.sfs"), full.names = T, all.files = T)

anoaNOuCO<-list.files(path="input", pattern =glob2rx("anoa_NO*unfolded.sfs"), full.names = T, all.files = T)
anoaSEuCO<-list.files(path="input", pattern =glob2rx("anoa_SE*unfolded.sfs"), full.names = T, all.files = T)
anoaBTuCO<-list.files(path="input", pattern =glob2rx("anoa_BT*unfolded.sfs"), full.names = T, all.files = T)

bbrsNO_df<-derivedAF_df(bbrsNOuCO,"babirusa","NO") # 24 individuals
bbrsSE_df<-derivedAF_df(bbrsSEuCO,"babirusa","SE") # 14 individuals
bbrsTO_df<-derivedAF_df(bbrsTOuCO,"babirusa","TO") # 8 individuals

bbrsNO_df<-bbrsNO_df %>% mutate(derivedAlleleFrequencyCategory=derivedAllele/48)
bbrsSE_df<-bbrsSE_df %>% mutate(derivedAlleleFrequencyCategory=derivedAllele/28)
bbrsTO_df<-bbrsTO_df %>% mutate(derivedAlleleFrequencyCategory=derivedAllele/16)
bbrs_df<-rbind(bbrsTO_df,
               bbrsSE_df,
               bbrsNO_df)

anoaNO_df<-derivedAF_df(anoaNOuCO,"anoa","NO")
anoaSE_df<-derivedAF_df(anoaSEuCO,"anoa","SE") # 15 individuals
anoaBT_df<-derivedAF_df(anoaBTuCO,"anoa","BT") # 12 individuals
anoaNO_df<-anoaNO_df %>% mutate(derivedAlleleFrequencyCategory=derivedAllele/50)
anoaSE_df<-anoaSE_df %>% mutate(derivedAlleleFrequencyCategory=derivedAllele/30)
anoaBT_df<-anoaBT_df %>% mutate(derivedAlleleFrequencyCategory=derivedAllele/46)

anoa_df<-rbind(anoaBT_df,
               anoaSE_df,
               anoaNO_df)

SFSder<-rbind(anoa_df,bbrs_df)
SFSder<-SFSder %>% mutate(island=factor(island, levels = 
                                              c("NO","SE","BT","TO")),
                          mutation=factor(mutation, levels =
                                            c("LOW_pos","MODIFIER_pos","MODERATE_pos","HIGH_pos")),
                          pop=case_when(island=='SE'~"Southeast\nSulawesi (2n=26)",
                                        island=='BT'~"Buton (2n=46)",
                                        island=='NO'~"North\nSulawesi (2n=50)",
                                        island=='TO'~"Togean (2n=16)"))

head(SFSder)

# Binned derived allele -----
SFSder_wBin<-SFSder %>% mutate(derivedAlleleCategoryBin=case_when(
  derivedAlleleFrequencyCategory>0 & derivedAlleleFrequencyCategory<=0.15 ~ "0.1",
  derivedAlleleFrequencyCategory>0.15 & derivedAlleleFrequencyCategory<=0.25 ~ "0.2",
  derivedAlleleFrequencyCategory>0.25 & derivedAlleleFrequencyCategory<=0.35 ~ "0.3",
  derivedAlleleFrequencyCategory>0.35 & derivedAlleleFrequencyCategory<=0.45 ~ "0.4",
  derivedAlleleFrequencyCategory>0.45 & derivedAlleleFrequencyCategory<=0.55 ~ "0.5",
  derivedAlleleFrequencyCategory>0.55 & derivedAlleleFrequencyCategory<=0.65 ~ "0.6",
  derivedAlleleFrequencyCategory>0.65 & derivedAlleleFrequencyCategory<=0.75 ~ "0.7",
  derivedAlleleFrequencyCategory>0.75 & derivedAlleleFrequencyCategory<=0.85 ~ "0.8",
  derivedAlleleFrequencyCategory>0.85 & derivedAlleleFrequencyCategory<=0.95 ~ "0.9",
  derivedAlleleFrequencyCategory>0.95 ~ "1")) %>%
  mutate(derivedAlleleCategoryBin=as.numeric(derivedAlleleCategoryBin))

SFSder_wBin2<-SFSder_wBin %>% group_by(derivedAlleleCategoryBin,taxa,mutation,island) %>%
  summarise(binnedFrequency=sum(derivedAlleleFrequency)) %>% as.data.frame()
SFSder2<-full_join(SFSder,SFSder_wBin2)