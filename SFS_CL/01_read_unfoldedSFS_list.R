rm(list = ls())
source("00_read_angsdSFS_functions.R")
# reading anoa and count data excluding invariant and fixed sites
bbrsSEu<-list.files(path="input", pattern =glob2rx("bbrs_SE*CDS*unfolded.sfs"), full.names = T, all.files = T)
bbrsTOu<-list.files(path="input", pattern =glob2rx("bbrs_TO*CDS*unfolded.sfs"), full.names = T, all.files = T)
bbrsNOu<-list.files(path="input", pattern =glob2rx("bbrs_NO*CDS*unfolded.sfs"), full.names = T, all.files = T)

anoaSEu<-list.files(path="input", pattern =glob2rx("anoa_SE*unfolded.sfs"), full.names = T, all.files = T)
anoaBTu<-list.files(path="input", pattern =glob2rx("anoa_BT*unfolded.sfs"), full.names = T, all.files = T)
anoaNOu<-list.files(path="input", pattern =glob2rx("anoa_NO*unfolded.sfs"), full.names = T, all.files = T)

# read count
bbrsSE_count<-readList_sfsDerivedCount(bbrsSEu)
bbrsTO_count<-readList_sfsDerivedCount(bbrsTOu)
bbrsNO_count<-readList_sfsDerivedCount(bbrsNOu)

anoaSE_count<-readList_sfsDerivedCount(anoaSEu)
anoaBT_count<-readList_sfsDerivedCount(anoaBTu)
anoaNO_count<-readList_sfsDerivedCount(anoaNOu)

# read freq
bbrsSE_freq<-readList_sfsDerivedFreq(bbrsSEu)
bbrsTO_freq<-readList_sfsDerivedFreq(bbrsTOu)
bbrsNO_freq<-readList_sfsDerivedFreq(bbrsNOu)

anoaSE_freq<-readList_sfsDerivedFreq(anoaSEu)
anoaBT_freq<-readList_sfsDerivedFreq(anoaBTu)
anoaNO_freq<-readList_sfsDerivedFreq(anoaNOu)

rowSums(anoaNO_count)
rowSums(anoaSE_count)
rowSums(anoaBT_count)

rowSums(bbrsNO_count)
rowSums(bbrsSE_count)
rowSums(bbrsTO_count)
