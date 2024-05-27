source("02_clr_bootstrap.R")
source("01_read_SFS_anoa_babirusa_derivedAllele.R")
library(ggh4x)
library(patchwork)
theme_sb<-theme_classic() +
  theme(axis.text = element_text(size=14,colour = "black"),
        axis.title = element_text(size=16, colour = "black", face="bold"),
        axis.title.x = element_text(margin = margin(t=10,r=0,b=0,l=0)),
        axis.title.y = element_text(margin = margin(t=0,r=10,b=0,l=0)),
        strip.text = element_text(size=14),
        title = element_text(size=12),
        legend.title = element_blank(),
        legend.margin = margin(2, 1, 1, 1, unit = 'mm'),
        legend.text = element_text(size=14),
        legend.key.size = unit(7, unit = 'mm'),
        plot.margin=unit(c(0.5,0.5,1,1),"lines"))
area_col<-c("#55b7a6","#f25f5c","#f0b185")

clr_df<-clr_df %>% mutate(area=case_when(island=="NO"~"North Sulawesi",
                                         island=="SE"~"Southeast Sulawesi",
                                         island=="BT"~"Buton/Togean",
                                         island=="TO"~"Buton/Togean",
                                         TRUE~"Other")) %>% 
  mutate(area=factor(area, levels=c("North Sulawesi",
                                    "Southeast Sulawesi",
                                    "Buton/Togean",
                                    "Other")))

norm_clr_df<-clr_df %>% 
  group_by(impact,island,taxa,area) %>%
  mutate(bootstrap_n = as.character(seq(1,1000))) %>%
  ungroup() %>%
  spread(impact,clr) %>%
  mutate(normalisedHigh=high/low,
         normalisedModerate=moderate/low,
         normalisedModifier=modifier/low,
         normalisedLow=low/low)

n_clr_df<-norm_clr_df %>% gather(key="normalisedImpact",value="normalisedClr",
                                 -c(island,taxa,area,low,high,moderate,modifier,bootstrap_n)) %>%
  gather(key="impact",value="clr",
         -c(island,taxa,area,bootstrap_n,normalisedImpact,normalisedClr)) %>%
  mutate(normalisedImpact=factor(normalisedImpact, levels=c("normalisedLow",
                                                            "normalisedModifier",
                                                            "normalisedModerate",
                                                            "normalisedHigh")))

# Normalised CL ----
# anoa
a_nCL<-n_clr_df %>%
  filter(island %in% c("NO","SE","TO","BT")) %>%
  filter(taxa=="anoa") %>%
  ggplot() +
  geom_density(aes(x=normalisedClr,
                   fill=island,
                   alpha=normalisedImpact), 
               color="black",size=0.2)+
  scale_fill_manual(values=area_col, guide="none")+
  scale_color_manual(values=area_col, guide="none")+
  scale_alpha_discrete(labels=c("Low","Modifier","Moderate","High"),guide="none")+
  facet_wrap(~area, scales = "free_y")+
  scale_x_continuous(limits=c(0.5,2), breaks = c(0.5,1,1.5,2))+
  theme_sb+
  labs(x=expression(bold(normalised~mean~bootstrapped~CL)))
# babirusa
b_nCL<-n_clr_df %>%
  filter(island %in% c("NO","SE","TO","BT")) %>%
  filter(taxa=="babirusa") %>%
  ggplot() +
  geom_density(aes(x=normalisedClr,
                   fill=island,
                   alpha=normalisedImpact),
               color="black",size=0.2)+
  scale_fill_manual(values=area_col, guide="none")+
  scale_color_manual(values=area_col, guide="none")+
  scale_alpha_discrete(labels=c("Low","Modifier","Moderate","High"),guide="none")+
  facet_wrap(~area, scales = "free_y")+
  scale_x_continuous(limits=c(0.5,2), breaks = c(0.5,1,1.5,2))+
  theme_sb+
  labs(x=expression(bold(normalised~mean~bootstrapped~CL)))

# SFS plot binned -----
# anoa
a_binned<-SFSder_wBin2 %>%
  filter(island %in% c("NO","SE","BT","TO")) %>%
  filter(taxa=="anoa") %>%
  ggplot() +
  geom_col(aes(x=derivedAlleleCategoryBin,
               y=binnedFrequency,
               fill=island,
               alpha=mutation),
           color="black",
           size=0.2,
           position="dodge")+
  scale_fill_manual(values=area_col, guide="none")+
  scale_alpha_discrete(range = c(0,1),labels=c("Low","Modifier","Moderate","High"))+ 
  labs(x="derived allele frequency",
       y="proportion of SNPs")+
  facet_wrap(~(island),labeller = as_labeller(popA), ncol=3)+
  scale_x_continuous(breaks=c(0.1,0.5,1))+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
  theme_sb
# babirusa
b_binned<-SFSder_wBin2 %>%
  filter(island %in% c("NO","SE","BT","TO")) %>%
  filter(taxa=="babirusa") %>%
  ggplot() +
  geom_col(aes(x=derivedAlleleCategoryBin,
               y=binnedFrequency,
               fill=island,
               alpha=mutation),
           color="black",
           size=0.2,
           position="dodge")+
  scale_fill_manual(values=area_col,guide="none")+ 
  scale_alpha_discrete(range = c(0,1),labels=c("Low","Modifier","Moderate","High"))+
  labs(x="derived allele frequency",
       y="proportion of SNPs")+
  facet_wrap(~(island),labeller = as_labeller(popB), ncol=3)+
  scale_x_continuous(breaks=c(0.1,0.5,1))+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2))+
  theme_sb

# from https://stackoverflow.com/questions/70187493/how-can-one-control-the-number-of-axis-ticks-within-facet-wrap
equal_breaks <- function(n = 3, s = 0.05, r = 0,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    if(seq[2]-seq[1] < 10^(-r)) seq else round(seq, r)
  }
}

pdf("paper_SFSbinnedAndCL_AB_areaCol_3ticks.pdf", width=15, height = 10)
design <- 
  "
AABB
CCDD
"
(a_binned + theme(strip.text.x = element_text(size=10))) + 
  (b_binned + theme(strip.text.x = element_text(size=10))) + 
  (a_nCL+scale_y_continuous(breaks=equal_breaks(n=3))+theme(strip.text = element_blank())) + 
  (b_nCL+scale_y_continuous(breaks=equal_breaks(n=3))+theme(strip.text = element_blank())) +
  plot_layout(guides = "collect", design=design) +
  plot_annotation(tag_levels = "A")
dev.off()