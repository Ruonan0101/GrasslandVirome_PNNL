###fig2:barchat with sig
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(scales)
setwd("~/Desktop/Project/P1_Ruonan_Virome/Revision_replicate_dir/")
d <-read.csv("WA_KS_IA_abd_count_normalized104_gramSoil_log_forR.csv", header = TRUE)

stat.test <- d %>%
  group_by(Variable) %>%
  t_test(Value ~ Location) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
stat.test <- stat.test %>% add_xy_position(x = "Location")

cbbPalette <- c("#FF7D21", "#8CC4F3", "#C8ACD3")
myplot <- ggboxplot(
  d, x = "Location", y = "Value",
  fill = "Location", palette = cbbPalette, legend = "none",add = "jitter",
  ggtheme = theme_pubr(border = TRUE)
) +
  facet_wrap(~Variable,scales = "free")+
  ylab('Log transformed value/(10^4 reads*g of soil)')+
  theme(axis.text = element_text(size=10),
        strip.background = element_blank(),
        strip.text.x = element_blank())
p.adj_new<- scientific(stat.test$p.adj, digits = 3)
myplot + stat_pvalue_manual(stat.test, label = "{p.adj_new}{p.adj.signif}", size=3, label.y=c(5.0,5.1,5.2), tip.length = 0.01)+
  scale_y_continuous(expand = expansion(mult = c(0,0.05)))

###fig 3b, viral composition heatmap clustering
library(pheatmap)
setwd('/Users/wuru978/Desktop/Project/P1_Ruonan_Virome/Revision_replicate_dir/')
data<-read.delim("ViralClusterComposition.txt",header=T, row.names="gene")
cal_z_score <- function(x){
  (log(x+1))
}
data_norm <- t(apply(data, 1, cal_z_score))
my_sample_col <- data.frame(sample = rep(c("IA", "KS","WA"), c(3,3,3)))
row.names(my_sample_col) <- colnames(data_norm)
my_sample_col
sample        <- c("#FF7D21", "#8CC4F3", "#C8ACD3")
names(sample) <- c("IA","KS","WA")
anno_colors <- list(sample = sample)
anno_colors
my_sample_col
pheatmap(data_norm,
        annotation=my_sample_col,
        annotation_colors = anno_colors,
        cutree_cols = 6,fontsize = 6,cluster_rows=FALSE)

###fig3c: crispr hit, %contigs 
setwd("~/Desktop/Project/P1_Ruonan_Virome/Revision_replicate_dir/")
d<- read.csv("Common_SiteSpecific_oercentGotCRISPR_forR_removingSingleton.csv", header = T)
head (d)
#colnames(host)<-c("vcluster",'percent_gotCRISPR',"group")
library(ggplot2)
ggplot(d, aes(x=group, y=percent_gotCRISPR, color=group)) + 
  geom_boxplot()+
  geom_jitter(position=position_jitter(0.2)) + 
  theme_bw()+
  theme(axis.title.x=element_blank())+
  theme(text = element_text(size=13))+
  ylab('% of contigs with CRISPR hit per viral cluster')

###fig4 host VC pairing
setwd('/Users/wuru978/Desktop/Project/P1_Ruonan_Virome/Revision_replicate_dir/')
library(ggalluvial)
library(ggplot2)
library(ggrepel)
library(grid)
d<-read.csv("host_VC_location_CommonOrNot_Alluvial_forR.csv",check.names=FALSE)

p<-ggplot(data = d,aes(axis1 = Host, axis3=VC)) +
  scale_x_discrete(limits = c("Host","Viral cluster"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = Host), alpha=0.9) +
  geom_stratum(width=1/12, alpha=.5, color='white') 
library(RColorBrewer)
library(colorRamps)
p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  facet_wrap(Location~Cat, scale='free')+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(27))+
  theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"), 
        axis.text.x = element_text(size = 8, color='black'), 
        panel.spacing=unit(0.1, "lines"))+
  guides(fill=guide_legend(ncol=1))

###fig5b:barchat with sig
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(scales)
setwd("~/Desktop/Project/P1_Ruonan_Virome/Revision_replicate_dir/")
d <-read.csv("Percent_Spacer_hit_VC_normalized.csv", header = TRUE)

stat.test <- d %>%
  t_test(Value ~ Location) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
stat.test <- stat.test %>% add_xy_position(x = "Location")

cbbPalette <- c("#FF7D21", "#8CC4F3", "#C8ACD3")
myplot <- ggboxplot(
  d, x = "Location", y = "Value",
  fill = "Location", palette = cbbPalette, legend = "none",add = "jitter",
  ggtheme = theme_pubr(border = TRUE)
) +
  ylab('Normalized %spacer hit to viral contigs')+
  theme(axis.text = element_text(size=10),
        strip.text = element_text(size = 10), 
        axis.title.x = element_blank())
p.adj_new<- scientific(stat.test$p.adj, digits = 3)
myplot + stat_pvalue_manual(stat.test, label = "{p.adj_new}{p.adj.signif}", size=3, label.y=c(5.0,5.1,5.2), tip.length = 0.01)+
  scale_y_continuous(expand = expansion(mult = c(0,0.05)), labels = scales::percent, breaks=c(0.00,0.05,0.10,0.15,0.20,0.25,0.30))

#######fig6 AMG composition heatmap clustering
library(pheatmap)
setwd('/Users/wuru978/Desktop/Project/P1_Ruonan_Virome/Revision_replicate_dir/')
data<-read.delim("AMG_replicatedMG_normalized_finalForPheatmap.txt",header=T, row.names="gene")
cal_z_score <- function(x){
  (log(x+1))
}
data_norm <- t(apply(data, 1, cal_z_score))
my_sample_col <- data.frame(sample = rep(c("IA", "KS","WA"), c(3,3,3)))
my_gene_col <- data.frame('Function'= rep(c('Antibiotics biosynthesis','Aromatic compound degradation','Biofilm formation',
                                     'Cellulose_synthase','Chemcal metabolism','Chitin/chitosan degradation','Fatty acid biosynthesis',
                                     'Folate biosynthesis','Fructose and mannose metabolism','Glutamine metabolism',
                                     'Glutathione metabolism','Oxidative phosphorylation',
                                     'Protein metabolism','Radical-based catalysis','Sporulation','Sulfur metabolism','UDP-glocuse production','Xylan degradation'),
                          c(2,1,1,1,1,2,1,6,5,1,1,4,1,1,1,1,3,3)))
row.names(my_gene_col) <- rownames(data_norm)
row.names(my_sample_col) <- colnames(data_norm)
my_sample_col
sample        <- c("#FF7D21", "#8CC4F3", "#C8ACD3")
names(sample) <- c("IA","KS","WA")
anno_colors <- list(sample = sample)
anno_colors
my_sample_col
pheatmap(data_norm,annotation_row=my_gene_col,annotation_col = my_sample_col,annotation_colors = anno_colors, cutree_cols = 6,fontsize = 8,cellwidth=12, cellheight=8, cluster_rows = F)

####16s composition heatmap
library(pheatmap)
setwd('/Users/wuru978/Desktop/Project/P1_Ruonan_Virome/16s_18s_dir/')
data<-read.delim("16s_abd_masterTable_forPheatmap.txt",header=T, row.names="lineage")
cal_z_score <- function(x){
  (log(x+1))
}
data_norm <- t(apply(data, 1, cal_z_score))
my_sample_col <- data.frame(sample = rep(c("WA", "KS","IA"), c(3,3,3)))
row.names(my_sample_col) <- colnames(data_norm)
my_sample_col
sample        <- c("#FF7D21", "#8CC4F3", "#C8ACD3")
names(sample) <- c("IA","KS","WA")
anno_colors <- list(sample = sample)
anno_colors
my_sample_col
pheatmap(data_norm,
         annotation=my_sample_col,
         annotation_colors = anno_colors,
         cutree_cols = 6,fontsize = 6,cellwidth=12, cellheight=6,cluster_rows=FALSE)