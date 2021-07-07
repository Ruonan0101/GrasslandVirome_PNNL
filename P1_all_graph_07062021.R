######Fig2
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(scales)
d <-read.csv("~/Desktop/2021_five_Manuscripts/PNNL_Manuscript_Revision/P1_ThreeGrassland_replicateMG_revision/IA_KS_WA_boxplot_normalized_more.csv", header = TRUE)
stat.test <- d %>%
  group_by(Variable) %>%
  t_test(Value ~ Location) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test_new <- stat.test %>% add_xy_position(x = "Location")
p.adj_new<- scientific(stat.test_new$p.adj, digits = 3)
M_abd<-read.csv("~/Desktop/2021_five_Manuscripts/PNNL_Manuscript_Revision/P1_ThreeGrassland_replicateMG_revision/MicrobialBiomass.csv", header = TRUE)
plt1<-ggplot(M_abd, aes(x = Location, y = log(ValueG))) +
  geom_boxplot(aes(y = log(ValueG))) +
  geom_point(aes(y = log(ValueG), color=Identity), position=position_jitterdodge())+
  scale_y_continuous(name = "Microbial biomass(ng DNA/g soil)")+
  scale_fill_manual(values = c("#0000e7"))+
  scale_color_manual(values = c("#0000e7"))+
  theme_bw()+
  theme(axis.text = element_text(size=9),
        axis.title.x = element_blank(),
        axis.title.y=element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = 10)),
        legend.title = element_blank(), legend.position = 'NA',
        plot.margin = unit(c(1, 1, 0, 1), "lines"))
plt1


M_div<-read.csv("~/Desktop/2021_five_Manuscripts/PNNL_Manuscript_Revision/P1_ThreeGrassland_replicateMG_revision/MicrobialDiversity.csv", header = TRUE)
plt2<-ggplot(M_div, aes(x = Location, y = log(ValueG))) +
  geom_boxplot(aes(y = log(ValueG))) +
  geom_point(aes(y = log(ValueG), color=Identity), position=position_jitterdodge())+
  scale_y_continuous(name = "NO. of 16S rRNA gene clusters")+
  scale_fill_manual(values = c("#0000e7"))+
  scale_color_manual(values = c("#0000e7"))+
  theme_bw()+
  theme(axis.text = element_text(size=9),
        axis.title.x = element_blank(),
        axis.title.y=element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = 10)),
        legend.title = element_blank(), legend.position = 'NA',
        plot.margin = unit(c(1, 1, 0, 1), "lines"))
plt2

d_activity <-read.csv("~/Desktop/Project/P1_Ruonan_Virome/Revision_replicate_dir/Percent_Spacer_hit_VC_normalized_updated.csv", header = TRUE)
plt3<-ggplot(d_activity, aes(x = Location, y = Value*100)) +
  geom_boxplot(aes(y = Value*100)) +
  geom_point(aes(y = Value*100, color=Identity), position=position_jitterdodge())+
  scale_y_continuous(name = "%spacers that hit viral contigs")+
  scale_fill_manual(values = c("#0000e7"))+
  scale_color_manual(values = c("#0000e7"))+
  theme_bw()+
  theme(axis.text = element_text(size=9),
        axis.title.x = element_blank(),
        axis.title.y=element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = 10)),
        legend.title = element_blank(), legend.position = 'NA',
        plot.margin = unit(c(1, 1, 0, 1), "lines"))

plt3

V_abd<-read.csv("~/Desktop/2021_five_Manuscripts/PNNL_Manuscript_Revision/P1_ThreeGrassland_replicateMG_revision/ViralAbundance.csv", header = TRUE)
plt4<-ggplot(V_abd, aes(x = Location, y = log(ValueG))) +
  geom_boxplot(aes(y = log(ValueG))) +
  geom_point(aes(y = log(ValueG), color=Identity), position=position_jitterdodge())+
  scale_y_continuous(name = "Read coverage of viral contigs/g soil")+
  scale_fill_manual(values = c("#E7B800"))+
  scale_color_manual(values = c("#E7B800"))+
  theme_bw()+
  theme(axis.text = element_text(size=9),
        axis.title.x = element_blank(),
        axis.title.y=element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = 10)),
        legend.title = element_blank(), legend.position = 'NA',
        plot.margin = unit(c(1, 1, 0, 1), "lines"))
plt4


V_div<-read.csv("~/Desktop/2021_five_Manuscripts/PNNL_Manuscript_Revision/P1_ThreeGrassland_replicateMG_revision/ViralDiversity.csv", header = TRUE)
plt5<-ggplot(M_div, aes(x = Location, y = log(ValueG))) +
  geom_boxplot(aes(y = log(ValueG))) +
  geom_point(aes(y = log(ValueG), color=Identity), position=position_jitterdodge())+
  scale_y_continuous(name = "NO. of viral clusters")+
  scale_fill_manual(values = c("#E7B800"))+
  scale_color_manual(values = c("#E7B800"))+
  theme_bw()+
  theme(axis.text = element_text(size=9),
        axis.title.x = element_blank(),
        axis.title.y=element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = 10)),
        legend.title = element_blank(), legend.position = 'NA',
        plot.margin = unit(c(1, 1, 0, 1), "lines"))
plt5

d_lysogeny <-read.csv("~/Desktop/2021_five_Manuscripts/PNNL_Manuscript_Revision/P1_ThreeGrassland_replicateMG_revision/lysogeny_normalized.csv", header = TRUE)
plt6<-ggplot(d_lysogeny, aes(x = Location, y = log(ValueG))) +
  geom_boxplot(aes(y = log(ValueG))) +
  geom_point(aes(y = log(ValueG), color=Identity), position=position_jitterdodge())+
  scale_y_continuous(name = "Read coverage of lysogenic markers/g soil")+
  scale_fill_manual(values = c("#E7B800"))+
  scale_color_manual(values = c("#E7B800"))+
  theme_bw()+
  theme(axis.text = element_text(size=9),
        axis.title.x = element_blank(),
        axis.title.y=element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = 10)),
        legend.title = element_blank(), legend.position = 'NA',
        plot.margin = unit(c(1, 1, 0, 1), "lines"))
plt6

library(gridExtra)
gridExtra::grid.arrange(plt1, plt2,plt3, plt4, plt5, plt6, nrow=2)


###fig 3b, viral composition heatmap clustering
library(pheatmap)
data<-read.delim("/Users/wuru978/Desktop/Project/P1_Ruonan_Virome/Revision_replicate_dir/ViralClusterComposition.txt",header=T, row.names="gene")
cal_z_score <- function(x){
  (log(x+1))
}
data_norm <- t(apply(data, 1, cal_z_score))
my_sample_col <- data.frame(sample = rep(c("IA", "KS","WA"), c(3,3,3)))
row.names(my_sample_col) <- colnames(data_norm)
sample        <- c("#FF7D21", "#8CC4F3", "#C8ACD3")
names(sample) <- c("IA","KS","WA")
anno_colors <- list(sample = sample)
pheatmap(data_norm,
         annotation=my_sample_col,
         annotation_colors = anno_colors,
         cutree_cols = 6,fontsize = 6,cluster_rows=FALSE, show_colnames =FALSE, show_rownames = FALSE)

##fig3c.differential abundances of VC-common/shared cluster, pie chart
library(plotrix)
plt1_pieval<-c(27,18)
plt1_pielabels<-
  c("Not significant\n60.0%","Significant\n40.0%")
plt1<-pie3D(plt1_pieval,radius=0.9,labels=plt1_pielabels,explode=0.1, 
            col=c("#999999","#FF3333"), labelcex = 0.6
            # ,labelpos=lp
)

plt2_pieval<-c(29,26)
plt2_pielabels<-
  c("Not significant\n52.7%","Significant\n47.3%")
plt2<-pie3D(plt2_pieval,radius=0.9,labels=plt2_pielabels,explode=0.1, 
            col=c("#999999","#FF3333"), labelcex = 0.6)

###fig3d: crispr hit, %contigs 
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(scales)
d<- read.csv("~/Desktop/Project/P1_Ruonan_Virome/Revision_replicate_dir/Common_SiteSpecific_oercentGotCRISPR_forR_removingSingleton_outlier.csv", header = T)
stat.test <- d %>%
  t_test(percent_gotCRISPR ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test_new <- stat.test %>% add_xy_position(x = "group")
p.adj_new<- scientific(stat.test_new$p.adj, digits = 3)
ggplot(d, aes(x=group, y=percent_gotCRISPR, color=group)) + 
  geom_boxplot()+
  theme_bw()+
  theme(axis.title.x = element_blank(), legend.position = 'none', axis.title.y=element_text(size=9), axis.text.x=element_text(size=9))+
  ylab('% of contigs with CRISPR\nspacer hits per viral cluster')+
  stat_pvalue_manual(stat.test_new, hide.ns = TRUE, y.position = 1)

###fig4 host VC pairing
library(ggalluvial)
library(ggplot2)
library(ggrepel)
library(grid)
library(RColorBrewer)
library(colorRamps)
d<-read.csv("/Users/wuru978/Desktop/Project/P1_Ruonan_Virome/Revision_replicate_dir/host_VC_location_CommonOrNot_Alluvial_forR_new.csv",check.names=FALSE)
p<-ggplot(data = d,aes(axis1 = Host, axis3=VC)) +
  scale_x_discrete(limits = c("Host","Viral cluster"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = Host), alpha=0.9) +
  geom_stratum(width=1/12, alpha=.5, color='white') 
p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  facet_wrap(Location~Cat, scale='free')+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(27))+
  theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"), 
        axis.text.x = element_text(size = 8, color='black'), 
        panel.spacing=unit(0.1, "lines"))+
  guides(fill=guide_legend(ncol=1))

#######fig5 AMG composition heatmap clustering
library(pheatmap)
library(ggplot2)
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
p<-pheatmap(data_norm,annotation_row=my_gene_col,
            annotation_col = my_sample_col,
            annotation_colors = anno_colors, 
            cutree_cols = 6,fontsize = 8,
            cellwidth=12, cellheight=8, 
            cluster_rows = F)
p
p+theme(legend.position='bottom')

####fig.6 SEM
library(lavaan)
library(semPlot)
d<-read.table("~/Desktop/2021_five_Manuscripts/PNNL_Manuscript_Revision/P1_ThreeGrassland_replicateMG_revision/SEM_selectedF_renamed.txt", sep="\t", header=T)
mydata<-na.omit(d)
Model3<-'
AMG+Lysogeny~Precipitation
Vabd+AMG~Lysogeny
Vabd~VHI
VHI~AMG
Mbio~AMG+VHI
OM~Mbio+AMG
'
#  chisq     df pvalue    gfi    cfi    rmr   srmr  rmsea 
#7.164 10.000  0.710  0.834  1.000  3.252  0.042  0.000 
Model3Fit <- sem(Model3, data=mydata, fixed.x=F, check.gradient = FALSE)
fitMeasures(Model3Fit, c("chisq","df","pvalue","gfi","cfi","rmr","srmr","rmsea"))
semPaths(Model3Fit, what='std',layout='tree3', ask = FALSE,
         edge.label.cex=0.9,sizeMan=10,label.font=30,
         posCol = c("blue", "red"),
         fade=F, nCharNodes=0, intercepts = F,residuals=F, thresholds = F, legend=FALSE)
table<-parameterEstimates(Model3Fit,standardized=TRUE)
table<-table[!table$lhs==table$rhs,]
b<-gettextf('%.2f\np=%.2f', table$std.all, digits=table$pvalue)
lbls <- c("Precipitation", "NH4", "Lysogeny", "pH", 
          "VHI", "VADF", 'Biomass','OM', 'AMG','ENV','Biomass_adj','Vabd','Mbio')
grps <- list(Abiotic = c("NH4", "pH", 'OM', "Precipitation",'ENV'), 
             microbe=c('Biomass','Biomass_adj','VHI','Mbio'),
             #FreeV=c('VHI'),
             #LysoV=c("Lysogeny", 'AMG'),
             virus = c("VADF",'Vabd',"Lysogeny", 'AMG')
)
semPaths(Model3Fit, what='std',layout = 'tree2', ask = FALSE,
         edge.label.cex=0.8,sizeMan=9,
         edgeLabels = b,label.font=50,edge.label.position=0.5,
         groups = grps, color = c("#dfc27d", 
                                  '#e66101',
                                  # '#6495ED',
                                  '#99b8f2'),
         posCol = c("blue", "red"),
         fade=F, nCharNodes=0, intercepts = F,residuals=F, thresholds = F, legend=FALSE)

####Suppl_file_1:boxplot_soil chem
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(scales)
setwd("~/Desktop/Project/P1_Ruonan_Virome/Revision_replicate_dir/")
d <-read.csv("WA_IA_KS_soilChem.csv", header = TRUE)
stat.test <- d %>%
  group_by(Variable) %>%
  t_test(Value ~ Location) %>%
  adjust_pvalue(method = "BY") %>%
  ####"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"
  add_significance()%>%
  filter(p.adj < 0.05)
stat.test <- d %>%
  group_by(Variable) %>%
  pairwise_t_test(
    Value ~ Location, paired = TRUE, 
    p.adjust.method = "bonferroni"
  ) %>%
  select(-df, -statistic, -p) # Remove details
stat.test
write.csv(stat.test, 'SoilChem_stat_pairt-test.csv')
stat.test
stat.test <- stat.test %>% add_xy_position(x = "Location")
stat.test
stat_write <- apply(stat.test,2,as.character)
write.csv(stat_write, 'SoilChem_stat.csv')
cbbPalette <- c("#FF7D21", "#8CC4F3", "#C8ACD3")

bxp <- ggboxplot(d, x = "Location", y = "Value", fill = 'Location', ggtheme = theme_pubr(border = TRUE))+
  facet_wrap(~Variable,scales = "free", strip.position='left', ncol = 5)+
  scale_fill_manual(values=cbbPalette)

bxp
p.adj_new<- scientific(stat.test$p.adj, digits = 3)
bxp + 
  #stat_pvalue_manual(stat.test, label = "{p.adj.signif}", hide.ns = TRUE,size=5,bracket.nudge.y = -155, tip.length = 0.005)+
  #scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))+
  ylab(NULL)+
  theme(axis.text = element_text(size=10),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "bottom",
        strip.placement='outside',
        axis.title.x = element_blank(),
        axis.text.x = element_blank())


####Suppl_file_2: 16s composition heatmap
library(pheatmap)
setwd('/Users/wuru978/Desktop/Project/P1_Ruonan_Virome/16s_18s_dir/')
data<-read.delim("/Users/wuru978/Desktop/Project/P1_Ruonan_Virome/16s_18s_dir/16s_abd_masterTable_forPheatmap.txt",header=T, row.names="lineage")
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