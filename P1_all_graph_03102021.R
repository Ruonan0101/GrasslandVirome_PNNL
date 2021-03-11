######Fig2, used
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
d_abd <-read.csv("Abundance_normalized.csv", header = TRUE)
plt1<-ggplot(d_abd, aes(x = Location, y = log(ValueG))) +
  geom_boxplot(aes(y = log(ValueG), fill=Identity)) +
  geom_point(aes(y = log(ValueG), color=Identity), position=position_jitterdodge())+
  scale_y_continuous(name = "Read coverage of viral contigs", sec.axis = sec_axis(~.*10, name = "ng DNA yield"))+
  scale_fill_manual(values = c("#FFFFFF", "#FFFFFF"))+
  scale_color_manual(values = c("#000099", "#E7B800"))+
  theme_bw()+
  theme(axis.text = element_text(size=9),
        axis.title.x = element_blank(),
        axis.title.y=element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = 10)),
        legend.title = element_blank(), legend.position = 'bottom',
        plot.margin = unit(c(1, 1, 0, 1), "lines"))
d_div <-read.csv("Diversity_normalized.csv", header = TRUE)
plt2<-ggplot(d_div, aes(x = Location, y = log(ValueG))) +
  geom_boxplot(aes(y = log(ValueG), fill=Identity)) +
  geom_point(aes(y = log(ValueG), color=Identity), position=position_jitterdodge())+
  scale_y_continuous(name = "NO. of viral clusters", sec.axis = sec_axis(~.*10, name = "NO. of 16S clusters"))+
  scale_fill_manual(values = c("#FFFFFF", "#FFFFFF"))+
  scale_color_manual(values = c("#000099", "#E7B800"))+
  theme_bw()+
  theme(axis.text = element_text(size=9),
        axis.title.x = element_blank(),
        axis.title.y=element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = 10)),
        legend.title = element_blank(), legend.position = 'bottom',
        plot.margin = unit(c(1, 1, 0, 1), "lines"))
d_lysogeny <-read.csv("lysogeny_normalized.csv", header = TRUE)
plt3<-ggplot(d_lysogeny, aes(x = Location, y = log(ValueG))) +
  geom_boxplot(aes(y = log(ValueG))) +
  geom_point(aes(y = log(ValueG), color=Identity), position=position_jitterdodge())+
  scale_y_continuous(name = "Abundance of lysogenic markers")+
  scale_fill_manual(values = c("#E7B800"))+
  scale_color_manual(values = c("#E7B800"))+
  theme_bw()+
  theme(axis.text = element_text(size=9),
        axis.title.x = element_blank(),
        axis.title.y=element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = 10)),
        legend.title = element_blank(), legend.position = 'bottom',
        plot.margin = unit(c(1, 1, 0, 1), "lines"))
d_activity <-read.csv("~/Desktop/Project/P1_Ruonan_Virome/Revision_replicate_dir/Percent_Spacer_hit_VC_normalized_updated.csv", header = TRUE)
plt4<-ggplot(d_activity, aes(x = Location, y = Value*100)) +
  geom_boxplot(aes(y = Value*100)) +
  geom_point(aes(y = Value*100, color=Identity), position=position_jitterdodge())+
  scale_y_continuous(name = "%spacers that hit viral contigs")+
  scale_fill_manual(values = c("#E7B800"))+
  scale_color_manual(values = c("#E7B800"))+
  theme_bw()+
  theme(axis.text = element_text(size=9),
        axis.title.x = element_blank(),
        axis.title.y=element_text(size=10, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = 10)),
        legend.title = element_blank(), legend.position = 'bottom',
        plot.margin = unit(c(1, 1, 0, 1), "lines"))
library(gridExtra)
gridExtra::grid.arrange(plt1, plt2,plt4, plt3, nrow=1)


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

##fig3c.differential abundances of VC-common/shared cluster
library(ggplot2)
library(DESeq2)
library(dplyr)
library(tidyr)
abund_table<-read.csv("/Users/wuru978/Desktop/2021_five_Manuscripts/PNNL_Manuscript_Revision/P1_ThreeGrassland_replicateMG_revision/ViralClusterComposition_ForCorrelation_CommonOrNot_RemovingSiteSpecific.csv",row.names=1,check.names=FALSE)
abund_table<-t(abund_table)
grouping_info<-data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))
countData = round(as(abund_table, "matrix"), digits = 0)
countData<-(t(countData+1)) 
dds <- DESeqDataSetFromMatrix(countData, grouping_info, as.formula(~ X1))
#Reference:https://github.com/MadsAlbertsen/ampvis/blob/master/R/amp_test_species.R
#Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
data_deseq_test = DESeq(dds, test="Wald", fitType="parametric")
res = results(data_deseq_test, cooksCutoff = FALSE)
res_tax = cbind(as.data.frame(res), as.matrix(countData[rownames(res), ]), OTU = rownames(res))
sig = 0.05
fold = 0
plot.point.size = 2
label=T
tax.display = NULL
tax.aggregate = "OTU"
res_tax_sig = subset(res_tax, padj < sig & fold < abs(log2FoldChange))
res_tax_sig <- res_tax_sig[order(res_tax_sig$padj),]
res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig) , "Yes", "No")
res_tax$Significant[is.na(res_tax$Significant)] <- "No"
res_tax_new<-res_tax %>%separate(OTU, c("Category", "col2"), "Ruonan")
head (res_tax_new)
p1 <- ggplot(data = res_tax_new, aes(x = baseMean, y = log2FoldChange, color = Significant, shape=Category)) +
  geom_jitter(position = position_jitter(width = 0.5, height=0.5), size=1)+
  scale_x_log10() +
  scale_color_manual(values=c("black", "red")) +
  labs(x = "Mean abundance", y = "Log2 fold change")+theme_bw()
if(label == T){
  if (!is.null(tax.display)){
    rlab <- data.frame(res_tax, Display = apply(res_tax[,c(tax.display, tax.aggregate)], 1, paste, collapse="; "))
  } else {
    rlab <- data.frame(res_tax, Display = res_tax[,tax.aggregate])
  }
#  p1 <- p1 + geom_text(data = subset(rlab, Significant == "Yes"), aes(label = Display), size = 2, vjust = 1)
}
p1+theme(axis.title=element_text(size=10), legend.text = element_text(size=8), legend.title = element_text(size=8))

res_tax_sig_abund = cbind(as.data.frame(countData[rownames(res_tax_sig), ]), OTU = rownames(res_tax_sig), padj = res_tax[rownames(res_tax_sig),"padj"]) 
write.table (res_tax_sig_abund, file='/Users/wuru978/Desktop/2021_five_Manuscripts/PNNL_Manuscript_Revision/P1_ThreeGrassland_replicateMG_revision/VC_DF_sigOut_removeSiteSpecific.txt', 
             row.name=T, sep = "\t", quote = F)

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
         edge.label.cex=0.9,sizeMan=9,label.font=30,
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
                                  '#80cdc1'),
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
  adjust_pvalue(method = "BH") %>%
  add_significance()%>%
  filter(p.adj < 0.05)
stat.test
stat.test <- stat.test %>% add_xy_position(x = "Location")
stat.test
stat_write <- apply(stat.test,2,as.character)
write.csv(stat_write, 'SoilChem_stat.csv')
cbbPalette <- c("#FF7D21", "#8CC4F3", "#C8ACD3")
bxp <- ggboxplot(d, x = "Location", y = "Value", fill = 'Location', add = "jitter", ggtheme = theme_pubr(border = TRUE))+
  facet_wrap(~Variable,scales = "free", strip.position='left', ncol = 5)+
  scale_fill_manual(values=cbbPalette)
bxp
p.adj_new<- scientific(stat.test$p.adj, digits = 3)
bxp + stat_pvalue_manual(stat.test, label = "{p.adj.signif}", hide.ns = TRUE,size=5,bracket.nudge.y = -95, tip.length = 0.005)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))+
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