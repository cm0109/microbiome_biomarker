# Using Caries data for researching methods for microbiome-based biomarker discovery

# Load required libraries
library(reshape2)
library(vegan)
library(ggplot2)
library(plyr)
library(phyloseq)
library(DESeq2); packageVersion("DESeq2") # ‘1.28.1’
library(metagenomeSeq); packageVersion("metagenomeSeq") # ‘1.28.1’


# Data source: 16S V1-V3 Sequencing of Severe Early Childhood Caries Affected Subjects
# Data subset for only Caries Cavity and Intact Enamel Samples

# Read species-level counts data
caries_counts <- readRDS("data/v13_car_ed_counts.rds") #Load counts data
caries_meta <- readRDS("data/meta_c_ed.rds")  #Load metadata

dim(caries_counts) # 68 365



# Visualize Diversity 
# Convert to rel abundance
caries_pct <- decostand(caries_counts, method = "total")

# Compute NMDS
set.seed(12345); capture.output(caries_pct.mds <- metaMDS(caries_pct, trymax = 200, autotransform = F, wascores = F))

# Permanova
set.seed(12345); caries_pct.perm <- adonis(formula = caries_pct ~ caries_meta$group, strata = caries_meta$subject) # R2 0.19607, p < 0.001

# Making dataframe for plotting
caries_pct.mds.df <- data.frame(scores(caries_pct.mds, display = 'sites'))
caries_pct.mds.df$Groups <- caries_meta$group[match(row.names(caries_pct.mds.df), caries_meta$sample)]

# Colors
car_cols <- c("forestgreen", "firebrick4")

# Plot NMDS
ggplot(caries_pct.mds.df, aes(x=NMDS1, y=NMDS2)) + stat_ellipse(alpha=0.8, aes(color=Groups), show.legend = F, lwd=0.2) + 
  geom_point(alpha=0.9, aes(fill = Groups), size=3, color="black", pch=21, stroke=0.2) + scale_fill_manual(values=car_cols) + scale_color_manual(values=car_cols) + 
  labs(title = "Species Level: Beta Diversity Comparison", subtitle = "Relative Abundance", fill="Subject Groups") +
  annotate("text", x = max(caries_pct.mds.df$NMDS1)-0.2, y = max(caries_pct.mds.df$NMDS2), 
           label = paste("p <", caries_pct.perm$aov.tab$`Pr(>F)`[1], "(PERMANOVA)", sep=" ")) + theme_classic() +
  theme(plot.title = element_text(size=15, face="bold", hjust=0.5), plot.subtitle = element_text(size=10, hjust=0.5), 
        axis.title = element_text(size=10, face="bold"), axis.text = element_text(size=8, face="bold"), legend.position="bottom", axis.line = element_line(size = 0.3),
        legend.title = element_text(size=13, face="bold"), legend.text = element_text(size = 11))
ggsave(file = "figs/caries_pct_mds.pdf", width = 10, height = 6, units = "in")


# Presence/Absence
# Rarefy Counts for Pr/Ab
set.seed(12345); caries_rar <- as.data.frame(rrarefy(caries_counts, min(rowSums(caries_counts))))

# Convert to Pr/Ab
caries_prab <- data.frame((caries_rar > 0)*1, check.names = F)

# Compute NMDS
set.seed(12345); capture.output(caries_prab.mds <- metaMDS(caries_prab, trymax = 200, autotransform = F, wascores = F))

# Permanova
set.seed(12345); caries_prab.perm <- adonis(formula = caries_prab ~ caries_meta$group, strata = caries_meta$subject) # R2 0.19607, p < 0.001

# Making dataframe for plotting
caries_prab.mds.df <- data.frame(scores(caries_prab.mds, display = 'sites'))
caries_prab.mds.df$Groups <- caries_meta$group[match(row.names(caries_prab.mds.df), caries_meta$sample)]

# Colors
car_cols <- c("forestgreen", "firebrick4")

# Plot NMDS
ggplot(caries_prab.mds.df, aes(x=NMDS1, y=NMDS2)) + 
  #stat_ellipse(alpha=0.8, aes(color=Groups), show.legend = F, lwd=0.2) + 
  geom_point(alpha=0.9, aes(fill = Groups), size=3, color="black", pch=21, stroke=0.2) + scale_fill_manual(values=car_cols) + scale_color_manual(values=car_cols) + 
  labs(title = "Species Level: Beta Diversity Comparison", subtitle = "Presence-Absence", fill="Subject Groups") +
  annotate("text", x = max(caries_prab.mds.df$NMDS1)-0.2, y = max(caries_prab.mds.df$NMDS2), 
           label = paste("p <", caries_prab.perm$aov.tab$`Pr(>F)`[1], "(PERMANOVA)", sep=" ")) + theme_classic() +
  theme(plot.title = element_text(size=15, face="bold", hjust=0.5), plot.subtitle = element_text(size=10, hjust=0.5), 
        axis.title = element_text(size=10, face="bold"), axis.text = element_text(size=8, face="bold"), legend.position="bottom", axis.line = element_line(size = 0.3),
        legend.title = element_text(size=13, face="bold"), legend.text = element_text(size = 11))
ggsave(file = "figs/caries_prab_mds.pdf", width = 10, height = 6, units = "in")



# Methods to identify differentially abundance species


# I. Wilcoxon on each species + FDR ("BH")
# Wilcoxon test
caries_wilc <- sapply(colnames(caries_pct),function(s){
  #filters for analysis
  x <- caries_pct[,s]
  hf <- caries_meta$group %in% c("Intact Enamel","Dentin");
  #Health Plaque vs Caries Plaque
  f <- x[hf] ~ caries_meta$group[hf]
  # Wilcox test
  wilc.out <- wilcox.test(f, paired = T, exact = FALSE)
  return(list(
    pval=wilc.out$p.value
  ))
},simplify=FALSE)

# Tabulate and subset significant wilcox test results
caries_wilc_tab <- ldply(lapply(caries_wilc, function(x){data.frame(x)}))
length(which(caries_wilc_tab$pval < 0.05)) # 178

#Sig table
caries_wilc_tab$qval <- p.adjust(caries_wilc_tab$pval, method = "BH") # Apply false discovery correction
length(which(caries_wilc_tab$qval < 0.05)) # 154

# Filter pct table for sig fried
caries_pct_sig <- caries_pct[,colnames(caries_pct) %in% caries_wilc_tab$.id[caries_wilc_tab$qval < 0.05]]
dim(caries_pct_sig) # 68 x 154

# Melt relative abundance table
caries_pct_sig_m <- melt(as.matrix(caries_pct_sig))

# Assign colnames
colnames(caries_pct_sig_m) <- c("Sample","otu", "Levels")

# Add group field
caries_pct_sig_m$group <- caries_meta$group[match(caries_pct_sig_m$Sample, caries_meta$sample)]

# Add q value
caries_pct_sig_m$qvalue <- signif(caries_wilc_tab$qval[match(caries_pct_sig_m$otu, caries_wilc_tab$.id)],2)

# Levels to %
caries_pct_sig_m$Levels <- caries_pct_sig_m$Levels*100


# Separate increasing/decreasing

# Split by group
caries_caries_pct_split <- split(caries_pct_sig, droplevels(caries_meta$group))

# Sum for each group
caries_caries_pct_split_sum <- ldply(lapply(caries_caries_pct_split, colSums))
row.names(caries_caries_pct_split_sum) <- caries_caries_pct_split_sum$.id
caries_caries_pct_split_sum$.id <- NULL

# Transpose
caries_caries.wilc.t <- data.frame(t(caries_caries_pct_split_sum), check.names = F)
dim(caries_caries.wilc.t) # 72 x 2

# Increasing/Decreasing
caries_caries.wilc.t.disease <- caries_caries.wilc.t[caries_caries.wilc.t$Dentin > caries_caries.wilc.t$`Intact Enamel`,]
caries_caries.wilc.t.health <- caries_caries.wilc.t[caries_caries.wilc.t$`Intact Enamel` > caries_caries.wilc.t$Dentin,]

caries_caries.wilc.t.disease$diff <- caries_caries.wilc.t.disease$Dentin - caries_caries.wilc.t.disease$`Intact Enamel`

# Count
nrow(caries_caries.wilc.t.disease) # 23
nrow(caries_caries.wilc.t.health) # 131

# Rename
caries_wilc_species_markers <- caries_caries.wilc.t.disease








# Convert to Phyloseq
caries_ps <- phyloseq(otu_table(caries_counts, taxa_are_rows=FALSE), 
                      sample_data(caries_meta))

# Convert to DESeq2
caries_deseq <- phyloseq_to_deseq2(caries_ps, ~ group) # phyloseq data is converted to the relevant DESeqDataSet object, and design is included

# Compute DESeq2
caries_deseq_wald_par <- DESeq(caries_deseq, test="Wald", fitType="parametric")

# Summarize results
res <- results(caries_deseq_wald_par, cooksCutoff = FALSE)
alpha <- 0.01
sigtab <- res[which(res$padj < alpha), ]

# Make data frame
car_deseq_res <- cbind.data.frame(sigtab)
car_deseq_res$otu <- row.names(car_deseq_res)
car_deseq_res <- car_deseq_res[order(car_deseq_res$log2FoldChange, decreasing = T), ]

# Separate health and disease associated biomarkers
car_deseq_res_dis <- car_deseq_res[car_deseq_res$log2FoldChange > 0, ]
nrow(car_deseq_res_dis) # 29
car_deseq_res_health <- car_deseq_res[car_deseq_res$log2FoldChange < 0, ]
nrow(car_deseq_res_health) # 65

# Sort by fold change
car_deseq_res_dis$otu <- factor(car_deseq_res_dis$otu, levels = car_deseq_res_dis$otu[order(car_deseq_res_dis$log2FoldChange)])

# Plot Disease associated biomarkers
ggplot(car_deseq_res_dis, aes(x=otu, y=log2FoldChange, color=otu)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5), legend.key.size = unit(0.05, "cm"), 
        legend.text=element_text(size=5))
ggsave(file = "figs/caries_pct_deseq.pdf", width = 10, height = 6, units = "in")




## MetagenomeSeq

# Convert Phyloseq object to Metagenomeseq
caries_metseq <- phyloseq_to_metagenomeSeq(caries_ps) # phyloseq data is converted to the relevant MRexperiment-class object

# Cumulative sum scaling percentile selection
caries_metseq.p <- cumNormStat(caries_metseq)

# Normalize using the calculated pth quantile value
caries_metseq.norm = cumNorm(caries_metseq, p = caries_metseq.p)

# Retreive pheno data
caries_metseq.pheno <- pData(caries_metseq.norm) 

# Create model from pheno data
caries_metseq.mod <- model.matrix(~1 + group, data = caries_metseq.pheno) # creates a modelmatrix by expanding factors to a set of dummy variables 

# Computes differential abundance analysis using a zero-inflated log-normal mode 
caries_metseq.res <- fitFeatureModel(caries_metseq, caries_metseq.mod)

# Making table for top features
caries_metseq.res_top <- MRcoefs(caries_metseq.res, number=100, adjustMethod="BH") # Using Benjamini & Hochberg method for p-value adjustment for multiple comparisons

# Setting threshold as logFC > 2 and adjPvalue < 0.05
caries_metseq.res_sig <- caries_metseq.res_top[abs(caries_metseq.res_top$logFC) > 2 & caries_metseq.res_top$adjPvalues < 0.05, ] # q threshold can be < 0.1 ?

# Separate health and disease associated biomarkers
caries_metseq.res_sig_dis <- caries_metseq.res_sig[caries_metseq.res_sig$logFC > 0, ]
nrow(caries_metseq.res_sig_dis) # 30
caries_metseq.res_sig_health <- caries_metseq.res_sig[caries_metseq.res_sig$logFC < 0, ]
nrow(caries_metseq.res_sig_health) # 8
