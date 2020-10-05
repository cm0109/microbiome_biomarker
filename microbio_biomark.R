# Using Caries data for researching methods for microbiome-based biomarker discovery


library(reshape2)
library(vegan)
library(ggplot2)

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
ggplot(caries_prab.mds.df, aes(x=NMDS1, y=NMDS2)) + stat_ellipse(alpha=0.8, aes(color=Groups), show.legend = F, lwd=0.2) + 
  geom_point(alpha=0.9, aes(fill = Groups), size=3, color="black", pch=21, stroke=0.2) + scale_fill_manual(values=car_cols) + scale_color_manual(values=car_cols) + 
  labs(title = "Species Level: Beta Diversity Comparison", subtitle = "Presence-Absence", fill="Subject Groups") +
  annotate("text", x = max(caries_prab.mds.df$NMDS1)-0.2, y = max(caries_prab.mds.df$NMDS2), 
           label = paste("p <", caries_prab.perm$aov.tab$`Pr(>F)`[1], "(PERMANOVA)", sep=" ")) + theme_classic() +
  theme(plot.title = element_text(size=15, face="bold", hjust=0.5), plot.subtitle = element_text(size=10, hjust=0.5), 
        axis.title = element_text(size=10, face="bold"), axis.text = element_text(size=8, face="bold"), legend.position="bottom", axis.line = element_line(size = 0.3),
        legend.title = element_text(size=13, face="bold"), legend.text = element_text(size = 11))
ggsave(file = "figs/caries_prab_mds.pdf", width = 10, height = 6, units = "in")
