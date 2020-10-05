# Using Caries data for researching methods for microbiome-based biomarker discovery

# Data source: 16S V1-V3 Sequencing of Severe Early Childhood Caries Affected Subjects
# Data subset for only Caries Cavity and Intact Enamel Samples

# Read species-level counts data
caries_counts <- readRDS("data/v13_car_ed_counts.rds") #Load counts data
caries_metadata <- readRDS("data/meta_c_ed.rds")  #Load metadata

dim(caries_counts) # 68 365

