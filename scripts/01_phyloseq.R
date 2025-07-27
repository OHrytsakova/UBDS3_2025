# # Install {phyloseq} package
# Installation instructions: https://joey711.github.io/phyloseq/install.html
# source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
#        local = TRUE)
# 
# # Check version
# packageVersion('phyloseq')

# Set the Environment ####
rm(list = ls()) # Reset R's brain

# Load necessary libraries
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(tibble)

# Load Ukraine admin boundaries from a file
load("./data/Ukraine_poly_adm1.Rdata")

# Define a default theme for ggplot graphics.
theme_set(theme_bw())

# Choose a marker
# marker <- "ITS2"
marker <- "SSU"


# Read data ####
if(marker == "ITS2") {
    # EcM only data (ITS2)
    physeq <- readRDS("./data/lotus2_ITS2/ecm_physeq.Rdata")
  } else {
    # AMF only (SSU)
    physeq <- readRDS("./data/lotus2_SSU_ASVs/amf_physeq.Rdata")
  }

# Adjust the whole phyloseq object ####
## Add Read_depth as a variable
sample_data(physeq)$Read_depth <- sample_sums(physeq)

#List read depth
get_variable(physeq, "Read_depth")

# Add Site_name to sample table
sample_data(physeq)$Site_name <- rownames(sample_data(physeq))


# Add reference sequences
# library(BiocManager)
# BiocManager::install("Biostrings")
library(Biostrings)

# Read the FASTA file
if (marker == "ITS2") {           # ITS2
  ref_seqs <- readDNAStringSet(
    "./data/lotus2_ITS2/OTU.fna", 
    format = "fasta")
} else {                          # SSU
  ref_seqs <- readDNAStringSet(
    "./data/lotus2_SSU_ASVs/OTU.fna", 
    format = "fasta")
}

# Remove extra descriptions if present
names(ref_seqs) <- gsub(" .*", "", names(ref_seqs))

# Convert reference sequences to a phyloseq object
refseq_physeq <- phyloseq(ref_seqs)

# Merge with the existing phyloseq object
physeq <- merge_phyloseq(physeq, refseq_physeq)

# Cleaning up
rm(ref_seqs, refseq_physeq)


# Extract main data blocks ####
# OTU/ASV abundance data
otu_table(physeq)       # View
otu_table <- physeq@otu_table@.Data

# taxonomy assignments
tax_table(physeq)       # View
tax_table <- physeq@tax_table@.Data

# sample metadata
sample_data(physeq)     # View
sample_data <- sample_data(physeq)

# the phylogenetic tree
phy_tree(physeq)        # View

# reference sequences
refseq(physeq)          # View

# Add reference sequences as variable to the taxonomy table
# Extract reference sequences as a character vector
ref_seqs_char <- as.character(refseq(physeq))

# Get the taxonomy table
tax_table_df <- as.data.frame(tax_table(physeq))

# Add reference sequences as a new column
tax_table_df$Reference_Sequence <- ref_seqs_char

# Convert back to a phyloseq-compatible taxonomy table
tax_table(physeq) <- as.matrix(tax_table_df)

# Check the updated taxonomy table
head(tax_table(physeq))




# Data cleaning ####
# Sample table
colnames(sample_data(physeq))

# Check columns for duplicates/errors
table(sample_data(physeq)$Site_name)

table(sample_data(physeq)$Country)

table(sample_data(physeq)$Land_use)

table(sample_data(physeq)$Ecosystem)

table(sample_data(physeq)$Vegetation)

sample_data(physeq)$Vegetation[
  sample_data(physeq)$Vegetation == "deciduoud_forest"
] <- "forest_deciduous"


# Filter Fungi only
fungi <- subset_taxa(physeq, Domain == "Fungi")
# There's an issue in the original data for SSU - 
# in the Domain variable only Eukaryota.
# Change it manually in filter above if you are using SSU marker.
# fungi <- physeq # for SSU only!

# Explore data ####
# See unique Phyla
unique(fungi@tax_table@.Data[,"Phylum"])

# Create table, number of features for each phyla
table(tax_table(fungi)[, "Phylum"], exclude = NULL)

# Basic visualizations ####

## Taxon histogram ####
plot_bar(fungi, fill = "Family")


p_phylum <- plot_bar(fungi, fill = "Phylum") +
  coord_flip()

# plot(p_phylum)

png(
  "./figures/phylum_sample_hist.png",
  width = 25, height = 40,
  units = "cm", res = 300
)
plot(p_phylum)
dev.off()


## Network Representation ####
# Tutorial: https://joey711.github.io/phyloseq/plot_network-examples.html

(p_net <- plot_net(fungi, 
                   # maxdist = 0.4, 
                   color="Ecosystem", 
                   # shape="Enterotype",
                   point_label = "Site_name"
))

png(
  "./figures/net.png",
  width = 25, height = 20,
  units = "cm", res = 300
)
plot(p_net)
dev.off()

# Create an igraph-based network based on the default distance 
# method, “Jaccard”, and a maximum distance between connected 
# nodes of 0.3
ig <- make_network(fungi, max.dist = 0.8)
(p_ig <- plot_network(ig, 
                      physeq, 
                      color="Ecosystem", 
                      # shape="Enterotype", 
                      line_weight=0.4))

png(
  "./figures/igraph_jac_0.8.png",
  width = 25, height = 20,
  units = "cm", res = 300
)
plot(p_ig)
dev.off()


## Tree ####
png(
  "./figures/tree.png",
  width = 35, height = 50,
  units = "cm", res = 300
)
plot_tree(fungi, 
          color="Phylum",
          # size="abundance"
)
dev.off()


# Diversity metrics ####
estimate_richness(fungi, measures = "Shannon")



# Subsets ####
### Selected OTUs ####

# Define the OTU name of interest
otu_id <- "OTU2362"

# Extract the reference sequence
otu_seq <- refseq(physeq)[otu_id]

# Print the sequence
otu_seq
as.character(otu_seq)

# Define the OTUs you want to keep
included_otu <- c("OTU2362", "OTU717", "OTU953")

# # Ensure only existing OTUs are selected
# included_otu <- included_otu[included_otu %in% taxa_names(physeq)]

# Subset phyloseq object
otu_subset <- prune_taxa(included_otu, physeq)

# Check result
otu_subset

tax_table(otu_subset)

### Selected taxa ####
# Filter entries with Order name
# Agaicales
agar <- subset_taxa(fungi, Order == "Agaricales")

agar_tax <- agar@tax_table@.Data

## Selected sites ####
# Subset samples by Site names
filtered_indices <- grep("^mazyk_", sample_data(physeq)$Site_name)
filtered_names <- sample_data(physeq)$Site_name[filtered_indices]
zyk <- prune_samples(filtered_names, physeq)

# Prune OTUs (taxa) that have zero counts across all samples
zyk <- prune_taxa(taxa_sums(zyk) > 0, zyk)

# Export to Excel
source("./functions/phyloseq_to_MDT_excel.R")
phyloseq_to_MDT_excel(zyk)

# End of the script ####
