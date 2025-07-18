# Exporting a PhyloSeq object for the MDT ####
# https://docs.gbif-uat.org/mdt-user-guide/en/#from-phyloseq-to-excel


# Adjust column names to Darwin Core terms

# colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# From PhyloSeq to Excel

# Load necessary libraries
library(phyloseq)
library(openxlsx)

# Function to save Excel workbook
phyloseq_to_MDT_excel <- function(physeq, filename = "Phyloseq_Tables.xlsx") {
  # Load neccessary library
  library(openxlsx)
  # Create a new Excel workbook
  wb <- createWorkbook()
  
  # Save OTU Table
  otu_table_df <- data.frame(otu_table(physeq))
  addWorksheet(wb, "OTU_table")
  writeData(wb, "OTU_table", cbind(id = rownames(otu_table_df), otu_table_df))
  
  # Save Sample Data
  sample_data_df <- data.frame(sample_data(physeq))
  addWorksheet(wb, "Samples")
  writeData(wb, "Samples", cbind(id = rownames(sample_data_df), sample_data_df))
  
  # Save Taxonomy Table
  tax_table_df <- data.frame(tax_table(physeq))
  addWorksheet(wb, "Taxonomy")
  writeData(wb, "Taxonomy", cbind(id = rownames(tax_table_df), tax_table_df))
  
  # Save the workbook
  saveWorkbook(wb, filename, overwrite = TRUE)
}

# Example usage
# Replace 'your_phyloseq_object' with your actual phyloseq object
# phyloseq_to_MDT_excel(your_phyloseq_object)