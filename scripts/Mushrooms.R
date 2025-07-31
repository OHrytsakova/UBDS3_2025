library("tidyverse")
library(iNEXT)

data_t<- read_csv("./data/sh_abundance.csv")
data_t <- subset(data_t, select = -c(1))
data_t <- column_to_rownames(data_t, var = "scientificName")

data_wide <- t(data_t) %>% 
  as.matrix()



# Incidence approach - sample coverage
data_t_pa <- ifelse(data_t>0, 1, 0)

m2 <- list(data_t_pa)
names(m2) <- "IncidenceData - Fungi"
out.pf <- iNEXT(m2, q = c(0,1,2), datatype = "incidence_raw")

ggiNEXT(out.pf, type = 2)

chao1 <- iNEXT::ChaoRichness(m2, datatype = "incidence_raw")

chao1


# Incidence approach - species richness
# Convert counts to presence/absence
data_inc <- data_wide > 0

# Count in how many samples each species occurs
species_incidence_freq <- colSums(data_inc)

# Totalnumber of sampling units
n_units <- nrow(data_wide)

# Create the input for iNEXT
incidence_input <- list(c(n_units, species_incidence_freq))
names(incidence_input) <- "IncidenceData - Fungi"

out_inc <- iNEXT(incidence_input, q = 0, datatype = "incidence_freq")

ggiNEXT(out_inc, type = 1)
