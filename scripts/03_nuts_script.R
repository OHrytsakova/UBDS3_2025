# Set the Environment ####
rm(list = ls()) # Reset R's brain

# Load necessary libraries
library("tidyverse")

# Load from csv file
data <- read_csv("./nuts/nuts_data_raw.csv")

# Check data
head(data)
summary(data)

# Transform raw data to abundance table and wide format
data %>%
  group_by(region, sample) %>%
  summarise(n_ids = n(), .groups = "drop")

data_abund <- data %>% 
  group_by(region, sample, species) %>% 
  summarise(abundance = n(), .groups = "drop")

head(data_abund)

# Convert to wide format
data_wide <- data_abund %>% 
  pivot_wider(names_from = species,
              values_from = abundance,
              values_fill = 0) %>% 
  unite("sample_id", region, sample, remove = T) %>% 
  column_to_rownames("sample_id")

write.csv(data_wide, "./nuts/data_wide_spe.csv")


# Let's explore data!
# Look at species richness

data_regions <- data_abund %>% 
  group_by(region, species) %>% 
  summarise(abundance = sum(abundance))

library(vegan) # vegan = "vegetation analysis"

# How many species per region?
SR_region <- data_regions %>% 
  group_by(region) %>% 
  summarise(SR = specnumber(abundance))

SR_region

# How many species per sample?
SR_sample <- data_abund %>% 
  group_by(region, sample) %>% 
  summarise(SR = specnumber(abundance))
  
SR_sample

# Mean species richness per region
mean_SR_sample <- SR_sample %>%
  summarise(mSR = mean(SR))
  
mean_SR_sample

# Add mean species richness as a variable
SR_all <- SR_region %>% 
  mutate(mSR_sample = mean_SR_sample$mSR)

# Create simple plot
library(ggplot2)

ggplot(SR_sample, aes(x = region, y = SR, colour = region)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(
    title = "Species richness",
    y = "species richness"
  ) +
  theme_bw()

?specnumber


# Abundance curves
dat_landscapes_sort <- data_regions |>  
  arrange(region, desc(abundance))

dat_landscapes_sort <- dat_landscapes_sort |> 
  group_by(region) |> 
  mutate(rank = row_number())

ggplot(dat_landscapes_sort, aes(x = rank, y = abundance, color = region, shape = region)) +
  geom_line() +
  geom_point(size = 3) +
  scale_y_log10() +
  ylab("abundance") +
  ggtitle("Rank-abundance-curves all regions") +
  theme_bw()

# End of the script ####

