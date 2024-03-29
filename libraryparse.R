#
#  ____            __  __ _      _
# |  _ \  __ _ ___|  \/  (_) ___(_)_ __   __ _
# | | | |/ _` / __| |\/| | |/ __| | '_ \ / _` |
# | |_| | (_| \__ \ |  | | | (__| | | | | (_| |
# |____/ \__,_|___/_|  |_|_|\___|_|_| |_|\__,_|
#
# Written by: Michael Ronzetti NIH/NCATS 2022
#
#   To Do
#
#   Need to make wide form compatible for Nate's analysis (wait for Sam meeting)
#
# Sorts out compound list to maximize difference in molecular weights of a screening run, generates
#   exportable ECHO file, control compound dispense file

# List of packages to check/install
packages <- c("tidyverse", "ggthemes", "rcdk", "readxl")

# Function to check if a package is installed
is_package_installed <- function(package_name) {
  is.element(package_name, installed.packages()[, "Package"])
}

# Iterate over the packages and install if not installed
for (package in packages) {
  if (!is_package_installed(package)) {
    install.packages(package)
  }
}
library(tidyverse)
library(ggthemes)
library(rcdk)
library(readxl)

#   Parameters for library parsing
#
#   numInGroup = number of compounds per well
#   echoDispVol = dispense volume on Echo (nL)
#   echoDispType = Type of dispensing protocol
#   numControlWells = # of wells to leave empty for controls
#   numTotalWells = Plate format for experiment
#
numInGroup <- 10
echoDispVol <- 20
echoDispType <- '1536LDV_DMSO'
numControlWells <- 5
numTotalWells <- 96

# Import csv file with sample id and smiles
# Requires:
#   plate.id, source.well,sample.id, smiles
raw.df <- read_xlsx(path = './data/AID_complete.xlsx')

# Setup result.df to store compound values
result.df <-
  data.frame(sample.id = character(),
             formula = character(),
             mass = numeric())

# Calculate the molecular formula and molecular weight with rcdk package for each
for (x in 1:nrow(raw.df)) {
  smiles_string <- raw.df$smiles[x]
  # Check if smiles value is NULL, if so, assign NA to columns
  if (is.null(smiles_string) | smiles_string == 'null') {
    cat("Warning: Missing SMILES data at row", x, "\n")
    formula <- NA
    mass <- 1
    sample.df <- data.frame(
      sample.id = raw.df$sample.id[x],
      formula = formula,
      mass = mass
    )
  } else {
    sp <- get.smiles.parser()
    molecule.test <-
      parse.smiles(smiles = raw.df$smiles[x], kekulise = FALSE)[[1]]
    convert.implicit.to.explicit(molecule.test)
    formula <- get.mol2formula(molecule.test, charge = 0)
    sample.df <- data.frame(
      sample.id = raw.df$sample.id[x],
      formula = formula@string,
      mass = formula@mass
    )
  }
  result.df <- rbind(result.df, sample.df)
}

# Sort the entire list by MW in ascending order
df <- arrange(result.df, mass)

# Find the # of groups needed for the sorting
numGroups <- as.integer(ceiling(nrow(df) / numInGroup))
groupCount <- 1

# Initialize groupID column
df$groupID <- as.integer(0)

# Attach the groupID num to each row
for (x in 1:nrow(df)) {
  df$groupID[x] <- groupCount
  groupCount = groupCount + 1
  if (groupCount == numGroups + 1)
    groupCount = 1
}

# Calculate the maximum number of groups that can be placed on a single plate
groupsPerPlate <- numTotalWells - numControlWells

# Create 'destination.plate' column
df$destination.plate <- ceiling(df$groupID / groupsPerPlate)
df$destination.plate <-
  paste0("Plate", sprintf("%03d", df$destination.plate))

# Construct a master file with addresses
master.df <- inner_join(raw.df, df)
master.df <-
  master.df %>% arrange(master.df$groupID, master.df$destination.plate)

# Create a vector of well IDs
fullplatewellIDs <- paste0(LETTERS[1:8], rep(1:12, each = 8))
wellIDs <- fullplatewellIDs[1:groupsPerPlate]

# Assign a unique well ID to each unique groupID within each plate
unique_plates <- unique(df$destination.plate)
wellIDs_all <- character()
for (plate in unique_plates) {
  df_plate <- df[df$destination.plate == plate, ]
  unique_groups_plate <- unique(df_plate$groupID)
  wellIDs_group <- wellIDs[1:length(unique_groups_plate)]
  names(wellIDs_group) <-
    paste(plate, unique_groups_plate, sep = "_")
  wellIDs_all <- c(wellIDs_all, wellIDs_group)
}

df$wellID <-
  wellIDs_all[paste(df$destination.plate, df$groupID, sep = "_")]

# Construct a master file with addresses
master.df <- inner_join(raw.df, df)
master.df <-
  master.df %>% arrange(master.df$groupID, master.df$destination.plate)

# Remove leading zeros from the Source Well column
master.df$source.well <-
  gsub("(?<=[A-Za-z])0+(?=[0-9])", "", master.df$source.well, perl = TRUE)

# Export the master df
write.csv(x = master.df,
          file = './output/masterList.csv',
          row.names = FALSE)

# Split the data frame by 'wellID' into a list of data frames
filtered.df <- master.df %>%
  select(sample.id, formula, wellID)
split_df <- split(filtered.df, filtered.df$wellID)

# Function to reshape the sub-data frames
reshape_sub_df <- function(df) {
  well_id <- unique(df$wellID)
  df <- df %>%
    mutate(row_num = row_number()) %>%
    pivot_wider(
      names_from = row_num,
      values_from = c(sample.id, formula),
      names_sep = " "
    ) %>%
    select(-wellID) %>%
    add_column(wellID = well_id, .before = 1)
  return(df)
}

# Apply the function to each sub-data frame in the list
reshaped_sub_dfs <- lapply(split_df, reshape_sub_df)

# Merge the reshaped data frames back together
result_df <- bind_rows(reshaped_sub_dfs)

# Generate the desired column order
num_pairs <- 10
col_order <- c("wellID",
               unlist(lapply(1:num_pairs, function(i)
                 c(
                   paste("sample.id", i), paste("formula", i)
                 ))))

# Arrange columns in the desired order
result_df <- result_df %>% select(all_of(col_order))

# Print the result
print(result_df)

# Export compound list
write.csv(x = result_df, file = './output/wideDF.csv', quote = FALSE)

# Create ECHO file for dispense
echo.df <- master.df %>%
  select(sample.id, source.well, wellID, plate.id, destination.plate) %>%
  mutate('Transfer Volume' = as.integer(echoDispVol)) %>%
  rename('Source Well' = source.well) %>%
  rename('Destination Well' = wellID) %>%
  rename('Destination Plate Name' = destination.plate) %>%
  mutate('Source Plate Type' = as.character(echoDispType)) %>%
  mutate('Source Plate Name' = as.character(plate.id)) %>%
  select(-plate.id)

# Remove leading zeros from the Source Well column
echo.df$`Source Well` <-
  gsub("(?<=[A-Za-z])0+(?=[0-9])", "", echo.df$`Source Well`, perl = TRUE)

# Export the ECHO df
write.csv(
  x = echo.df,
  file = './output/echodispense.csv',
  row.names = FALSE,
  quote = FALSE
)

# Plot out MWs vs. group ID as confirmation
df$groupID <- as.factor(df$groupID)
ggplot(df, aes(x = groupID, y = mass)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, shape = 21, aes(fill = groupID)) +
  theme_clean() +
  theme(legend.position = 'none', axis.text = element_text(size = 8)) +
  labs(title = 'AS-MS Compound Splitting')

# ggsave(
#   './output/compound_split.png',
#   dpi = 600,
#   scale = 2.5,
#   plot = last_plot()
# )

# # # # # # # # # #
# Count the number of compounds in each well
compounds_per_well <- df %>%
  group_by(wellID, destination.plate) %>%
  summarise(numCompounds = n(), .groups = "keep")

# Identify wells with fewer than numInGroup compounds, filter, and find DMSO backfill amount
compounds_per_well$needControl <-
  compounds_per_well$numCompounds < numInGroup
backfill_wells <-
  filter(compounds_per_well, needControl == TRUE) %>%
  mutate(backfillAmount = (numInGroup - numCompounds) * echoDispVol) %>%
  select(!c(needControl, numCompounds))

# Identify which wells are unused and for control compound
control_wells <- fullplatewellIDs[!fullplatewellIDs %in% wellIDs]
message(control_wells)

# Export control dispense dataframe
write.csv(
  x = backfill_wells,
  file = './output/control_dispense.csv',
  row.names = FALSE,
  quote = FALSE
)
