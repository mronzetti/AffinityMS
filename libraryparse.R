#
#  ____            __  __ _      _
# |  _ \  __ _ ___|  \/  (_) ___(_)_ __   __ _
# | | | |/ _` / __| |\/| | |/ __| | '_ \ / _` |
# | |_| | (_| \__ \ |  | | | (__| | | | | (_| |
# |____/ \__,_|___/_|  |_|_|\___|_|_| |_|\__,_|
#
# Written by: Michael Ronzetti NIH/NCATS 2022
#
# Sorts out compound list to maximize difference in molecular weights of a screening run.
library(tidyverse)
library(ggthemes)
library(rcdk)
library(readxl)

# Parameters for library parsing
# 
# numInGroup = number of compounds per well
# echoDispVol = dispense volume on Echo (nL)
numInGroup <- 10
echoDispVol <- 8
echoDispType <- '1536LDV_DMSO'

# Import a csv file with sample id and smiles
raw.df <- read_xlsx(path = './data/protease_compoundlist.xlsx')

# Comment previous line and uncomment this to specify groups by plate size.
#numSmolecules <- nrow(df)
#plateWells <- 30
#numInGroup <- round(numSmolecules / plateWells, digits = 0)

# Setup df to store values
result.df <-
  data.frame(sample.id = character(),
             formula = character(),
             mass = numeric())

# Calculate the molecular formula and molecular weight with rcdk package
for (x in 1:nrow(raw.df)) {
  smiles_string <- raw.df$smiles[x]
  
  # Check if smiles value is NULL, if so, assign NA to columns
  if (is.null(smiles_string) | smiles_string == 'null') {
    cat("Warning: Missing SMILES data at row", x, "\n")
    formula <- NA
    mass <- 1
    sample.df <- data.frame(sample.id = raw.df$sample.id[x],
                            formula = formula,
                            mass = mass)
  } else {
    sp <- get.smiles.parser()
    molecule.test <- parse.smiles(smiles = raw.df$smiles[x], kekulise = FALSE)[[1]]
    convert.implicit.to.explicit(molecule.test)
    formula <- get.mol2formula(molecule.test, charge = 0)
    message(paste(
      raw.df$sample.id[x],
      ' at row: ',
      x,
      '; formula: ',
      formula@string,
      '; mass: ',
      formula@mass,
      '\n',
      sep = ''
    ))
    sample.df <- data.frame(sample.id = raw.df$sample.id[x],
                            formula = formula@string,
                            mass = formula@mass)
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

# Export a master file with addresses
master.df <- inner_join(raw.df, df)
master.df <- master.df %>% arrange(master.df$groupID)

# create a vector of well IDs for a 96-well plate
wellIDs <- paste0(LETTERS[1:8], rep(1:12, each = 8))

# assign a unique well ID to each unique groupID
unique_groups <- unique(master.df$groupID)
wellIDs_group <- wellIDs[1:length(unique_groups)]
names(wellIDs_group) <- unique_groups
master.df$wellID <- wellIDs_group[match(master.df$groupID, names(wellIDs_group))]

# Export the master df
write.csv(x = master.df, file = './output/masterList.csv', row.names = FALSE)

# Create ECHO file for dispense
echo.df <- master.df %>%
  select(sample.id, source.well, wellID) %>%
  mutate('Transfer Volume' = echoDispVol) %>%
  rename('Source Well' = source.well) %>%
  rename('Destination Well' = wellID) %>%
  mutate('Source Plate Type' = echoDispType)

# Remove leading zeros from the Source Well column
echo.df$`Source Well` <- gsub("(?<=[A-Za-z])0+(?=[0-9])", "", echo.df$`Source Well`, perl = TRUE)

# Export the ECHO df
write.csv(x = echo.df, file = './output/echodispense.csv', row.names = FALSE)

# Sp  lit the df according to group ID
splitDF <- split(df, df$groupID, drop = false)

# Add NA to last row to equalize row numbers
for (x in 1:length(splitDF)) {
  if (nrow(splitDF[[x]]) < numInGroup)
    splitDF[[x]][nrow(splitDF[[x]]) + 1,] <- NA
}

# Plot out MWs vs. group ID as confirmation
df$groupID <- as.factor(df$groupID)
ggplot(df, aes(x = groupID, y = mass)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, shape = 21, aes(fill = groupID)) +
  theme_clean() +
  theme(legend.position = 'none', axis.text = element_text(size = 4)) +
  labs(title = 'AS-MS Compound Splitting')

ggsave(
  './output/compound_split.png',
  dpi = 600,
  scale = 2.5,
  plot = last_plot()
)

# Write out csv file for CoMa
write.csv(splitDF, './output/splitDF.csv', na = "")