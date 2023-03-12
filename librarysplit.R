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

# Import a csv file with sample id and molecular weight column (as 'MW')
df <- read.csv(file = './data/proteaselib.csv')

# Specify the number of smolecules per group here.
#numInGroup <- 20

# Comment previous line and uncomment this to specify groups by plate size.
numSmolecules <- nrow(df)
plateWells <- 30
numInGroup <- round(numSmolecules / plateWells, digits = 0)

# Sort the entire list by MW in ascending order
df <- arrange(df, MW)

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

# Split the df according to group ID
splitDF <- split(df, df$groupID, drop = false)

# Add NA to last row to equalize row numbers
for (x in 1:length(splitDF)) {
  if (nrow(splitDF[[x]]) < numInGroup)
    splitDF[[x]][nrow(splitDF[[x]]) + 1, ] <- NA
}

# Plot out MWs vs. group ID as confirmation
df$groupID <- as.factor(df$groupID)
ggplot(df, aes(x = groupID, y = MW)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, shape = 21, aes(fill = groupID)) +
  theme_clean() +
  theme(legend.position = 'none', axis.text = element_text(size = 4)) +
  labs(title = 'AS-MS Compound Splitting')

ggsave(
  'compound_split.png',
  dpi = 600,
  scale = 2.5,
  plot = last_plot()
)

# Write out csv file for CoMa
write.csv(splitDF, './splitDF.csv', na = "")
