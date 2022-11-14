library(tidyverse)

df <- read.csv(file = './data/proteaselib.csv')
numInGroup <- 20

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
  if(nrow(splitDF[[x]]) < numInGroup)
    splitDF[[x]][nrow(splitDF[[x]])+1,] <- NA
}

# Plot out MWs vs. group ID as confirmation
df$groupID <- as.factor(df$groupID)
ggplot(df, aes(x = groupID, y = MW)) +
  geom_boxplot() +
  geom_point()

# Write out csv file for CoMa
write.csv(splitDF, './splitDF.csv', na="")