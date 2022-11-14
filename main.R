library(tidyverse)

df <- read.csv(file = './data/compoundlist.csv')
numInGroup <- 15

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

splitDF <- split(df, df$groupID, drop = true)

df$groupID <- as.factor(df$groupID)
ggplot(df, aes(x = groupID, y = MW)) +
  geom_boxplot() +
  geom_point()

write.csv(splitDF, './splitDF.csv')