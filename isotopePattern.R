# Reads in sample SMILES and returns isotope pattern estimates

library(rcdk)
library(tidyverse)
library(openxlsx)

# Read in csv with sample.id in first column and structure (SMILES) in second column.
raw.df <- read.csv('./data/compound_info.csv')

# Create an empty data frame to store the results
result_df <-
  data.frame(sample.id = character(),
             mass = numeric(),
             abund = numeric())

# Read in each individual SMILES and find isotope patterns
for (x in 1:nrow(raw.df)) {
  sp <- get.smiles.parser()
  molecule.test <- parse.smiles(smiles = raw.df$structure[x])[[1]]
  convert.implicit.to.explicit(molecule.test)
  formula <- get.mol2formula(molecule.test, charge = 0)
  isotopes <- get.isotopes.pattern(formula, minAbund = 0.05)
  print(isotopes)
  
  # Setup png export and plot
  filename <-
    paste('./output/isotope_plots/', raw.df$sample.id[x], '_isotopes.png', sep = '')
  png(
    filename = filename,
    width = 1600,
    height = 900,
    res = 300,
    pointsize = 6
  )
  plot(
    isotopes,
    type = 'h',
    xlab = 'mass',
    ylab = 'intensity',
    main = raw.df$sample.id[x],
    lwd = 3,
    col = 'purple'
  )
  dev.off()
  
  # Create a data frame for the current sample with sample.id, mass, and abund
  sample_df <- data.frame(sample.id = raw.df$sample.id[x],
                          mass = isotopes[, 1],
                          abund = isotopes[, 2])
  
  # Append the sample_df to the result_df
  result_df <- rbind(result_df, sample_df)
}

write.xlsx(result_df, file = './output/isotopeDF.xlsx')