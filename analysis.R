#  ██      ▄▄▄▄▄   █▀▄▀█    ▄▄▄▄▄       █ ▄▄  ██   █▄▄▄▄   ▄▄▄▄▄   ▄███▄   
#  █ █    █     ▀▄ █ █ █   █     ▀▄     █   █ █ █  █  ▄▀  █     ▀▄ █▀   ▀  
#  █▄▄█ ▄  ▀▀▀▀▄   █ ▄ █ ▄  ▀▀▀▀▄       █▀▀▀  █▄▄█ █▀▀▌ ▄  ▀▀▀▀▄   ██▄▄    
#  █  █  ▀▄▄▄▄▀    █   █  ▀▄▄▄▄▀        █     █  █ █  █  ▀▄▄▄▄▀    █▄   ▄▀ 
#  █               █                  █       █   █             ▀███▀   
#  █               ▀                    ▀     █   ▀                                                 
# Written by Michael Ronzetti NIH/NCATS 2023
library(tidyverse)
library(ggthemes)

# Read in raw df
raw.df <- read.csv(file = 'data/asms_data.csv')

# Convert fraction from chr to numeric
raw.df$Fraction <- gsub('%', '', raw.df$Fraction)
raw.df$Fraction <- as.numeric(raw.df$Fraction)
raw.df$Sample.Name <- as.factor(raw.df$Sample.Name)
raw.df$Sample.ID <- as.factor(raw.df$Sample.ID)

# Group by sample name, now using df.clean
df.clean <- raw.df %>%
  group_by(Sample.Name)

# Plot all data points
all.plot <-
  ggplot(df.clean, aes(x = Sample.Name, y = Fraction, fill = Sample.ID)) +
  geom_point(shape = 21, size = 4) +
  scale_fill_colorblind() +
  ylim(0, 100) +
  labs(title = 'ASMS Output',
       x = 'Sample Name',
       y = 'Fraction (%)') +
  theme(axis.text.x = element_text(
    angle = 60,
    vjust = 1,
    hjust = 1
  ))
all.plot

# Plot of elution and unbound
bind.plot <-
  ggplot(
    subset(
      df.clean,
      df.clean$Sample.ID == "Bound/Elution" |
        df.clean$Sample.ID == "Unbound"
    ),
    aes(x = Sample.Name, y = Fraction, fill = Sample.ID)
  ) +
  geom_point(shape = 21, size = 4) +
  scale_fill_colorblind() +
  ylim(0, 100) +
  labs(title = 'ASMS Output',
       x = 'Sample Name',
       y = 'Fraction (%)') +
  theme(axis.text.x = element_text(
    angle = 60,
    vjust = 1,
    hjust = 1
  ))
bind.plot
