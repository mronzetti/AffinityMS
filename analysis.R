#  ██      ▄▄▄▄▄   █▀▄▀█    ▄▄▄▄▄       █ ▄▄  ██   █▄▄▄▄   ▄▄▄▄▄   ▄███▄
#  █ █    █     ▀▄ █ █ █   █     ▀▄     █   █ █ █  █  ▄▀  █     ▀▄ █▀   ▀
#  █▄▄█  ▄ ▀▀▀▀▄   █ ▄ █ ▄  ▀▀▀▀▄       █▀▀▀  █▄▄█ █▀▀▌  ▄▀▀▀▀▄   ██▄▄
#  █  █  ▀▄▄▄▄▀    █   █  ▀▄▄▄▄▀        █     █  █ █  █  ▀▄▄▄▄▀  █▄▄▀
#  █               █                    █     █    █             ▀███▀
#  █               ▀                    ▀     █    ▀
# Written by Michael Ronzetti NIH/NCATS 2023
# To Do
#  2x2 graph of ASMS output by fraction area and percentage

# Need to Do: Raw Area, change efficiency to % in elution.

library(tidyverse)
library(ggthemes)
library(cowplot)

# Read in raw df
raw.df <- read.csv(file = 'data/asms_data.csv')

# Convert fraction from chr to numeric
raw.df$Fraction <- gsub('%', '', raw.df$Fraction)
raw.df$Fraction <- as.numeric(raw.df$Fraction)
raw.df$Sample.Name <- as.character(raw.df$Sample.Name)
raw.df$Sample.ID <- as.character(raw.df$Sample.ID)

# Group by sample name, now using df.clean
df.clean <- raw.df %>%
  group_by(Sample.Name)

# Calculate binding efficiency
# Defined as (Elution fraction)/(Rest of fractions)
df.bindingEff <-
  tibble(sample.Name = unique(df.clean$Sample.Name),
         bindingEff = as.numeric(0))
for (x in 1:nrow(df.bindingEff)) {
  elution <-
    df.clean$Fraction[df.clean$Sample.Name == df.bindingEff$sample.Name[x] &
                            df.clean$Sample.ID == "Bound/Elution"]
  wash1 <-
    df.clean$Fraction[df.clean$Sample.Name == df.bindingEff$sample.Name[x] &
                            df.clean$Sample.ID == "Wash1"]
  wash2 <-
    df.clean$Fraction[df.clean$Sample.Name == df.bindingEff$sample.Name[x] &
                            df.clean$Sample.ID == "Wash2"]
  unbound <-
    df.clean$Fraction[df.clean$Sample.Name == df.bindingEff$sample.Name[x] &
                            df.clean$Sample.ID == "Unbound"]
  df.bindingEff$bindingEff[x] <- elution / (wash1 + wash2 + unbound)
}

# Prepare a df that is in wide format for Prism input
# Working with two grouping variables.
df.prism <-
  pivot_wider(
    data = df.clean,
    names_from = Sample.Name,
    values_from = Analyte.Area,
    id_cols = Sample.ID
  )

# Widen the binding efficiency tables
df.bindingEff.Temp <-
  pivot_wider(data = df.bindingEff,
              names_from = sample.Name,
              values_from = bindingEff)

# Attach the Binding Efficiency row and rename the NA
df.prism <- df.prism %>%
  add_row(df.bindingEff.Temp)
df.prism$Sample.ID <-
  df.prism$Sample.ID %>% replace_na("BindingEff")

# Export the Prism df as a csv
write.csv(x = df.prism,
          file = './output/ASMSPrism.csv',
          row.names = FALSE)

# Graph ASMS output
# Plot all data points
all.plot <-
  ggplot(df.clean, aes(x = Sample.Name, y = Fraction, fill = Sample.ID)) +
  geom_point(shape = 21, size = 2) +
  scale_fill_colorblind() +
  ylim(0, 100) +
  labs(title = 'ASMS Output',
       subtitle = 'All Fractions',
       y = 'Fraction (%)') +
  theme_clean() +
  theme(axis.text.x = element_text(
    angle = 60,
    vjust = 1,
    hjust = 1
  ),
  axis.title.x = element_blank())
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
  geom_point(shape = 21, size = 2) +
  scale_fill_colorblind() +
  ylim(0, 100) +
  labs(title = 'ASMS Output',
       subtitle = 'Binding Fractions',
       y = 'Fraction (%)') +
  theme_clean() +
  theme(axis.text.x = element_text(
    angle = 60,
    vjust = 1,
    hjust = 1
  ),
  axis.title.x = element_blank())
bind.plot

# Plot the ASMS Binding Efficiency
eff.plot <-
  ggplot(df.bindingEff, aes(x = sample.Name, y = bindingEff)) + geom_point(shape = 21,
                                                                           fill = 'blue',
                                                                           size = 4) +
  labs(title = 'ASMS Output',
       subtitle = 'Binding Efficiency') +
  theme_clean() +
  theme(axis.text.x = element_text(
    angle = 60,
    vjust = 1,
    hjust = 1
  ),
  axis.title.x = element_blank()) +
  labs(title = 'ASMS Output',
       subtitle = 'Binding Efficiency',
       y = 'Efficiency')
eff.plot

# Cowplot these together
fractionsLeft <- plot_grid(all.plot, bind.plot, ncol = 1)
fractionsPlot <- plot_grid(fractionsLeft, eff.plot, nrow = 1)
fractionsPlot

# Save fraction cowplot to output folder
ggsave(filename = './output/fractionsPlot.png', plot = fractionsPlot, dpi = 'retina', scale = 2)
