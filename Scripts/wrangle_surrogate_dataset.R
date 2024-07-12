# Wrangle surrogate calibration dataset
# Author: Mary Lofton
# Date: 12JUL24

# Purpose: rowbind all the individual model runs and also format observations for Parul

# load packages
library(tidyverse)
library(lubridate)

# rowbind model runs
runs <- list.files("./Output/", full.names = TRUE)
dat2 <- map_df(runs, read_csv) 

# plot model runs according to value of some parameter
ggplot(data = dat2, aes(x = datetime, y = prediction, group = R_growth_Nfixer, color = R_growth_Nfixer))+
  geom_line()+
  theme_classic()+
  scale_x_datetime(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))

# read in observation dataset
obs <- read_csv("./Data/EXOChla.csv") %>%
  filter(Depth == 1.6) %>%
  rename(observed = PHY_TCHLA) %>%
  select(-Depth)
write.csv(obs, file = "./Output/surrogate-observations_12JUL24.csv", row.names = FALSE)

# write final dataset to file
final <- dat2 %>%
  filter(datetime %in% obs$DateTime)
write.csv(final, file = "./Output/surrogate-GLM-AED-runs_12JUL24.csv",row.names = FALSE)

