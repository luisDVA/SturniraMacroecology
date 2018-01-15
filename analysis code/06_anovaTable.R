# reporting LMM anova
library(dplyr)
library(purrr)
library(here)
library(forcats)

# read the exported effects for reporting
Allaov <- 
  list.files(here("lmmEffects"), pattern = "^aov", full.names = TRUE) %>% 
  map_df(read.csv) 

# relevel and arrange
Allaov$measurement <-  fct_relevel(Allaov$measurement,"forearm length","body mass","skull isosize","head length")
Allaov <- Allaov %>% arrange(measurement)

# rename column and drop DF
Allaov <- Allaov %>% rename(model=measurement) %>% select(-df)
# for copying into document
Allaov %>% mutate_if(is.numeric, funs(formatC(.,format="g",digits = 3))) %>% 
  write.csv(quote = FALSE,row.names = FALSE)
