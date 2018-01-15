# reporting LMM estimates
library(dplyr)
library(purrr)
library(here)
library(forcats)

# read the exported effects for reporting
Alleffects <- 
  list.files(here("lmmEffects"), pattern = "output.csv$", full.names = TRUE) %>% 
  map_df(read.csv) %>% select(-est.)

# relevel and arrange
Alleffects$body.measurement <-  fct_relevel(Alleffects$body.measurement,"forearm length","body mass","skull isosize","head length")
Alleffects <- Alleffects %>% arrange(body.measurement)

# rename column
Alleffects <- Alleffects %>% rename(variable="term",SE=std.error)
# round ICC
Alleffects$ICC <- round(Alleffects$ICC,3)

# for reporting
# drop SE, rename CI columns
Alleffects <- Alleffects %>% select(-SE) %>% rename("CIlower"=lower,"CIupper"="upper")
# pasteICC to model name to save space
Alleffects <- Alleffects %>% mutate(model=paste0(body.measurement," ","(",ICC,")")) %>% 
  select(model,everything()) %>% select(-body.measurement) %>% select(-ICC)

# for copying into document
Alleffects%>% mutate_if(is.numeric, funs(formatC(.,format="g",digits = 3))) %>% 
    write.csv(quote = FALSE,row.names = FALSE)

