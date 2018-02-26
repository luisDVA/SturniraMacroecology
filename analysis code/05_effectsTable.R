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
# drop SE and ICC, rename CI columns
Alleffects <- Alleffects %>% select(-SE,-ICC) %>% rename("CIlower"=lower,"CIupper"="upper")

# number formatting
Alleffects <-  Alleffects %>% mutate_if(is.numeric, funs(formatC(.,format="g",digits = 2))) 

# merge estimate and CIs to save space

Alleffects <-
  Alleffects %>% mutate(Estimate=paste0(estimate," ","(",CIlower," ","-"," ",CIupper,")")) %>% select(Model=body.measurement,Variable=variable,'Estimate (95CI)'=Estimate) 
Alleffects <-
Alleffects %>% 
    mutate(Variable = recode(Variable,
                        "specieshondurensis"= "species (hondurensis)",
                        "speciesparvidens"  = "species (parvidens)",
                        "sexM" = "sex (male)",
                        "specieshondurensis:sympatryboth" = "species (hondurensis): occurrence (sympatry)",
                        "speciesparvidens:sympatryboth" = "species (parvidens): occurrence (sympatry)"))


# for pasting in .doc
Alleffects%>% 
  write.csv(row.names = FALSE,quote = FALSE)

