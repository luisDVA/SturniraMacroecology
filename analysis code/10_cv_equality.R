# LMM CV 'variance' testing for forearm length with bioclimatic variables as covariates

# Requires
library(dplyr)
library(purrr)
library(nlme)
library(here)
# devtools::install_github("benmarwick/cvequality")
library(cvequality)
library(forcats)

# read the data
# measurements and environment
sturnfull <- read.csv(here("data","flwbiocl.csv"),stringsAsFactors = FALSE, encoding = "UTF-8") 
occPatterns <- read.csv(here("data","spOccPatterns.csv"),stringsAsFactors = FALSE,encoding = "UTF-8") %>% select(locality,sympatry)

# join to characterize localities
sturnfull <- left_join(sturnfull,occPatterns)
# subset for testing
FL_hond <- sturnfull %>% filter(species=="hondurensis")
FL_parv <- sturnfull %>% filter(species=="parvidens")

# test
FL_parv_cvTest <- 
  with(FL_parv, 
       mslr_test(nr = 1e4, 
                 FL,
                 sympatry))

FL_hond_cvTest <- 
  with(FL_hond, 
       mslr_test(nr = 1e4, 
                 FL,
                 sympatry))

#### Body mass
sturnBMsubset <- sturnfull %>% filter(!is.na(W))
# subset for testing
W_hond <- sturnBMsubset %>% filter(species=="hondurensis")
W_parv <- sturnBMsubset %>% filter(species=="parvidens")
# test
W_parv_cvTest <- 
    with(W_parv, 
         mslr_test(nr = 1e4, 
                   W,
                   sympatry))
W_hond_cvTest <- 
    with(W_hond, 
         mslr_test(nr = 1e4, 
                   W,
                   sympatry))

# head length
# read the data
batskulls <- read.csv(here("data","cranBiocl.csv"),stringsAsFactors = FALSE)
# occurrence in the localities
occPatterns <- read.csv(here("data","spOccPatterns.csv"),stringsAsFactors = FALSE,encoding = "UTF-8") %>% select(locality,sympatry)

# join to characterize localities
batskulls <- left_join(batskulls,occPatterns)
# subset for testing
CBL_hond <- batskulls %>% filter(species=="hondurensis")
CBL_parv <- batskulls %>% filter(species=="parvidens")
# test
CBL_parv_cvTest <- 
    with(CBL_parv, 
         mslr_test(nr = 1e4, 
                   CBL,
                   sympatry))
CBL_hond_cvTest <- 
    with(CBL_hond, 
         mslr_test(nr = 1e4, 
                   CBL,
                   sympatry))

# organize output
testsOut <- mget(ls(pattern = "cvTest$"))
cvtests <- 
map(testsOut,simplify2array) %>% as_tibble() %>% t() %>% as.data.frame %>% 
  tibble::rownames_to_column() %>% select("meassp"=1,"M-SLRT"=2,"p-value"=3) %>% 
  mutate(species=if_else(grepl("hond",meassp),"_S._ _hondurensis_","_S._ _parvidens_")) %>% 
  mutate(measurement=case_when(
    grepl("CBL",meassp)~"head length",
    grepl("FL",meassp)~"forearm length",
    TRUE~"body mass")) %>% select(measurement,species,everything()) %>% select(-meassp)

# reorder for reporting
cvtestsFin <- cvtests %>% arrange(forcats::fct_relevel(measurement,c("forearm length","body mass","head length")))

# write to disk
#fs::dir_create("variability_tests")
write.csv(cvtestsFin,file = here("variability_tests/","vartests.csv"),row.names = FALSE)















