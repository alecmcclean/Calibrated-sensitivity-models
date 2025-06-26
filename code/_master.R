####################################################
### Data analyses
####################################################

set.seed(20231203)

###########################################
### Load npcausal for ATE estimation

library(npcausal) ### see github/ehkennedy for installation instructions

###########################################
### Load sensemakr for Darfur data

library(sensemakr)

#########################################
### Package for loading smoking data 

library(foreign)

############################################################
### Load packages for data manipulation and plots

library(tidyverse)
library(magrittr)
library(janitor)
library(ggthemes) 
library(latex2exp)
library(scales)
library(foreach)
library(doParallel)


############################
### Darfur data analysis

source("./darfur_analysis/1_clean_data.R")
source("./darfur_analysis/2_fx_diff.R")
source("./darfur_analysis/3_bdd_odds_ratio.R")
source("./darfur_analysis/4_bdd_odds_ratio_M.R")

source("./1_darfur_figures.R")

############################
### Smoking data analysis

source("./smoking_analysis/1_clean_data.R")
source("./smoking_analysis/2_fx_diff.R")

source("./2_smoking_figures.R")
