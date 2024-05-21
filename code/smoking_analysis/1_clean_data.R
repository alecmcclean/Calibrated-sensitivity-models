####################################################
### Smoking data analysis cleaning
####################################################

########################
### Load data

smoking <- read.dta("../../data/C_2010_JOE-dataRS5K.dta")

### Dichotomize treatment
smoking$is_smoke <- ifelse(smoking$T > 0, 1, 0)
smoking <- select(smoking, -T) 

### Remove individuals with age zero
smoking %<>% filter(dfage > 0)

### Remove constant
smoking <- select(smoking, -const)
smoking_covariates <- 
  as.data.frame(
    model.matrix(~ -1 + ., data = smoking %>% select(-dbirwt, -is_smoke))
  ) %>%
  clean_names()

save(smoking, smoking_covariates, file = "../../intermediate/smoking/smoking_data.RData")

rm(list = ls(all = T))
gc()
