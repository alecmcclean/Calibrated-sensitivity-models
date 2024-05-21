####################################################
### Darfur data analysis cleaning
####################################################

###################################
### Clean data

### Load Darfur data
data("darfur")

### Bucket all villages with less than 10 individuals in "Other"
darfur$village <- as.character(darfur$village)
darfur <- darfur %>%
  mutate(village = as.character(village)) %>%
  group_by(village) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  mutate(village = ifelse(count > 10, village, "Other"))

### Create model matrix
darfur_covariates <- 
  as.data.frame(
    model.matrix(~ -1 + female + age + farmer_dar + herder_dar + pastvoted + 
                   hhsize_darfur + village, data = darfur)
  ) %>%
  clean_names()

darfur$peacefactor <- darfur$peacefactor[,1]

save(darfur, darfur_covariates, file = "../../intermediate/darfur/darfur_data.RData")

rm(list = ls(all = T))
gc()
