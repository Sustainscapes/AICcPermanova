12/01, 2023

-   [1 AICcPermanova](#1-aiccpermanova)
-   [2 Vegetation models](#2-vegetation-models)
    -   [2.1 Presence-absence data](#21-presence-absence-data)
        -   [2.1.1 Generation of all possible
            models](#211-generation-of-all-possible-models)
        -   [2.1.2 Model fitting](#212-model-fitting)
    -   [2.2 Abundance data](#22-abundance-data)
        -   [2.2.1 Generation of all possible
            models](#221-generation-of-all-possible-models)
        -   [2.2.2 Model fitting](#222-model-fitting)
-   [3 Baterial models](#3-baterial-models)
    -   [3.1 Abundance data](#31-abundance-data)
        -   [3.1.1 Generation of all possible
            models](#311-generation-of-all-possible-models)
        -   [3.1.2 Model fitting](#312-model-fitting)
    -   [3.2 Presence absence data](#32-presence-absence-data)
        -   [3.2.1 Generation of all possible
            models](#321-generation-of-all-possible-models)
        -   [3.2.2 Model fitting](#322-model-fitting)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# 1 AICcPermanova

<!-- badges: start -->
<!-- badges: end -->

The goal of AICcPermanova is to evaluate the best models for plant
communities and bacterial communities in Denmark in order to do that we
require the following packages

``` r
library(vegan)
library(ampvis2)
library(dplyr)
library(readxl)
library(parallel)
library(foreach)
library(doParallel)
library(tidyr)
library(car)
```

# 2 Vegetation models

## 2.1 Presence-absence data

### 2.1.1 Generation of all possible models

After that we read in the datasets for environmental layers and generate
all possible models to fit, in this case we will limit ourselves to only
using at most one variable per ten observations, in this case that means
up to 5 variables per model. The code for generating all possible models
can be expanded bellow

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Model generator

</summary>

``` r
METADATAS <- list.files(pattern = "PERMANOVA_VEGETATION_", full.names = T)

AllForms <- list()

for(x in 1:length(METADATAS)){
  meta.data = read_excel(METADATAS[x]) %>% 
    janitor::clean_names()
  vegetation_data = read_excel("Presence_absence_vegetation_AC_Danielsen.xlsx")%>% 
    janitor::clean_names()
  vegetation_data_no_ID = subset(vegetation_data, select = -plot)
  env.data = subset(meta.data, select = -c(order))
  env.data <- env.data %>% tidyr::drop_na()  
  
  Vars <- colnames(env.data)
  Dataset <- "JaccardDistance"
  Response = env.data
  
  Forms <- list()
  
  Models <- for(i in 1:floor(nrow(env.data)/10)){
    Test <- combn(Vars, i, simplify = F)
    cl <- makeCluster(21)
    registerDoParallel(cl)
    Formulas <- foreach(j = 1:length(Test), .combine = "rbind", .packages = c("dplyr")) %dopar% {
      Dataset <- "JaccardDistance"
      DF <- data.frame(Form = NA, AICc = NA)
      Temp <- paste(Dataset,"~", paste(Test[[j]], collapse = " + ")) 
      DF$Form <- Temp
      DF <- DF %>% 
        mutate(Dataset = METADATAS[x])
      gc()
      DF 
    }
    stopCluster(cl)
    message(paste(i, "of", floor(nrow(env.data)/10), "ready", Sys.time()))
    Forms[[i]] <- Formulas
  }
  
  AllForms[[x]] <- Forms %>% 
    purrr::reduce(bind_rows) 
  print(paste(x, "of", length(METADATAS), "ready", Sys.time()))
  
  Dataset <- "JaccardDistance"
}
#> [1] "1 of 3 ready 2023-01-12 11:17:50"
#> [1] "2 of 3 ready 2023-01-12 11:18:58"
#> [1] "3 of 3 ready 2023-01-12 11:20:05"


AllForms <- AllForms %>% 
  purrr::reduce(bind_rows) %>% 
  dplyr::distinct(Form, AICc, .keep_all = T) %>%
  dplyr::mutate(Max_VIF = NA)

NullMod <- data.frame(Form = paste(Dataset, "~ 1", collapse = ""), AICc = NA) %>% 
  mutate(Dataset = METADATAS[1]) %>%
  dplyr::mutate(Max_VIF = NA)

AllForms <- AllForms %>% 
  bind_rows(NullMod)
  

saveRDS(AllForms, "AllForms.rds")
openxlsx::write.xlsx(AllForms, "AllForms.xlsx")
```

</details>

This generate up to 3,504 models to evaluate, which can be downloaded as
an excel file
[here](https://github.com/Sustainscapes/AICcPermanova/raw/master/AllForms.xlsx)
an rds
[here](https://github.com/Sustainscapes/AICcPermanova/blob/master/AllForms.rds).

### 2.1.2 Model fitting

Then in the following code each model is fitted and AICc is calculated
to order it

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Model fitting code

</summary>

``` r
if(file.exists("log.txt")){
  file.remove("log.txt")
}
#> [1] TRUE

cl <- makeCluster(21)
registerDoParallel(cl)

Fs <- foreach(x = 1:nrow(AllForms), .packages = c("vegan", "dplyr", "tidyr", "readxl", "car", "janitor"), .combine = bind_rows) %dopar% {
  meta.data = read_excel(AllForms$Dataset[x]) %>% 
    janitor::clean_names()
  
  vegetation_data = read_excel("Presence_absence_vegetation_AC_Danielsen.xlsx")%>% 
    janitor::clean_names()
  # Remocing the first column (ID) in vegetation-dataset and cheking there are no NAs
  
  vegetation_data_no_ID = subset(vegetation_data, select = -plot)
  
  JaccardDistance <- vegan::vegdist(vegetation_data_no_ID, method = "jaccard")
  
  
  # Removing columns from env.data that is not used in the analysis 
  
  env.data = subset(meta.data, select = -c(order))
  
  
  
  env.data <- env.data %>% tidyr::drop_na()  
  
  Vars <- colnames(env.data)
  Dataset <- "JaccardDistance"
  
  Response = env.data
  

  gc()
  AICc.PERMANOVA2 <- function(adonis2.model) {
    
    # check to see if object is an adonis2 model...
    
    if (is.na(adonis2.model$SumOfSqs[1]))
      stop("object not output of adonis2 {vegan} ")
    
    # Ok, now extract appropriate terms from the adonis model Calculating AICc
    # using residual sum of squares (RSS or SSE) since I don't think that adonis
    # returns something I can use as a likelihood function... maximum likelihood
    # and MSE estimates are the same when distribution is gaussian See e.g.
    # https://www.jessicayung.com/mse-as-maximum-likelihood/;
    # https://towardsdatascience.com/probability-concepts-explained-maximum-likelihood-estimation-c7b4342fdbb1
    # So using RSS or MSE estimates is fine as long as the residuals are
    # Gaussian https://robjhyndman.com/hyndsight/aic/ If models have different
    # conditional likelihoods then AIC is not valid. However, comparing models
    # with different error distributions is ok (above link).
    
    
    RSS <- adonis2.model$SumOfSqs[ length(adonis2.model$SumOfSqs) - 1 ]
    MSE <- RSS / adonis2.model$Df[ length(adonis2.model$Df) - 1 ]
    
    nn <- adonis2.model$Df[ length(adonis2.model$Df) ] + 1
    
    k <- nn - adonis2.model$Df[ length(adonis2.model$Df) - 1 ]
    
    
    # AIC : 2*k + n*ln(RSS/n)
    # AICc: AIC + [2k(k+1)]/(n-k-1)
    
    # based on https://en.wikipedia.org/wiki/Akaike_information_criterion;
    # https://www.statisticshowto.datasciencecentral.com/akaikes-information-criterion/ ;
    # https://www.researchgate.net/post/What_is_the_AIC_formula;
    # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf;
    # https://medium.com/better-programming/data-science-modeling-how-to-use-linear-regression-with-python-fdf6ca5481be 
    
    
    AIC <- 2*k + nn*log(RSS/nn)
    AICc <- AIC + (2*k*(k + 1))/(nn - k - 1)
    
    output <- data.frame(AICc = AICc, k = k, N = nn)
    
    return(output)   
    
  }
  Temp <- AllForms[x,]
  Temp$AICc <-  try(AICc.PERMANOVA2(adonis2(as.formula(AllForms$Form[x]), data = Response, by = "margin"))$AICc, silent = T)
  
  Response$y <- rnorm(n = nrow(Response))
  
  VIF <- function(model) {
    tryCatch({
        vif <- car::vif(model)
        max(vif)
    }, error = function(e) {
        if (grepl("aliased coefficients", e$message)) {
            20000
        } else if (grepl("model contains fewer than 2 terms", e$message)) {
            0
        } else {
            stop(e)
        }
    })
}
  
  Temp$Max_VIF <- VIF(lm(as.formula(stringr::str_replace_all(AllForms$Form[x], "JaccardDistance ", "y")), data = Response))
  
  Rs <- broom::tidy(adonis2(as.formula(AllForms$Form[x]), data = Response, by = "margin")) %>% dplyr::filter(!(term %in% c("Residual", "Total"))) %>% dplyr::select(term, R2) %>%  pivot_wider(names_from = term, values_from = R2)
  if((x %% 100) == 0){
    sink("log.txt", append = T)
    cat(paste("finished", x, "number of models", Sys.time(), "of",  nrow(AllForms)))
    cat("\n")
    sink()
  }
  
  bind_cols(Temp, Rs)
}

stopCluster(cl)

saveRDS(Fs, "FS.rds")

Fs <- Fs %>% arrange(AICc)

saveRDS(Fs, "FS.rds")
openxlsx::write.xlsx(Fs, "FS.xlsx")
```

</details>

As seen in table <a href="#tab:SummaryPlantPA">2.1</a> there are 35
models within 2 AICc where the max VIF is lower or equal than 6 of each
other, you can see there how many times a variable has been selected

| Variable          | Number_of_models | Full_Akaike_Adjusted_RSq | Subset_Akaike_Adjusted_RSq |
|:------------------|-----------------:|-------------------------:|---------------------------:|
| habitat_type      |               35 |                    0.372 |                      0.372 |
| fine_silt         |               11 |                    0.009 |                      0.028 |
| ec                |               10 |                    0.009 |                      0.030 |
| oc_beregnet       |                9 |                    0.008 |                      0.030 |
| clay              |                6 |                    0.003 |                      0.021 |
| wr                |                5 |                    0.003 |                      0.023 |
| nitrat_nitrit     |                5 |                    0.002 |                      0.017 |
| dexter_n          |                5 |                    0.002 |                      0.018 |
| finesilt_and_clay |                4 |                    0.003 |                      0.025 |
| coarse_silt_sand  |                3 |                    0.002 |                      0.022 |
| water_content     |                2 |                    0.001 |                      0.024 |
| ammonium          |                1 |                    0.000 |                      0.009 |

Table 2.1: Number of selected models were variables are present and
their Akaike Weighted R squared for the Marginal effect of the terms

Now we can see the top models that have a delta AICc within 2 of the
best model in table <a href="#tab:BestPLantModels">2.2</a> if expanded

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Show table of selected models for vegetation presence absence

</summary>

| Form                                                              |    AICc | DeltaAICc | Max_VIF | habitat_type |    ec | water_content |    wr | ammonium | nitrat_nitrit | oc_beregnet |  clay | fine_silt | coarse_silt_sand | finesilt_and_clay | dexter_n |
|:------------------------------------------------------------------|--------:|----------:|--------:|-------------:|------:|--------------:|------:|---------:|--------------:|------------:|------:|----------:|-----------------:|------------------:|---------:|
| JaccardDistance \~ habitat_type                                   | -60.591 |     0.000 |       0 |        0.411 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec                              | -60.371 |     0.220 |       6 |        0.351 | 0.029 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + fine_silt                  | -60.298 |     0.293 |       6 |        0.337 | 0.032 |            NA |    NA |       NA |            NA |          NA |    NA |     0.030 |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + fine_silt                       | -60.217 |     0.374 |       6 |        0.385 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.027 |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet                     | -60.055 |     0.536 |       6 |        0.403 |    NA |            NA |    NA |       NA |            NA |       0.025 |    NA |        NA |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content + oc_beregnet     | -60.002 |     0.589 |       6 |        0.382 |    NA |         0.030 |    NA |       NA |            NA |       0.040 |    NA |        NA |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + fine_silt         | -59.990 |     0.601 |       6 |        0.379 |    NA |            NA |    NA |       NA |            NA |       0.029 |    NA |     0.030 |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + finesilt_and_clay               | -59.913 |     0.678 |       6 |        0.387 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |             0.024 |       NA |
| JaccardDistance \~ habitat_type + wr                              | -59.843 |     0.747 |       6 |        0.368 |    NA |            NA | 0.023 |       NA |            NA |          NA |    NA |        NA |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + finesilt_and_clay          | -59.789 |     0.802 |       6 |        0.336 | 0.030 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |             0.025 |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + coarse_silt_sand  | -59.765 |     0.826 |       6 |        0.382 |    NA |            NA |    NA |       NA |            NA |       0.035 |    NA |        NA |            0.028 |                NA |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + finesilt_and_clay | -59.765 |     0.826 |       6 |        0.382 |    NA |            NA |    NA |       NA |            NA |       0.030 |    NA |        NA |               NA |             0.028 |       NA |
| JaccardDistance \~ habitat_type + clay                            | -59.562 |     1.029 |       6 |        0.394 |    NA |            NA |    NA |       NA |            NA |          NA | 0.020 |        NA |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + wr + fine_silt                  | -59.492 |     1.098 |       6 |        0.356 |    NA |            NA | 0.023 |       NA |            NA |          NA |    NA |     0.027 |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + coarse_silt_sand                | -59.425 |     1.165 |       6 |        0.391 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.018 |                NA |       NA |
| JaccardDistance \~ habitat_type + dexter_n                        | -59.409 |     1.182 |       6 |        0.409 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |                NA |    0.018 |
| JaccardDistance \~ habitat_type + nitrat_nitrit                   | -59.345 |     1.246 |       6 |        0.391 |    NA |            NA |    NA |       NA |         0.017 |          NA |    NA |        NA |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + clay              | -59.332 |     1.259 |       6 |        0.388 |    NA |            NA |    NA |       NA |            NA |       0.029 | 0.023 |        NA |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + clay + fine_silt           | -59.314 |     1.277 |       6 |        0.335 | 0.032 |            NA |    NA |       NA |            NA |          NA | 0.021 |     0.032 |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + clay + fine_silt                | -59.260 |     1.331 |       6 |        0.383 |    NA |            NA |    NA |       NA |            NA |          NA | 0.021 |     0.028 |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + clay                       | -59.258 |     1.333 |       6 |        0.339 | 0.028 |            NA |    NA |       NA |            NA |          NA | 0.019 |        NA |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + dexter_n                   | -59.223 |     1.368 |       6 |        0.351 | 0.029 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |                NA |    0.019 |
| JaccardDistance \~ habitat_type + wr + oc_beregnet                | -59.192 |     1.399 |       6 |        0.362 |    NA |            NA | 0.022 |       NA |            NA |       0.024 |    NA |        NA |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content                   | -59.183 |     1.408 |       6 |        0.384 |    NA |         0.016 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + nitrat_nitrit              | -59.141 |     1.450 |       6 |        0.347 | 0.029 |            NA |    NA |       NA |         0.018 |          NA |    NA |        NA |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + coarse_silt_sand           | -59.136 |     1.454 |       6 |        0.326 | 0.028 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.018 |                NA |       NA |
| JaccardDistance \~ habitat_type + wr + finesilt_and_clay          | -59.120 |     1.471 |       6 |        0.355 |    NA |            NA | 0.023 |       NA |            NA |          NA |    NA |        NA |               NA |             0.024 |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + clay + fine_silt  | -59.044 |     1.547 |       6 |        0.377 |    NA |            NA |    NA |       NA |            NA |       0.029 | 0.021 |     0.028 |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + fine_silt + dexter_n       | -59.030 |     1.561 |       6 |        0.337 | 0.033 |            NA |    NA |       NA |            NA |          NA |    NA |     0.029 |               NA |                NA |    0.018 |
| JaccardDistance \~ habitat_type + fine_silt + dexter_n            | -58.878 |     1.713 |       6 |        0.386 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.026 |               NA |                NA |    0.017 |
| JaccardDistance \~ habitat_type + nitrat_nitrit + oc_beregnet     | -58.800 |     1.791 |       6 |        0.382 |    NA |            NA |    NA |       NA |         0.018 |       0.026 |    NA |        NA |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + nitrat_nitrit + fine_silt  | -58.782 |     1.809 |       6 |        0.338 | 0.032 |            NA |    NA |       NA |         0.016 |          NA |    NA |     0.028 |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + nitrat_nitrit + fine_silt       | -58.698 |     1.893 |       6 |        0.382 |    NA |            NA |    NA |       NA |         0.015 |          NA |    NA |     0.025 |               NA |                NA |       NA |
| JaccardDistance \~ habitat_type + wr + dexter_n                   | -58.602 |     1.989 |       6 |        0.367 |    NA |            NA | 0.023 |       NA |            NA |          NA |    NA |        NA |               NA |                NA |    0.018 |
| JaccardDistance \~ habitat_type + ammonium                        | -58.594 |     1.997 |       6 |        0.389 |    NA |            NA |    NA |    0.009 |            NA |          NA |    NA |        NA |               NA |                NA |       NA |

Table 2.2: Best models

</details>

## 2.2 Abundance data

### 2.2.1 Generation of all possible models

After that we read in the datasets for environmental layers and generate
all possible models to fit, in this case we will limit ourselves to only
using at most one variable per ten observations, in this case that means
up to 5 variables per model. The code for generating all possible models
can be expanded bellow

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Model generator

</summary>

``` r
METADATAS <- list.files(pattern = "PERMANOVA_VEGETATION_", full.names = T)

AllForms <- list()

for(x in 1:length(METADATAS)){
  meta.data = read_excel(METADATAS[x]) %>% 
    janitor::clean_names()
  vegetation_data = read_excel("Pinpoint-data til ordination_minus_Lønstrup.xlsx")%>% 
    janitor::clean_names()
  vegetation_data_no_ID = vegetation_data
  env.data = subset(meta.data, select = -c(order))
  env.data <- env.data %>% tidyr::drop_na()  
  
  Vars <- colnames(env.data)
  Dataset <- "BrayDistance"
  Response = env.data
  
  Forms <- list()
  
  Models <- for(i in 1:floor(nrow(env.data)/10)){
    Test <- combn(Vars, i, simplify = F)
    cl <- makeCluster(21)
    registerDoParallel(cl)
    Formulas <- foreach(j = 1:length(Test), .combine = "rbind", .packages = c("dplyr")) %dopar% {
      Dataset <- "BrayDistance"
      DF <- data.frame(Form = NA, AICc = NA)
      Temp <- paste(Dataset,"~", paste(Test[[j]], collapse = " + ")) 
      DF$Form <- Temp
      DF <- DF %>% 
        mutate(Dataset = METADATAS[x])
      gc()
      DF  
    }
    stopCluster(cl)
    message(paste(i, "of", floor(nrow(env.data)/10), "ready", Sys.time()))
    Forms[[i]] <- Formulas
  }
  
  AllForms[[x]] <- Forms %>% 
    purrr::reduce(bind_rows) 
  print(paste(x, "of", length(METADATAS), "ready", Sys.time()))
  
  Dataset <- "BrayDistance"
}
#> [1] "1 of 3 ready 2023-01-12 11:29:27"
#> [1] "2 of 3 ready 2023-01-12 11:30:34"
#> [1] "3 of 3 ready 2023-01-12 11:31:41"


AllForms <- AllForms %>% 
  purrr::reduce(bind_rows) %>% 
  dplyr::distinct(Form, AICc, .keep_all = T) %>%
  dplyr::mutate(Max_VIF = NA)

NullMod <- data.frame(Form = paste(Dataset, "~ 1", collapse = ""), AICc = NA) %>% 
  mutate(Dataset = METADATAS[1]) %>%
  dplyr::mutate(Max_VIF = NA)

AllFormsVegAbund <- AllForms %>% 
  bind_rows(NullMod)
  

saveRDS(AllFormsVegAbund, "AllFormsVegAbund.rds")
openxlsx::write.xlsx(AllFormsVegAbund, "AllFormsVegAbund.xlsx")
```

</details>

This generate up to 3,504 models to evaluate, which can be downloaded as
an excel file
[here](https://github.com/Sustainscapes/AICcPermanova/raw/master/AllFormsVegAbund.xlsx)
an rds
[here](https://github.com/Sustainscapes/AICcPermanova/blob/master/AllFormsVegAbund.rds).

### 2.2.2 Model fitting

Then in the following code each model is fitted and AICc is calculated
to order it

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Model fitting code for vegetation abundance

</summary>

``` r
if(file.exists("logVegAbund.txt")){
  file.remove("logVegAbund.txt")
}
#> [1] TRUE

cl <- makeCluster(21)
registerDoParallel(cl)

Fs <- foreach(x = 1:nrow(AllFormsVegAbund), .packages = c("vegan", "dplyr", "tidyr", "readxl", "car", "janitor"), .combine = bind_rows) %dopar% {
  meta.data = read_excel(AllFormsVegAbund$Dataset[x]) %>% 
    janitor::clean_names()
  
  vegetation_data = read_excel("Pinpoint-data til ordination_minus_Lønstrup.xlsx")%>% 
    janitor::clean_names()
  # Remocing the first column (ID) in vegetation-dataset and cheking there are no NAs
  
  vegetation_data_no_ID = vegetation_data
  
  BrayDistance <- vegan::vegdist(vegetation_data_no_ID, method = "bray")
  
  
  # Removing columns from env.data that is not used in the analysis 
  
  env.data = subset(meta.data, select = -c(order))
  
  
  
  env.data <- env.data %>% tidyr::drop_na()  
  
  Vars <- colnames(env.data)
  Dataset <- "BrayDistance"
  
  Response = env.data
  

  gc()
  AICc.PERMANOVA2 <- function(adonis2.model) {
    
    # check to see if object is an adonis2 model...
    
    if (is.na(adonis2.model$SumOfSqs[1]))
      stop("object not output of adonis2 {vegan} ")
    
    # Ok, now extract appropriate terms from the adonis model Calculating AICc
    # using residual sum of squares (RSS or SSE) since I don't think that adonis
    # returns something I can use as a likelihood function... maximum likelihood
    # and MSE estimates are the same when distribution is gaussian See e.g.
    # https://www.jessicayung.com/mse-as-maximum-likelihood/;
    # https://towardsdatascience.com/probability-concepts-explained-maximum-likelihood-estimation-c7b4342fdbb1
    # So using RSS or MSE estimates is fine as long as the residuals are
    # Gaussian https://robjhyndman.com/hyndsight/aic/ If models have different
    # conditional likelihoods then AIC is not valid. However, comparing models
    # with different error distributions is ok (above link).
    
    
    RSS <- adonis2.model$SumOfSqs[ length(adonis2.model$SumOfSqs) - 1 ]
    MSE <- RSS / adonis2.model$Df[ length(adonis2.model$Df) - 1 ]
    
    nn <- adonis2.model$Df[ length(adonis2.model$Df) ] + 1
    
    k <- nn - adonis2.model$Df[ length(adonis2.model$Df) - 1 ]
    
    
    # AIC : 2*k + n*ln(RSS/n)
    # AICc: AIC + [2k(k+1)]/(n-k-1)
    
    # based on https://en.wikipedia.org/wiki/Akaike_information_criterion;
    # https://www.statisticshowto.datasciencecentral.com/akaikes-information-criterion/ ;
    # https://www.researchgate.net/post/What_is_the_AIC_formula;
    # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf;
    # https://medium.com/better-programming/data-science-modeling-how-to-use-linear-regression-with-python-fdf6ca5481be 
    
    
    AIC <- 2*k + nn*log(RSS/nn)
    AICc <- AIC + (2*k*(k + 1))/(nn - k - 1)
    
    output <- data.frame(AICc = AICc, k = k, N = nn)
    
    return(output)   
    
  }
  Temp <- AllFormsVegAbund[x,]
  Temp$AICc <-  try(AICc.PERMANOVA2(adonis2(as.formula(AllFormsVegAbund$Form[x]), data = Response, by = "margin"))$AICc, silent = T)
  
  Response$y <- rnorm(n = nrow(Response))
  
  VIF <- function(model) {
    tryCatch({
        vif <- car::vif(model)
        max(vif)
    }, error = function(e) {
        if (grepl("aliased coefficients", e$message)) {
            20000
        } else if (grepl("model contains fewer than 2 terms", e$message)) {
            0
        } else {
            stop(e)
        }
    })
}
  
  Temp$Max_VIF <- VIF(lm(as.formula(stringr::str_replace_all(AllFormsVegAbund$Form[x], "BrayDistance", "y")), data = Response))
  
  Rs <- broom::tidy(adonis2(as.formula(AllFormsVegAbund$Form[x]), data = Response, by = "margin")) %>% dplyr::filter(!(term %in% c("Residual", "Total"))) %>% dplyr::select(term, R2) %>%  pivot_wider(names_from = term, values_from = R2)
  if((x %% 100) == 0){
    sink("logVegAbund.txt", append = T)
    cat(paste("finished", x, "number of models", Sys.time(), "of",  nrow(AllFormsVegAbund)))
    cat("\n")
    sink()
  }
  
  bind_cols(Temp, Rs)
}

stopCluster(cl)

saveRDS(Fs, "FSVegAbund.rds")

Fs <- Fs %>% arrange(AICc)

saveRDS(Fs, "FSVegAbund.rds")
openxlsx::write.xlsx(Fs, "FSVegAbund.xlsx")
```

</details>

As seen in table <a href="#tab:SummaryVegAbund">2.3</a> there are 12
models within 2 AICc of each other, you can see there how many times a
variable has been selected

| Variable      | Number_of_models | Full_Akaike_Adjusted_RSq | Subset_Akaike_Adjusted_RSq |
|:--------------|-----------------:|-------------------------:|---------------------------:|
| habitat_type  |               12 |                    0.472 |                      0.472 |
| ec            |                4 |                    0.009 |                      0.025 |
| nitrat_nitrit |                3 |                    0.003 |                      0.014 |
| ammonium      |                2 |                    0.002 |                      0.013 |
| oc_beregnet   |                2 |                    0.003 |                      0.018 |
| fine_silt     |                2 |                    0.001 |                      0.011 |
| water_content |                1 |                    0.001 |                      0.010 |
| wr            |                1 |                    0.001 |                      0.014 |

Table 2.3: Number of selected models were variables are present and
their Akaike Weighted R squared for the Marginal effect of the terms

Now we can see the top models that have a delta AICc within 2 of the
best model in table <a href="#tab:BestVegAbundModels">2.4</a> if
expanded

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Show table of selected models for vegetation presence absence

</summary>

| Form                                                       |    AICc | DeltaAICc | Max_VIF | habitat_type |    ec | water_content |    wr | ammonium | nitrat_nitrit | oc_beregnet | fine_silt |
|:-----------------------------------------------------------|--------:|----------:|--------:|-------------:|------:|--------------:|------:|---------:|--------------:|------------:|----------:|
| BrayDistance \~ habitat_type + ec                          | -75.760 |     0.000 |       6 |        0.448 | 0.025 |            NA |    NA |       NA |            NA |          NA |        NA |
| BrayDistance \~ habitat_type                               | -75.690 |     0.070 |       0 |        0.538 |    NA |            NA |    NA |       NA |            NA |          NA |        NA |
| BrayDistance \~ habitat_type + oc_beregnet                 | -74.950 |     0.810 |       6 |        0.480 |    NA |            NA |    NA |       NA |            NA |       0.018 |        NA |
| BrayDistance \~ habitat_type + nitrat_nitrit               | -74.557 |     1.202 |       6 |        0.510 |    NA |            NA |    NA |       NA |         0.015 |          NA |        NA |
| BrayDistance \~ habitat_type + wr                          | -74.475 |     1.284 |       6 |        0.446 |    NA |            NA | 0.014 |       NA |            NA |          NA |        NA |
| BrayDistance \~ habitat_type + ec + ammonium               | -74.463 |     1.297 |       6 |        0.407 | 0.027 |            NA |    NA |    0.014 |            NA |          NA |        NA |
| BrayDistance \~ habitat_type + ec + nitrat_nitrit          | -74.461 |     1.298 |       6 |        0.445 | 0.024 |            NA |    NA |       NA |         0.014 |          NA |        NA |
| BrayDistance \~ habitat_type + ec + fine_silt              | -74.223 |     1.537 |       6 |        0.431 | 0.026 |            NA |    NA |       NA |            NA |          NA |     0.012 |
| BrayDistance \~ habitat_type + ammonium                    | -74.199 |     1.561 |       6 |        0.497 |    NA |            NA |    NA |    0.012 |            NA |          NA |        NA |
| BrayDistance \~ habitat_type + fine_silt                   | -74.123 |     1.637 |       6 |        0.499 |    NA |            NA |    NA |       NA |            NA |          NA |     0.011 |
| BrayDistance \~ habitat_type + water_content               | -74.022 |     1.737 |       6 |        0.475 |    NA |          0.01 |    NA |       NA |            NA |          NA |        NA |
| BrayDistance \~ habitat_type + nitrat_nitrit + oc_beregnet | -73.807 |     1.953 |       6 |        0.450 |    NA |            NA |    NA |       NA |         0.015 |       0.018 |        NA |

Table 2.4: Best models for vegetation abundance

</details>

# 3 Baterial models

## 3.1 Abundance data

### 3.1.1 Generation of all possible models

We read in the datasets for environmental layers and generate all
possible models to fit, in this case we will limit ourselves to only
using at most one variable per ten observations, in this case that means
up to 5 variables per model. The code for generating all possible models
can be expanded bellow

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Model generator bacterial abundance

</summary>

``` r
METADATAS <- list.files(pattern = "PERMANOVA_BACTERIA_", full.names = T)

AllForms <- list()

for(x in 1:length(METADATAS)){
  meta.data = read_excel(METADATAS[x]) %>% 
    janitor::clean_names()
  d <- amp_load(
  otutable = "AC_otutale_new.txt",
  metadata = "Metadata_mix-samples_AC_Danielsen_final.xlsx")
  bac.data.subset = amp_subset_taxa(d, "d__Bacteria")
  bacterial_data = amp_subset_samples(d, Investigator == "AC")
  bacterial_data = as.data.frame(t(bacterial_data$abund)) %>% 
    janitor::clean_names()

  env.data = subset(meta.data, select = -c(order))
  env.data <- env.data %>% tidyr::drop_na()  
  
  Vars <- colnames(env.data)
  Dataset <- "BrayDistance"
  Response = env.data
  
  Forms <- list()
  
  Models <- for(i in 1:floor(nrow(env.data)/10)){
    Test <- combn(Vars, i, simplify = F)
    cl <- makeCluster(21)
    registerDoParallel(cl)
    Formulas <- foreach(j = 1:length(Test), .combine = "rbind", .packages = c("dplyr")) %dopar% {
      Dataset <- "BrayDistance"
      DF <- data.frame(Form = NA, AICc = NA)
      Temp <- paste(Dataset,"~", paste(Test[[j]], collapse = " + ")) 
      DF$Form <- Temp
      DF <- DF %>% 
        mutate(Dataset = METADATAS[x])
      gc()
      DF 
    }
    stopCluster(cl)
    message(paste(i, "of", floor(nrow(env.data)/10), "ready", Sys.time()))
    Forms[[i]] <- Formulas
  }
  
  AllForms[[x]] <- Forms %>% 
    purrr::reduce(bind_rows) 
  print(paste(x, "of", length(METADATAS), "ready", Sys.time()))
  
  Dataset <- "BrayDistance"
}
#> [1] "1 of 3 ready 2023-01-12 11:53:44"
#> [1] "2 of 3 ready 2023-01-12 11:55:02"
#> [1] "3 of 3 ready 2023-01-12 11:56:20"


AllForms <- AllForms %>% 
  purrr::reduce(bind_rows) %>% 
  dplyr::distinct(Form, AICc, .keep_all = T) %>%
  dplyr::mutate(Max_VIF = NA)

NullMod <- data.frame(Form = paste(Dataset, "~ 1", collapse = ""), AICc = NA) %>% 
  mutate(Dataset = METADATAS[1]) %>%
  dplyr::mutate(Max_VIF = NA)

AllForms <- AllForms %>% 
  bind_rows(NullMod)
  

saveRDS(AllForms, "AllFormsBacterialAbund.rds")
openxlsx::write.xlsx(AllForms, "AllFormsBacterialAbund.xlsx")
```

</details>

This generate up to 5,061 models to evaluate, which can be downloaded as
an excel file
[here](https://github.com/Sustainscapes/AICcPermanova/raw/master/AllFormsBacterialAbund.xlsx)
an rds
[here](https://github.com/Sustainscapes/AICcPermanova/blob/master/AllFormsBacterialAbund.rds).

### 3.1.2 Model fitting

Then in the following code each model is fitted and AICc is calculated
to order it

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Model fitting code bacterial abund

</summary>

``` r
if(file.exists("logBacterialAbund.txt")){
  file.remove("logBacterialAbund.txt")
}
#> [1] TRUE

cl <- makeCluster(21)
registerDoParallel(cl)

Fs <- foreach(x = 1:nrow(AllForms), .packages = c("vegan", "dplyr", "tidyr", "readxl", "car", "janitor", "ampvis2"), .combine = bind_rows) %dopar% {
  AllForms <- readRDS("AllFormsBacterialAbund.rds")
  meta.data = read_excel(AllForms$Dataset[x]) %>% 
    janitor::clean_names()
  
  d <- amp_load(
  otutable = "AC_otutale_new.txt",
  metadata = "Metadata_mix-samples_AC_Danielsen_final.xlsx")
  bac.data.subset = amp_subset_taxa(d, "d__Bacteria")
  bacterial_data = amp_subset_samples(d, Investigator == "AC")
  bacterial_data = as.data.frame(t(bacterial_data$abund)) %>% 
    janitor::clean_names()
  # Remocing the first column (ID) in vegetation-dataset and cheking there are no NAs
  meta.data <- meta.data[match(rownames(bacterial_data), meta.data$seq_id),]
  BrayDistance <- vegan::vegdist(bacterial_data, method = "bray")
  
  
  # Removing columns from env.data that is not used in the analysis 
  
  env.data = subset(meta.data, select = -c(order))
  
  
  
  env.data <- env.data %>% tidyr::drop_na()  
  
  Vars <- colnames(env.data)
  Dataset <- "BrayDistance"
  
  Response = env.data
  

  gc()
  AICc.PERMANOVA2 <- function(adonis2.model) {
    
    # check to see if object is an adonis2 model...
    
    if (is.na(adonis2.model$SumOfSqs[1]))
      stop("object not output of adonis2 {vegan} ")
    
    # Ok, now extract appropriate terms from the adonis model Calculating AICc
    # using residual sum of squares (RSS or SSE) since I don't think that adonis
    # returns something I can use as a likelihood function... maximum likelihood
    # and MSE estimates are the same when distribution is gaussian See e.g.
    # https://www.jessicayung.com/mse-as-maximum-likelihood/;
    # https://towardsdatascience.com/probability-concepts-explained-maximum-likelihood-estimation-c7b4342fdbb1
    # So using RSS or MSE estimates is fine as long as the residuals are
    # Gaussian https://robjhyndman.com/hyndsight/aic/ If models have different
    # conditional likelihoods then AIC is not valid. However, comparing models
    # with different error distributions is ok (above link).
    
    
    RSS <- adonis2.model$SumOfSqs[ length(adonis2.model$SumOfSqs) - 1 ]
    MSE <- RSS / adonis2.model$Df[ length(adonis2.model$Df) - 1 ]
    
    nn <- adonis2.model$Df[ length(adonis2.model$Df) ] + 1
    
    k <- nn - adonis2.model$Df[ length(adonis2.model$Df) - 1 ]
    
    
    # AIC : 2*k + n*ln(RSS/n)
    # AICc: AIC + [2k(k+1)]/(n-k-1)
    
    # based on https://en.wikipedia.org/wiki/Akaike_information_criterion;
    # https://www.statisticshowto.datasciencecentral.com/akaikes-information-criterion/ ;
    # https://www.researchgate.net/post/What_is_the_AIC_formula;
    # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf;
    # https://medium.com/better-programming/data-science-modeling-how-to-use-linear-regression-with-python-fdf6ca5481be 
    
    
    AIC <- 2*k + nn*log(RSS/nn)
    AICc <- AIC + (2*k*(k + 1))/(nn - k - 1)
    
    output <- data.frame(AICc = AICc, k = k, N = nn)
    
    return(output)   
    
  }
  Temp <- AllForms[x,]
  Temp$AICc <-  try(AICc.PERMANOVA2(adonis2(as.formula(AllForms$Form[x]), data = Response, by = "margin"))$AICc, silent = T)
  
  Response$y <- rnorm(n = nrow(Response))
  
  VIF <- function(model) {
    tryCatch({
        vif <- car::vif(model)
        max(vif)
    }, error = function(e) {
        if (grepl("aliased coefficients", e$message)) {
            20000
        } else if (grepl("model contains fewer than 2 terms", e$message)) {
            0
        } else {
            stop(e)
        }
    })
}
  
  Temp$Max_VIF <- VIF(lm(as.formula(stringr::str_replace_all(AllForms$Form[x], "BrayDistance", "y")), data = Response))
  
  Rs <- broom::tidy(adonis2(as.formula(AllForms$Form[x]), data = Response, by = "margin")) %>% dplyr::filter(!(term %in% c("Residual", "Total"))) %>% dplyr::select(term, R2) %>%  pivot_wider(names_from = term, values_from = R2)
  if((x %% 100) == 0){
    sink("logBacterialAbund.txt", append = T)
    cat(paste("finished", x, "number of models", Sys.time(), "of",  nrow(AllForms)))
    cat("\n")
    sink()
  }
  
  bind_cols(Temp, Rs)
}

stopCluster(cl)

saveRDS(Fs, "FSBacterialAbund.rds")

Fs <- Fs %>% arrange(AICc)

saveRDS(Fs, "FSBacterialAbund.rds")
openxlsx::write.xlsx(Fs, "FSBacterialAbund.xlsx")
```

</details>

As seen in table <a href="#tab:SummaryBacterialAbund">3.1</a> there are
26 models within 2 AICc of each other, you can see there how many times
a variable has been selected

| Variable          | Number_of_models | Full_Akaike_Adjusted_RSq | Subset_Akaike_Adjusted_RSq |
|:------------------|-----------------:|-------------------------:|---------------------------:|
| habitat_type      |               26 |                    0.456 |                      0.456 |
| oc_beregnet       |                9 |                    0.010 |                      0.025 |
| water_content     |                5 |                    0.005 |                      0.022 |
| ammonium          |                5 |                    0.002 |                      0.013 |
| shannon_veg       |                5 |                    0.004 |                      0.021 |
| coarse_silt_sand  |                4 |                    0.002 |                      0.016 |
| wr                |                3 |                    0.003 |                      0.020 |
| clay              |                2 |                    0.001 |                      0.012 |
| fine_silt         |                2 |                    0.001 |                      0.012 |
| finesilt_and_clay |                2 |                    0.001 |                      0.012 |
| ec                |                1 |                    0.000 |                      0.015 |
| dexter_n          |                1 |                    0.000 |                      0.012 |

Table 3.1: Number of selected models were variables are present and
their Akaike Weighted R squared for the Marginal effect of the terms

Now we can see the top models that have a delta AICc within 2 of the
best model in table <a href="#tab:BestBacterialAbundModels">3.2</a> if
expanded

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Show table of selected models for bacterial abundance

</summary>

| Form                                                            |    AICc | DeltaAICc | Max_VIF | habitat_type |    ec | water_content |    wr | ammonium | oc_beregnet |  clay | fine_silt | coarse_silt_sand | shannon_veg | finesilt_and_clay | dexter_n |
|:----------------------------------------------------------------|--------:|----------:|--------:|-------------:|------:|--------------:|------:|---------:|------------:|------:|----------:|-----------------:|------------:|------------------:|---------:|
| BrayDistance \~ habitat_type + oc_beregnet                      | -86.086 |     0.000 |       6 |        0.476 |    NA |            NA |    NA |       NA |       0.028 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + oc_beregnet + shannon_veg        | -85.789 |     0.297 |       6 |        0.454 |    NA |            NA |    NA |       NA |       0.028 |    NA |        NA |               NA |       0.022 |                NA |       NA |
| BrayDistance \~ habitat_type                                    | -85.730 |     0.355 |       0 |        0.529 |    NA |            NA |    NA |       NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + water_content                    | -85.591 |     0.495 |       6 |        0.489 |    NA |         0.024 |    NA |       NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + water_content + oc_beregnet      | -85.582 |     0.504 |       6 |        0.449 |    NA |         0.020 |    NA |       NA |       0.025 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + wr                               | -85.547 |     0.539 |       6 |        0.392 |    NA |            NA | 0.023 |       NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + shannon_veg                      | -85.331 |     0.755 |       6 |        0.511 |    NA |            NA |    NA |       NA |          NA |    NA |        NA |               NA |       0.021 |                NA |       NA |
| BrayDistance \~ habitat_type + wr + oc_beregnet                 | -85.149 |     0.937 |       6 |        0.359 |    NA |            NA | 0.017 |       NA |       0.021 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + water_content + shannon_veg      | -85.047 |     1.039 |       6 |        0.466 |    NA |         0.022 |    NA |       NA |          NA |    NA |        NA |               NA |       0.020 |                NA |       NA |
| BrayDistance \~ habitat_type + coarse_silt_sand                 | -84.869 |     1.216 |       6 |        0.490 |    NA |            NA |    NA |       NA |          NA |    NA |        NA |            0.017 |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + ammonium                         | -84.683 |     1.403 |       6 |        0.498 |    NA |            NA |    NA |    0.016 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + water_content + coarse_silt_sand | -84.656 |     1.430 |       6 |        0.449 |    NA |         0.023 |    NA |       NA |          NA |    NA |        NA |            0.017 |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + ec                               | -84.556 |     1.530 |       6 |        0.354 | 0.015 |            NA |    NA |       NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + oc_beregnet + fine_silt          | -84.527 |     1.559 |       6 |        0.416 |    NA |            NA |    NA |       NA |       0.027 |    NA |     0.012 |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + oc_beregnet + coarse_silt_sand   | -84.491 |     1.594 |       6 |        0.418 |    NA |            NA |    NA |       NA |       0.022 |    NA |        NA |            0.011 |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + oc_beregnet + finesilt_and_clay  | -84.491 |     1.594 |       6 |        0.418 |    NA |            NA |    NA |       NA |       0.026 |    NA |        NA |               NA |          NA |             0.011 |       NA |
| BrayDistance \~ habitat_type + oc_beregnet + clay               | -84.448 |     1.638 |       6 |        0.431 |    NA |            NA |    NA |       NA |       0.026 | 0.011 |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + coarse_silt_sand + shannon_veg   | -84.442 |     1.644 |       6 |        0.472 |    NA |            NA |    NA |       NA |          NA |    NA |        NA |            0.017 |       0.021 |                NA |       NA |
| BrayDistance \~ habitat_type + finesilt_and_clay                | -84.339 |     1.746 |       6 |        0.474 |    NA |            NA |    NA |       NA |          NA |    NA |        NA |               NA |          NA |             0.013 |       NA |
| BrayDistance \~ habitat_type + clay                             | -84.328 |     1.757 |       6 |        0.487 |    NA |            NA |    NA |       NA |          NA | 0.013 |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + fine_silt                        | -84.314 |     1.771 |       6 |        0.472 |    NA |            NA |    NA |       NA |          NA |    NA |     0.012 |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + dexter_n                         | -84.266 |     1.820 |       6 |        0.517 |    NA |            NA |    NA |       NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.012 |
| BrayDistance \~ habitat_type + water_content + ammonium         | -84.265 |     1.820 |       6 |        0.476 |    NA |         0.022 |    NA |    0.014 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + wr + ammonium                    | -84.208 |     1.878 |       6 |        0.373 |    NA |            NA | 0.021 |    0.014 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + ammonium + shannon_veg           | -84.175 |     1.911 |       6 |        0.471 |    NA |            NA |    NA |    0.015 |          NA |    NA |        NA |               NA |       0.021 |                NA |       NA |
| BrayDistance \~ habitat_type + ammonium + oc_beregnet           | -84.102 |     1.983 |       6 |        0.471 |    NA |            NA |    NA |    0.008 |       0.020 |    NA |        NA |               NA |          NA |                NA |       NA |

Table 3.2: Best models

</details>

## 3.2 Presence absence data

### 3.2.1 Generation of all possible models

We read in the datasets for environmental layers and generate all
possible models to fit, in this case we will limit ourselves to only
using at most one variable per ten observations, in this case that means
up to 5 variables per model. The code for generating all possible models
can be expanded bellow

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Model generator bacterial abundance

</summary>

``` r
METADATAS <- list.files(pattern = "PERMANOVA_BACTERIA_", full.names = T)

AllForms <- list()

for(x in 1:length(METADATAS)){
  meta.data = read_excel(METADATAS[x]) %>% 
    janitor::clean_names()
  d <- amp_load(
  otutable = "AC_otutale_new.txt",
  metadata = "Metadata_mix-samples_AC_Danielsen_final.xlsx")
  bac.data.subset = amp_subset_taxa(d, "d__Bacteria")
  bacterial_data = amp_subset_samples(d, Investigator == "AC")
  bacterial_data = as.data.frame(t(bacterial_data$abund)) %>% 
    janitor::clean_names()

  env.data = subset(meta.data, select = -c(order))
  env.data <- env.data %>% tidyr::drop_na()  
  
  Vars <- colnames(env.data)
  Dataset <- "JaccardDistance"
  Response = env.data
  
  Forms <- list()
  
  Models <- for(i in 1:floor(nrow(env.data)/10)){
    Test <- combn(Vars, i, simplify = F)
    cl <- makeCluster(21)
    registerDoParallel(cl)
    Formulas <- foreach(j = 1:length(Test), .combine = "rbind", .packages = c("dplyr")) %dopar% {
      Dataset <- "JaccardDistance"
      DF <- data.frame(Form = NA, AICc = NA)
      Temp <- paste(Dataset,"~", paste(Test[[j]], collapse = " + ")) 
      DF$Form <- Temp
      DF <- DF %>% 
        mutate(Dataset = METADATAS[x])
      gc()
      DF 
    }
    stopCluster(cl)
    message(paste(i, "of", floor(nrow(env.data)/10), "ready", Sys.time()))
    Forms[[i]] <- Formulas
  }
  
  AllForms[[x]] <- Forms %>% 
    purrr::reduce(bind_rows) 
  print(paste(x, "of", length(METADATAS), "ready", Sys.time()))
  
  Dataset <- "JaccardDistance"
}


AllForms <- AllForms %>% 
  purrr::reduce(bind_rows) %>% 
  dplyr::distinct(Form, AICc, .keep_all = T) %>%
  dplyr::mutate(Max_VIF = NA)

NullMod <- data.frame(Form = paste(Dataset, "~ 1", collapse = ""), AICc = NA) %>% 
  mutate(Dataset = METADATAS[1]) %>%
  dplyr::mutate(Max_VIF = NA)

AllForms <- AllForms %>% 
  bind_rows(NullMod)
  

saveRDS(AllForms, "AllFormsBacterialPA.rds")
openxlsx::write.xlsx(AllForms, "AllFormsBacterialPA.xlsx")
```

</details>

This generate up to 5,061 models to evaluate, which can be downloaded as
an excel file
[here](https://github.com/Sustainscapes/AICcPermanova/raw/master/AllFormsBacterialPA.xlsx)
an rds
[here](https://github.com/Sustainscapes/AICcPermanova/blob/master/AllFormsBacterialPA.rds).

### 3.2.2 Model fitting

Then in the following code each model is fitted and AICc is calculated
to order it

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Model fitting code bacterial presence absence

</summary>

``` r
if(file.exists("logBacterialPA.txt")){
  file.remove("logBacterialPA.txt")
}

cl <- makeCluster(21)
registerDoParallel(cl)

Fs <- foreach(x = 1:nrow(AllForms), .packages = c("vegan", "dplyr", "tidyr", "readxl", "ampvis2", "car", "janitor"), .combine = bind_rows) %dopar% {
  AllForms <- readRDS("AllFormsBacterialPA.rds")
  meta.data = read_excel(AllForms$Dataset[x]) %>% 
    janitor::clean_names()
  
  d <- amp_load(
  otutable = "AC_otutale_new.txt",
  metadata = "Metadata_mix-samples_AC_Danielsen_final.xlsx")
  bac.data.subset = amp_subset_taxa(d, "d__Bacteria")
  bacterial_data = amp_subset_samples(d, Investigator == "AC")
  bacterial_data = as.data.frame(t(bacterial_data$abund)) %>% 
    janitor::clean_names() %>% 
    mutate_if(is.numeric, ~ifelse(.x > 0, 1, 0))
  # Remocing the first column (ID) in vegetation-dataset and cheking there are no NAs
  meta.data <- meta.data[match(rownames(bacterial_data), meta.data$seq_id),]
  JaccardDistance <- vegan::vegdist(bacterial_data, method = "jaccard")
  
  
  # Removing columns from env.data that is not used in the analysis 
  
  env.data = subset(meta.data, select = -c(order))
  
  
  
  env.data <- env.data %>% tidyr::drop_na()  
  
  Vars <- colnames(env.data)
  Dataset <- "JaccardDistance"
  
  Response = env.data
  

  gc()
  AICc.PERMANOVA2 <- function(adonis2.model) {
    
    # check to see if object is an adonis2 model...
    
    if (is.na(adonis2.model$SumOfSqs[1]))
      stop("object not output of adonis2 {vegan} ")
    
    # Ok, now extract appropriate terms from the adonis model Calculating AICc
    # using residual sum of squares (RSS or SSE) since I don't think that adonis
    # returns something I can use as a likelihood function... maximum likelihood
    # and MSE estimates are the same when distribution is gaussian See e.g.
    # https://www.jessicayung.com/mse-as-maximum-likelihood/;
    # https://towardsdatascience.com/probability-concepts-explained-maximum-likelihood-estimation-c7b4342fdbb1
    # So using RSS or MSE estimates is fine as long as the residuals are
    # Gaussian https://robjhyndman.com/hyndsight/aic/ If models have different
    # conditional likelihoods then AIC is not valid. However, comparing models
    # with different error distributions is ok (above link).
    
    
    RSS <- adonis2.model$SumOfSqs[ length(adonis2.model$SumOfSqs) - 1 ]
    MSE <- RSS / adonis2.model$Df[ length(adonis2.model$Df) - 1 ]
    
    nn <- adonis2.model$Df[ length(adonis2.model$Df) ] + 1
    
    k <- nn - adonis2.model$Df[ length(adonis2.model$Df) - 1 ]
    
    
    # AIC : 2*k + n*ln(RSS/n)
    # AICc: AIC + [2k(k+1)]/(n-k-1)
    
    # based on https://en.wikipedia.org/wiki/Akaike_information_criterion;
    # https://www.statisticshowto.datasciencecentral.com/akaikes-information-criterion/ ;
    # https://www.researchgate.net/post/What_is_the_AIC_formula;
    # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf;
    # https://medium.com/better-programming/data-science-modeling-how-to-use-linear-regression-with-python-fdf6ca5481be 
    
    
    AIC <- 2*k + nn*log(RSS/nn)
    AICc <- AIC + (2*k*(k + 1))/(nn - k - 1)
    
    output <- data.frame(AICc = AICc, k = k, N = nn)
    
    return(output)   
    
  }
  Temp <- AllForms[x,]
  Temp$AICc <-  try(AICc.PERMANOVA2(adonis2(as.formula(AllForms$Form[x]), data = Response, by = "margin"))$AICc, silent = T)
  
  Response$y <- rnorm(n = nrow(Response))
  
  VIF <- function(model) {
    tryCatch({
        vif <- car::vif(model)
        max(vif)
    }, error = function(e) {
        if (grepl("aliased coefficients", e$message)) {
            20000
        } else if (grepl("model contains fewer than 2 terms", e$message)) {
            0
        } else {
            stop(e)
        }
    })
}
  
  Temp$Max_VIF <- VIF(lm(as.formula(stringr::str_replace_all(AllForms$Form[x], "JaccardDistance", "y")), data = Response))
  
  Rs <- broom::tidy(adonis2(as.formula(AllForms$Form[x]), data = Response, by = "margin")) %>% dplyr::filter(!(term %in% c("Residual", "Total"))) %>% dplyr::select(term, R2) %>%  pivot_wider(names_from = term, values_from = R2)
  if((x %% 100) == 0){
    sink("logBacterialPA.txt", append = T)
    cat(paste("finished", x, "number of models", Sys.time(), "of",  nrow(AllForms)))
    cat("\n")
    sink()
  }
  
  bind_cols(Temp, Rs)
}

stopCluster(cl)

saveRDS(Fs, "FSBacterialPA.rds")

Fs <- Fs %>% arrange(AICc)

saveRDS(Fs, "FSBacterialPA.rds")
openxlsx::write.xlsx(Fs, "FSBacterialPA.xlsx")
```

</details>
