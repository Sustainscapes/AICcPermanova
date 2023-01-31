31/01, 2023

- <a href="#1-aiccpermanova" id="toc-1-aiccpermanova">1 AICcPermanova</a>
- <a href="#2-vegetation-models" id="toc-2-vegetation-models">2 Vegetation
  models</a>
  - <a href="#21-presence-absence-data"
    id="toc-21-presence-absence-data">2.1 Presence-absence data</a>
    - <a href="#211-generation-of-all-possible-models"
      id="toc-211-generation-of-all-possible-models">2.1.1 Generation of all
      possible models</a>
    - <a href="#212-model-fitting" id="toc-212-model-fitting">2.1.2 Model
      fitting</a>
  - <a href="#22-abundance-data" id="toc-22-abundance-data">2.2 Abundance
    data</a>
    - <a href="#221-generation-of-all-possible-models"
      id="toc-221-generation-of-all-possible-models">2.2.1 Generation of all
      possible models</a>
    - <a href="#222-model-fitting" id="toc-222-model-fitting">2.2.2 Model
      fitting</a>
- <a href="#3-baterial-models" id="toc-3-baterial-models">3 Baterial
  models</a>
  - <a href="#31-abundance-data" id="toc-31-abundance-data">3.1 Abundance
    data</a>
    - <a href="#311-generation-of-all-possible-models"
      id="toc-311-generation-of-all-possible-models">3.1.1 Generation of all
      possible models</a>
    - <a href="#312-model-fitting" id="toc-312-model-fitting">3.1.2 Model
      fitting</a>
  - <a href="#32-presence-absence-data"
    id="toc-32-presence-absence-data">3.2 Presence absence data</a>
    - <a href="#321-generation-of-all-possible-models"
      id="toc-321-generation-of-all-possible-models">3.2.1 Generation of all
      possible models</a>
    - <a href="#322-model-fitting" id="toc-322-model-fitting">3.2.2 Model
      fitting</a>

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
  env.data <- env.data %>% tidyr::drop_na() |> 
    dplyr::select(-habitat_type)
  
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
#> [1] "1 of 3 ready 2023-01-31 14:30:26"
#> [1] "2 of 3 ready 2023-01-31 14:31:33"
#> [1] "3 of 3 ready 2023-01-31 14:32:39"


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

This generate up to 2,358 models to evaluate, which can be downloaded as
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
  
  
  
  env.data <- env.data %>% tidyr::drop_na() |> 
    dplyr::select(-habitat_type)
  
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

As seen in table <a href="#tab:SummaryPlantPA">2.1</a> there are 1
models within 2 AICc where the max VIF is lower or equal than 6 of each
other, you can see there how many times a variable has been selected

| Variable | Number_of_models | Full_Akaike_Adjusted_RSq | Subset_Akaike_Adjusted_RSq |
|:---------|-----------------:|-------------------------:|---------------------------:|
| area     |                1 |                    0.579 |                      0.579 |

Table 2.1: Number of selected models were variables are present and
their Akaike Weighted R squared for the Marginal effect of the terms

Now we can see the top models that have a delta AICc within 2 of the
best model in table <a href="#tab:BestPLantModels">2.2</a> if expanded

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Show table of selected models for vegetation presence absence

</summary>

| Form                    |    AICc | DeltaAICc | Max_VIF |  area |
|:------------------------|--------:|----------:|--------:|------:|
| JaccardDistance \~ area | -55.246 |         0 |       0 | 0.579 |

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
  env.data <- env.data %>% tidyr::drop_na() |> 
    dplyr::select(-habitat_type)
  
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
#> [1] "1 of 3 ready 2023-01-31 14:39:30"
#> [1] "2 of 3 ready 2023-01-31 14:40:35"
#> [1] "3 of 3 ready 2023-01-31 14:41:40"


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

This generate up to 2,358 models to evaluate, which can be downloaded as
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
  
  
  
  env.data <- env.data %>% tidyr::drop_na() |> 
    dplyr::select(-habitat_type) 
  
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

As seen in table <a href="#tab:SummaryVegAbund">2.3</a> there are 1
models within 2 AICc of each other, you can see there how many times a
variable has been selected

| Variable | Number_of_models | Full_Akaike_Adjusted_RSq | Subset_Akaike_Adjusted_RSq |
|:---------|-----------------:|-------------------------:|---------------------------:|
| area     |                1 |                    0.681 |                      0.681 |

Table 2.3: Number of selected models were variables are present and
their Akaike Weighted R squared for the Marginal effect of the terms

Now we can see the top models that have a delta AICc within 2 of the
best model in table <a href="#tab:BestVegAbundModels">2.4</a> if
expanded

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Show table of selected models for vegetation presence absence

</summary>

| Form                 |    AICc | DeltaAICc | Max_VIF |  area |
|:---------------------|--------:|----------:|--------:|------:|
| BrayDistance \~ area | -72.174 |         0 |       0 | 0.681 |

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
  env.data <- env.data %>% tidyr::drop_na() |> 
    dplyr::select(-habitat_type)
  
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
#> [1] "1 of 3 ready 2023-01-31 14:46:05"
#> [1] "2 of 3 ready 2023-01-31 14:47:23"
#> [1] "3 of 3 ready 2023-01-31 14:48:45"


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

This generate up to 3,504 models to evaluate, which can be downloaded as
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
  
  
  
  env.data <- env.data %>% tidyr::drop_na()  |> 
    dplyr::select(-habitat_type)
  
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
96 models within 2 AICc of each other, you can see there how many times
a variable has been selected

| Variable          | Number_of_models | Full_Akaike_Adjusted_RSq | Subset_Akaike_Adjusted_RSq |
|:------------------|-----------------:|-------------------------:|---------------------------:|
| p_h\_water        |               96 |                    0.211 |                      0.211 |
| water_content     |               68 |                    0.033 |                      0.045 |
| shannon_veg       |               57 |                    0.021 |                      0.034 |
| oc_beregnet       |               52 |                    0.022 |                      0.040 |
| wr                |               27 |                    0.006 |                      0.023 |
| fine_silt         |               22 |                    0.006 |                      0.026 |
| ec                |               20 |                    0.003 |                      0.018 |
| coarse_silt_sand  |               17 |                    0.005 |                      0.027 |
| clay              |               16 |                    0.004 |                      0.022 |
| finesilt_and_clay |               15 |                    0.004 |                      0.027 |
| ammonium          |               10 |                    0.001 |                      0.014 |
| dexter_n          |                6 |                    0.001 |                      0.016 |
| nitrat_nitrit     |                3 |                    0.000 |                      0.012 |

Table 3.1: Number of selected models were variables are present and
their Akaike Weighted R squared for the Marginal effect of the terms

Now we can see the top models that have a delta AICc within 2 of the
best model in table <a href="#tab:BestBacterialAbundModels">3.2</a> if
expanded

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Show table of selected models for bacterial abundance

</summary>

| Form                                                                                       |    AICc | DeltaAICc | Max_VIF | p_h\_water |    ec | water_content |    wr | ammonium | nitrat_nitrit | oc_beregnet |  clay | fine_silt | coarse_silt_sand | shannon_veg | finesilt_and_clay | dexter_n |
|:-------------------------------------------------------------------------------------------|--------:|----------:|--------:|-----------:|------:|--------------:|------:|---------:|--------------:|------------:|------:|----------:|-----------------:|------------:|------------------:|---------:|
| BrayDistance \~ p_h\_water + water_content + oc_beregnet + fine_silt + shannon_veg         | -82.167 |     0.000 |   2.877 |      0.241 |    NA |         0.034 |    NA |       NA |            NA |       0.032 |    NA |     0.027 |               NA |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + oc_beregnet + shannon_veg                     | -82.119 |     0.048 |   2.341 |      0.270 |    NA |         0.034 |    NA |       NA |            NA |       0.033 |    NA |        NA |               NA |       0.035 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + oc_beregnet + coarse_silt_sand + shannon_veg  | -82.057 |     0.111 |   2.601 |      0.243 |    NA |         0.032 |    NA |       NA |            NA |       0.033 |    NA |        NA |            0.026 |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + oc_beregnet + finesilt_and_clay + shannon_veg | -82.057 |     0.111 |   2.601 |      0.243 |    NA |         0.032 |    NA |       NA |            NA |       0.033 |    NA |        NA |               NA |       0.034 |             0.026 |       NA |
| BrayDistance \~ p_h\_water + water_content + oc_beregnet + clay + shannon_veg              | -81.743 |     0.424 |   2.370 |      0.253 |    NA |         0.031 |    NA |       NA |            NA |       0.033 | 0.023 |        NA |               NA |       0.035 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + wr + shannon_veg                              | -81.664 |     0.503 |   2.658 |      0.164 |    NA |         0.059 | 0.028 |       NA |            NA |          NA |    NA |        NA |               NA |       0.036 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + fine_silt + shannon_veg                       | -81.649 |     0.518 |   1.717 |      0.246 |    NA |         0.061 |    NA |       NA |            NA |          NA |    NA |     0.028 |               NA |       0.033 |                NA |       NA |
| BrayDistance \~ p_h\_water + oc_beregnet + coarse_silt_sand + shannon_veg                  | -81.593 |     0.575 |   1.722 |      0.243 |    NA |            NA |    NA |       NA |            NA |       0.047 |    NA |        NA |            0.029 |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + oc_beregnet + finesilt_and_clay + shannon_veg                 | -81.593 |     0.575 |   1.514 |      0.243 |    NA |            NA |    NA |       NA |            NA |       0.060 |    NA |        NA |               NA |       0.034 |             0.029 |       NA |
| BrayDistance \~ p_h\_water + water_content + wr + coarse_silt_sand + shannon_veg           | -81.591 |     0.577 |   2.669 |      0.159 |    NA |         0.041 | 0.028 |       NA |            NA |          NA |    NA |        NA |            0.026 |       0.036 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + shannon_veg                                   | -81.568 |     0.599 |   1.055 |      0.301 |    NA |         0.067 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.033 |                NA |       NA |
| BrayDistance \~ p_h\_water + oc_beregnet + fine_silt + shannon_veg                         | -81.534 |     0.633 |   1.555 |      0.241 |    NA |            NA |    NA |       NA |            NA |       0.060 |    NA |     0.028 |               NA |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + wr + fine_silt + shannon_veg                  | -81.512 |     0.655 |   2.880 |      0.163 |    NA |         0.048 | 0.026 |       NA |            NA |          NA |    NA |     0.026 |               NA |       0.036 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + oc_beregnet + fine_silt                       | -81.508 |     0.659 |   2.857 |      0.249 |    NA |         0.033 |    NA |       NA |            NA |       0.031 |    NA |     0.028 |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + finesilt_and_clay + shannon_veg               | -81.502 |     0.665 |   1.546 |      0.252 |    NA |         0.059 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.033 |             0.027 |       NA |
| BrayDistance \~ p_h\_water + water_content + wr + finesilt_and_clay + shannon_veg          | -81.495 |     0.673 |   2.802 |      0.162 |    NA |         0.048 | 0.027 |       NA |            NA |          NA |    NA |        NA |               NA |       0.036 |             0.025 |       NA |
| BrayDistance \~ p_h\_water + water_content + coarse_silt_sand + shannon_veg                | -81.477 |     0.690 |   1.886 |      0.272 |    NA |         0.046 |    NA |       NA |            NA |          NA |    NA |        NA |            0.026 |       0.033 |                NA |       NA |
| BrayDistance \~ p_h\_water + oc_beregnet + shannon_veg                                     | -81.468 |     0.700 |   1.271 |      0.282 |    NA |            NA |    NA |       NA |            NA |       0.066 |    NA |        NA |               NA |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + oc_beregnet                                   | -81.451 |     0.716 |   2.325 |      0.278 |    NA |         0.034 |    NA |       NA |            NA |       0.032 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + oc_beregnet + coarse_silt_sand                | -81.391 |     0.776 |   2.588 |      0.252 |    NA |         0.032 |    NA |       NA |            NA |       0.032 |    NA |        NA |            0.027 |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + oc_beregnet + finesilt_and_clay               | -81.391 |     0.776 |   2.588 |      0.252 |    NA |         0.032 |    NA |       NA |            NA |       0.031 |    NA |        NA |               NA |          NA |             0.027 |       NA |
| BrayDistance \~ p_h\_water + oc_beregnet + clay + shannon_veg                              | -81.374 |     0.793 |   1.409 |      0.254 |    NA |            NA |    NA |       NA |            NA |       0.061 | 0.026 |        NA |               NA |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + oc_beregnet + shannon_veg                | -81.312 |     0.855 |   3.208 |      0.110 | 0.019 |         0.035 |    NA |       NA |            NA |       0.032 |    NA |        NA |               NA |       0.035 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + wr + clay + shannon_veg                       | -81.268 |     0.900 |   2.704 |      0.161 |    NA |         0.050 | 0.028 |       NA |            NA |          NA | 0.023 |        NA |               NA |       0.036 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + wr + oc_beregnet + shannon_veg                | -81.228 |     0.940 |   3.610 |      0.161 |    NA |         0.034 | 0.018 |       NA |            NA |       0.023 |    NA |        NA |               NA |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + clay + shannon_veg                            | -81.171 |     0.996 |   1.252 |      0.270 |    NA |         0.059 |    NA |       NA |            NA |          NA | 0.023 |        NA |               NA |       0.033 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + fine_silt                                     | -81.169 |     0.998 |   1.706 |      0.256 |    NA |         0.058 |    NA |       NA |            NA |          NA |    NA |     0.029 |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content                                                 | -81.083 |     1.084 |   1.010 |      0.312 |    NA |         0.063 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + oc_beregnet + clay                            | -81.055 |     1.112 |   2.353 |      0.261 |    NA |         0.031 |    NA |       NA |            NA |       0.032 | 0.023 |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + finesilt_and_clay                             | -81.033 |     1.134 |   1.542 |      0.262 |    NA |         0.056 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.027 |       NA |
| BrayDistance \~ p_h\_water + oc_beregnet + coarse_silt_sand                                | -81.015 |     1.153 |   1.685 |      0.251 |    NA |            NA |    NA |       NA |            NA |       0.045 |    NA |        NA |            0.029 |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + oc_beregnet + finesilt_and_clay                               | -81.015 |     1.153 |   1.453 |      0.251 |    NA |            NA |    NA |       NA |            NA |       0.056 |    NA |        NA |               NA |          NA |             0.029 |       NA |
| BrayDistance \~ p_h\_water + water_content + coarse_silt_sand                              | -80.997 |     1.170 |   1.885 |      0.282 |    NA |         0.045 |    NA |       NA |            NA |          NA |    NA |        NA |            0.027 |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + wr + shannon_veg                         | -80.984 |     1.183 |   4.606 |      0.087 | 0.020 |         0.055 | 0.028 |       NA |            NA |          NA |    NA |        NA |               NA |       0.036 |                NA |       NA |
| BrayDistance \~ p_h\_water + oc_beregnet + fine_silt                                       | -80.978 |     1.189 |   1.491 |      0.249 |    NA |            NA |    NA |       NA |            NA |       0.056 |    NA |     0.028 |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + oc_beregnet                                                   | -80.931 |     1.237 |   1.211 |      0.293 |    NA |            NA |    NA |       NA |            NA |       0.061 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + wr                                            | -80.921 |     1.246 |   2.381 |      0.179 |    NA |         0.056 | 0.026 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + shannon_veg                              | -80.913 |     1.254 |   2.269 |      0.126 | 0.020 |         0.058 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.033 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + fine_silt + shannon_veg                  | -80.898 |     1.270 |   2.588 |      0.129 | 0.019 |         0.058 |    NA |       NA |            NA |          NA |    NA |     0.027 |               NA |       0.033 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + clay + fine_silt + shannon_veg                | -80.832 |     1.335 |   4.127 |      0.247 |    NA |         0.064 |    NA |       NA |            NA |          NA | 0.019 |     0.024 |               NA |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + wr + coarse_silt_sand                         | -80.808 |     1.359 |   2.393 |      0.174 |    NA |         0.040 | 0.025 |       NA |            NA |          NA |    NA |        NA |            0.026 |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + oc_beregnet + clay                                            | -80.779 |     1.388 |   1.353 |      0.263 |    NA |            NA |    NA |       NA |            NA |       0.057 | 0.026 |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + ammonium + oc_beregnet + shannon_veg          | -80.752 |     1.415 |   2.894 |      0.264 |    NA |         0.033 |    NA |    0.013 |            NA |       0.031 |    NA |        NA |               NA |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + wr + fine_silt                                | -80.723 |     1.444 |   2.595 |      0.178 |    NA |         0.047 | 0.023 |       NA |            NA |          NA |    NA |     0.025 |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + wr + finesilt_and_clay                        | -80.693 |     1.474 |   2.515 |      0.177 |    NA |         0.046 | 0.024 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.025 |       NA |
| BrayDistance \~ p_h\_water + water_content + clay                                          | -80.688 |     1.479 |   1.252 |      0.280 |    NA |         0.056 |    NA |       NA |            NA |          NA | 0.023 |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + finesilt_and_clay + shannon_veg          | -80.671 |     1.496 |   2.567 |      0.128 | 0.019 |         0.056 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.033 |             0.025 |       NA |
| BrayDistance \~ p_h\_water + water_content + nitrat_nitrit + oc_beregnet + shannon_veg     | -80.664 |     1.503 |   2.353 |      0.239 |    NA |         0.033 |    NA |       NA |         0.012 |       0.032 |    NA |        NA |               NA |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + wr + oc_beregnet                              | -80.645 |     1.522 |   3.225 |      0.178 |    NA |         0.033 | 0.018 |       NA |            NA |       0.025 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + oc_beregnet                              | -80.639 |     1.528 |   3.146 |      0.110 | 0.018 |         0.034 |    NA |       NA |            NA |       0.030 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + wr + oc_beregnet + coarse_silt_sand + shannon_veg             | -80.627 |     1.540 |   3.639 |      0.158 |    NA |            NA | 0.017 |       NA |            NA |       0.031 |    NA |        NA |            0.028 |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + wr + oc_beregnet + finesilt_and_clay + shannon_veg            | -80.627 |     1.540 |   3.639 |      0.158 |    NA |            NA | 0.017 |       NA |            NA |       0.039 |    NA |        NA |               NA |       0.034 |             0.028 |       NA |
| BrayDistance \~ p_h\_water + wr + oc_beregnet + shannon_veg                                | -80.622 |     1.546 |   3.436 |      0.161 |    NA |            NA | 0.018 |       NA |            NA |       0.047 |    NA |        NA |               NA |       0.033 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + oc_beregnet + fine_silt + shannon_veg                    | -80.604 |     1.563 |   3.277 |      0.108 | 0.017 |            NA |    NA |       NA |            NA |       0.055 |    NA |     0.028 |               NA |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + oc_beregnet + shannon_veg                                | -80.598 |     1.569 |   3.187 |      0.110 | 0.018 |            NA |    NA |       NA |            NA |       0.055 |    NA |        NA |               NA |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + oc_beregnet + coarse_silt_sand + shannon_veg             | -80.594 |     1.573 |   3.215 |      0.108 | 0.017 |            NA |    NA |       NA |            NA |       0.048 |    NA |        NA |            0.028 |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + oc_beregnet + finesilt_and_clay + shannon_veg            | -80.594 |     1.573 |   3.215 |      0.108 | 0.017 |            NA |    NA |       NA |            NA |       0.055 |    NA |        NA |               NA |       0.034 |             0.028 |       NA |
| BrayDistance \~ p_h\_water + water_content + oc_beregnet + clay + fine_silt                | -80.549 |     1.618 |   4.362 |      0.250 |    NA |         0.035 |    NA |       NA |            NA |       0.031 | 0.017 |     0.022 |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + wr + dexter_n + shannon_veg                   | -80.513 |     1.654 |   2.721 |      0.164 |    NA |         0.059 | 0.028 |       NA |            NA |          NA |    NA |        NA |               NA |       0.033 |                NA |    0.015 |
| BrayDistance \~ p_h\_water + water_content + wr + oc_beregnet + fine_silt                  | -80.503 |     1.664 |   3.349 |      0.176 |    NA |         0.033 | 0.017 |       NA |            NA |       0.025 |    NA |     0.026 |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + wr + oc_beregnet + fine_silt + shannon_veg                    | -80.494 |     1.674 |   3.711 |      0.159 |    NA |            NA | 0.016 |       NA |            NA |       0.038 |    NA |     0.026 |               NA |       0.033 |                NA |       NA |
| BrayDistance \~ p_h\_water + oc_beregnet + clay + fine_silt + shannon_veg                  | -80.483 |     1.684 |   3.226 |      0.242 |    NA |            NA |    NA |       NA |            NA |       0.060 | 0.016 |     0.018 |               NA |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + wr + oc_beregnet + clay + shannon_veg                         | -80.478 |     1.689 |   3.520 |      0.158 |    NA |            NA | 0.018 |       NA |            NA |       0.042 | 0.026 |        NA |               NA |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + dexter_n + shannon_veg                        | -80.477 |     1.691 |   1.104 |      0.293 |    NA |         0.066 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.029 |                NA |    0.015 |
| BrayDistance \~ p_h\_water + water_content + wr + clay                                     | -80.462 |     1.705 |   2.420 |      0.176 |    NA |         0.048 | 0.025 |       NA |            NA |          NA | 0.023 |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + wr + oc_beregnet + coarse_silt_sand           | -80.449 |     1.718 |   3.302 |      0.175 |    NA |         0.032 | 0.017 |       NA |            NA |       0.024 |    NA |        NA |            0.026 |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + wr + oc_beregnet + finesilt_and_clay          | -80.449 |     1.718 |   3.302 |      0.175 |    NA |         0.032 | 0.017 |       NA |            NA |       0.025 |    NA |        NA |               NA |          NA |             0.026 |       NA |
| BrayDistance \~ p_h\_water + water_content + ammonium + shannon_veg                        | -80.446 |     1.721 |   2.434 |      0.285 |    NA |         0.039 |    NA |    0.015 |            NA |          NA |    NA |        NA |               NA |       0.032 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + coarse_silt_sand + shannon_veg           | -80.445 |     1.722 |   2.886 |      0.122 | 0.016 |         0.046 |    NA |       NA |            NA |          NA |    NA |        NA |            0.023 |       0.033 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + fine_silt + dexter_n + shannon_veg            | -80.443 |     1.724 |   1.956 |      0.246 |    NA |         0.059 |    NA |       NA |            NA |          NA |    NA |     0.027 |               NA |       0.029 |                NA |    0.015 |
| BrayDistance \~ p_h\_water + ec + water_content + fine_silt                                | -80.442 |     1.725 |   2.583 |      0.131 | 0.019 |         0.056 |    NA |       NA |            NA |          NA |    NA |     0.028 |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content                                            | -80.435 |     1.733 |   2.269 |      0.127 | 0.020 |         0.055 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + ammonium + oc_beregnet + coarse_silt_sand + shannon_veg       | -80.380 |     1.788 |   2.144 |      0.242 |    NA |            NA |    NA |    0.014 |            NA |       0.031 |    NA |        NA |            0.029 |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + ammonium + oc_beregnet + finesilt_and_clay + shannon_veg      | -80.380 |     1.788 |   1.848 |      0.242 |    NA |            NA |    NA |    0.014 |            NA |       0.035 |    NA |        NA |               NA |       0.034 |             0.029 |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + oc_beregnet + fine_silt                  | -80.358 |     1.810 |   3.436 |      0.105 | 0.015 |         0.031 |    NA |       NA |            NA |       0.027 |    NA |     0.025 |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + wr + coarse_silt_sand + shannon_veg                           | -80.347 |     1.820 |   2.652 |      0.159 |    NA |            NA | 0.034 |       NA |            NA |          NA |    NA |        NA |            0.044 |       0.035 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + dexter_n                                      | -80.345 |     1.822 |   1.041 |      0.307 |    NA |         0.063 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.019 |
| BrayDistance \~ p_h\_water + water_content + wr + ammonium + shannon_veg                   | -80.335 |     1.833 |   2.706 |      0.164 |    NA |         0.038 | 0.027 |    0.013 |            NA |          NA |    NA |        NA |               NA |       0.036 |                NA |       NA |
| BrayDistance \~ p_h\_water + ammonium + oc_beregnet + fine_silt + shannon_veg              | -80.318 |     1.849 |   1.853 |      0.240 |    NA |            NA |    NA |    0.014 |            NA |       0.035 |    NA |     0.029 |               NA |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + fine_silt + dexter_n                          | -80.295 |     1.872 |   1.956 |      0.256 |    NA |         0.057 |    NA |       NA |            NA |          NA |    NA |     0.027 |               NA |          NA |                NA |    0.018 |
| BrayDistance \~ p_h\_water + ec + water_content + clay + shannon_veg                       | -80.291 |     1.877 |   2.463 |      0.126 | 0.018 |         0.056 |    NA |       NA |            NA |          NA | 0.021 |        NA |               NA |       0.033 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + clay + fine_silt                              | -80.281 |     1.886 |   4.079 |      0.257 |    NA |         0.060 |    NA |       NA |            NA |          NA | 0.018 |     0.023 |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + oc_beregnet + clay + shannon_veg                         | -80.278 |     1.889 |   3.205 |      0.109 | 0.016 |            NA |    NA |       NA |            NA |       0.055 | 0.024 |        NA |               NA |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + coarse_silt_sand + dexter_n + shannon_veg     | -80.260 |     1.907 |   2.191 |      0.271 |    NA |         0.043 |    NA |       NA |            NA |          NA |    NA |        NA |            0.025 |       0.030 |                NA |    0.014 |
| BrayDistance \~ p_h\_water + ammonium + oc_beregnet + shannon_veg                          | -80.232 |     1.935 |   1.833 |      0.282 |    NA |            NA |    NA |    0.014 |            NA |       0.036 |    NA |        NA |               NA |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + finesilt_and_clay                        | -80.231 |     1.936 |   2.564 |      0.129 | 0.019 |         0.054 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.026 |       NA |
| BrayDistance \~ p_h\_water + water_content + nitrat_nitrit + shannon_veg                   | -80.230 |     1.937 |   1.433 |      0.251 |    NA |         0.063 |    NA |       NA |         0.013 |          NA |    NA |        NA |               NA |       0.033 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + ammonium + fine_silt + shannon_veg            | -80.226 |     1.942 |   3.194 |      0.244 |    NA |         0.034 |    NA |    0.012 |            NA |          NA |    NA |     0.025 |               NA |       0.032 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + ammonium + coarse_silt_sand + shannon_veg     | -80.222 |     1.945 |   3.363 |      0.263 |    NA |         0.030 |    NA |    0.014 |            NA |          NA |    NA |        NA |            0.025 |       0.033 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + wr                                       | -80.217 |     1.950 |   4.161 |      0.091 | 0.020 |         0.052 | 0.025 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + ammonium + oc_beregnet                        | -80.217 |     1.950 |   2.874 |      0.275 |    NA |         0.034 |    NA |    0.014 |            NA |       0.030 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + wr + oc_beregnet                                              | -80.216 |     1.952 |   3.109 |      0.178 |    NA |            NA | 0.019 |       NA |            NA |       0.048 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + nitrat_nitrit + oc_beregnet + shannon_veg                     | -80.197 |     1.970 |   1.562 |      0.239 |    NA |            NA |    NA |       NA |         0.013 |       0.063 |    NA |        NA |               NA |       0.034 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + oc_beregnet + coarse_silt_sand           | -80.190 |     1.977 |   3.292 |      0.106 | 0.015 |         0.030 |    NA |       NA |            NA |       0.030 |    NA |        NA |            0.023 |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + oc_beregnet + finesilt_and_clay          | -80.190 |     1.977 |   3.292 |      0.106 | 0.015 |         0.030 |    NA |       NA |            NA |       0.027 |    NA |        NA |               NA |          NA |             0.023 |       NA |
| BrayDistance \~ p_h\_water + water_content + wr + oc_beregnet + clay                       | -80.179 |     1.988 |   3.247 |      0.174 |    NA |         0.030 | 0.018 |       NA |            NA |       0.025 | 0.023 |        NA |               NA |          NA |                NA |       NA |

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
  env.data <- env.data %>% tidyr::drop_na() |> 
    dplyr::select(-habitat_type) 
  
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

This generate up to 3,504 models to evaluate, which can be downloaded as
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
  
  
  
  env.data <- env.data %>% tidyr::drop_na() |> 
    dplyr::select(-habitat_type)  
  
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
