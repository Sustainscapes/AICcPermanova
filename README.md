23/03, 2023

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
  Vars <- Vars[Vars != "area"]
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
#> [1] "1 of 3 ready 2023-03-23 11:36:08"
#> [1] "2 of 3 ready 2023-03-23 11:37:12"
#> [1] "3 of 3 ready 2023-03-23 11:38:14"


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

This generate up to 1,536 models to evaluate, which can be downloaded as
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
  Vars <- Vars[Vars != "area"]
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
  Temp$AICc <-  try(AICc.PERMANOVA2(with(Response, adonis2(as.formula(AllForms$Form[x]), data = Response, by = "margin", strata = area)))$AICc, silent = T)
  
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
  
  Rs <- broom::tidy(with(Response,adonis2(as.formula(AllForms$Form[x]), data = Response, by = "margin"))) %>% dplyr::filter(!(term %in% c("Residual", "Total"))) %>% dplyr::select(term, R2) %>%  pivot_wider(names_from = term, values_from = R2)
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

As seen in table <a href="#tab:SummaryPlantPA">2.1</a> there are 191
models within 2 AICc where the max VIF is lower or equal than 6 of each
other, you can see there how many times a variable has been selected

| Variable          | Number_of_models | Full_Akaike_Adjusted_RSq | Subset_Akaike_Adjusted_RSq |
|:------------------|-----------------:|-------------------------:|---------------------------:|
| p_h\_water        |              121 |                    0.051 |                      0.066 |
| ec                |              120 |                    0.030 |                      0.050 |
| water_content     |              108 |                    0.021 |                      0.043 |
| fine_silt         |               49 |                    0.009 |                      0.035 |
| wr                |               47 |                    0.006 |                      0.028 |
| oc_beregnet       |               46 |                    0.007 |                      0.031 |
| ammonium          |               44 |                    0.005 |                      0.027 |
| shannon_bak       |               42 |                    0.018 |                      0.041 |
| finesilt_and_clay |               33 |                    0.006 |                      0.034 |
| coarse_silt_sand  |               32 |                    0.006 |                      0.031 |
| clay              |               31 |                    0.004 |                      0.027 |
| nitrat_nitrit     |               13 |                    0.001 |                      0.018 |
| dexter_n          |                9 |                    0.001 |                      0.019 |

Table 2.1: Number of selected models were variables are present and
their Akaike Weighted R squared for the Marginal effect of the terms

Now we can see the top models that have a delta AICc within 2 of the
best model in table <a href="#tab:BestPLantModels">2.2</a> if expanded

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Show table of selected models for vegetation presence absence

</summary>

| Form                                                                                          |    AICc | DeltaAICc | Max_VIF | p_h\_water |    ec | water_content |    wr | ammonium | nitrat_nitrit | oc_beregnet |  clay | fine_silt | coarse_silt_sand | shannon_bak | finesilt_and_clay | dexter_n |
|:----------------------------------------------------------------------------------------------|--------:|----------:|--------:|-----------:|------:|--------------:|------:|---------:|--------------:|------------:|------:|----------:|-----------------:|------------:|------------------:|---------:|
| JaccardDistance \~ p_h\_water + ec + water_content                                            | -52.326 |     0.000 |   2.269 |      0.056 | 0.045 |         0.048 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + shannon_bak                                              | -52.312 |     0.014 |   5.976 |      0.060 | 0.047 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.048 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + fine_silt + shannon_bak                                  | -52.133 |     0.194 |   6.375 |      0.050 | 0.045 |            NA |    NA |       NA |            NA |          NA |    NA |     0.034 |               NA |       0.039 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + oc_beregnet                              | -52.069 |     0.258 |   3.146 |      0.050 | 0.040 |         0.058 |    NA |       NA |            NA |       0.033 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + fine_silt                                | -52.058 |     0.269 |   2.583 |      0.056 | 0.043 |         0.038 |    NA |       NA |            NA |          NA |    NA |     0.033 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + water_content + shannon_bak                                           | -52.055 |     0.272 |   2.143 |         NA | 0.060 |         0.056 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.052 |                NA |       NA |
| JaccardDistance \~ p_h\_water + fine_silt                                                     | -52.004 |     0.322 |   1.161 |      0.105 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.043 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + fine_silt                                                | -51.996 |     0.330 |   2.470 |      0.062 | 0.037 |            NA |    NA |       NA |            NA |          NA |    NA |     0.043 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + finesilt_and_clay                        | -51.984 |     0.342 |   2.564 |      0.056 | 0.043 |         0.041 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.032 |       NA |
| JaccardDistance \~ p_h\_water + ec + finesilt_and_clay + shannon_bak                          | -51.941 |     0.385 |   6.351 |      0.051 | 0.045 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.040 |             0.031 |       NA |
| JaccardDistance \~ p_h\_water + water_content + oc_beregnet                                   | -51.888 |     0.438 |   2.325 |      0.108 |    NA |         0.056 |    NA |       NA |            NA |       0.038 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content                                                 | -51.863 |     0.464 |   1.010 |      0.112 |    NA |         0.041 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + finesilt_and_clay                                             | -51.784 |     0.543 |   1.140 |      0.106 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.039 |       NA |
| JaccardDistance \~ p_h\_water + ec + wr + shannon_bak                                         | -51.772 |     0.554 |   8.770 |      0.065 | 0.046 |            NA | 0.029 |       NA |            NA |          NA |    NA |        NA |               NA |       0.053 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + finesilt_and_clay                                        | -51.741 |     0.585 |   2.405 |      0.062 | 0.037 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.039 |       NA |
| JaccardDistance \~ p_h\_water + water_content + oc_beregnet + fine_silt                       | -51.733 |     0.594 |   2.857 |      0.093 |    NA |         0.048 |    NA |       NA |            NA |       0.038 |    NA |     0.035 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + clay                                     | -51.717 |     0.609 |   2.462 |      0.056 | 0.043 |         0.045 |    NA |       NA |            NA |          NA | 0.028 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + wr                                       | -51.715 |     0.611 |   4.161 |      0.043 | 0.043 |         0.053 | 0.028 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + water_content + fine_silt + shannon_bak                               | -51.709 |     0.617 |   2.463 |         NA | 0.057 |         0.044 |    NA |       NA |            NA |          NA |    NA |     0.032 |               NA |       0.051 |                NA |       NA |
| JaccardDistance \~ p_h\_water + shannon_bak                                                   | -51.707 |     0.620 |   5.305 |      0.064 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.038 |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + fine_silt                                     | -51.702 |     0.624 |   1.706 |      0.089 |    NA |         0.033 |    NA |       NA |            NA |          NA |    NA |     0.035 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water                                                                 | -51.678 |     0.648 |   0.000 |      0.115 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + coarse_silt_sand                         | -51.662 |     0.665 |   2.883 |      0.055 | 0.039 |         0.039 |    NA |       NA |            NA |          NA |    NA |        NA |            0.027 |          NA |                NA |       NA |
| JaccardDistance \~ ec + water_content + finesilt_and_clay + shannon_bak                       | -51.649 |     0.677 |   2.420 |         NA | 0.058 |         0.047 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.051 |             0.031 |       NA |
| JaccardDistance \~ p_h\_water + ec                                                            | -51.634 |     0.692 |   1.841 |      0.063 | 0.037 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + fine_silt + shannon_bak                                       | -51.619 |     0.708 |   6.160 |      0.061 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.036 |               NA |       0.031 |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + finesilt_and_clay                             | -51.615 |     0.711 |   1.542 |      0.093 |    NA |         0.035 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.033 |       NA |
| JaccardDistance \~ p_h\_water + ec + clay + shannon_bak                                       | -51.594 |     0.733 |   6.217 |      0.054 | 0.045 |            NA |    NA |       NA |            NA |          NA | 0.026 |        NA |               NA |       0.044 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + shannon_bak                              | -51.592 |     0.735 |   8.200 |      0.030 | 0.049 |         0.026 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.026 |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + oc_beregnet + coarse_silt_sand                | -51.581 |     0.745 |   2.588 |      0.096 |    NA |         0.049 |    NA |       NA |            NA |       0.038 |    NA |        NA |            0.033 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + oc_beregnet + finesilt_and_clay               | -51.581 |     0.745 |   2.588 |      0.096 |    NA |         0.049 |    NA |       NA |            NA |       0.037 |    NA |        NA |               NA |          NA |             0.033 |       NA |
| JaccardDistance \~ p_h\_water + ec + coarse_silt_sand                                         | -51.552 |     0.774 |   2.875 |      0.056 | 0.038 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.036 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + coarse_silt_sand + shannon_bak                           | -51.552 |     0.774 |   6.737 |      0.042 | 0.042 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.026 |       0.037 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + ammonium                                 | -51.548 |     0.778 |   2.510 |      0.056 | 0.045 |         0.044 |    NA |    0.025 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + coarse_silt_sand                              | -51.540 |     0.787 |   1.885 |      0.099 |    NA |         0.038 |    NA |       NA |            NA |          NA |    NA |        NA |            0.032 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + coarse_silt_sand                                              | -51.528 |     0.798 |   1.030 |      0.112 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.035 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + ammonium + fine_silt                                     | -51.492 |     0.834 |   2.551 |      0.059 | 0.038 |            NA |    NA |    0.030 |            NA |          NA |    NA |     0.043 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + finesilt_and_clay + shannon_bak                               | -51.450 |     0.876 |   6.108 |      0.062 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.032 |             0.033 |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + oc_beregnet + fine_silt                  | -51.448 |     0.878 |   3.436 |      0.049 | 0.033 |         0.042 |    NA |       NA |            NA |       0.028 |    NA |     0.028 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ammonium + fine_silt                                          | -51.442 |     0.885 |   1.277 |      0.098 |    NA |            NA |    NA |    0.028 |            NA |          NA |    NA |     0.043 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + ammonium + shannon_bak                                   | -51.438 |     0.889 |   6.496 |      0.049 | 0.047 |            NA |    NA |    0.024 |            NA |          NA |    NA |        NA |               NA |       0.042 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + oc_beregnet + shannon_bak                                | -51.432 |     0.894 |   8.520 |      0.052 | 0.042 |            NA |    NA |       NA |            NA |       0.024 |    NA |        NA |               NA |       0.048 |                NA |       NA |
| JaccardDistance \~ ec + water_content + coarse_silt_sand + shannon_bak                        | -51.417 |     0.909 |   2.650 |         NA | 0.054 |         0.040 |    NA |       NA |            NA |          NA |    NA |        NA |            0.028 |       0.051 |                NA |       NA |
| JaccardDistance \~ ec + water_content + clay + shannon_bak                                    | -51.417 |     0.910 |   2.301 |         NA | 0.059 |         0.052 |    NA |       NA |            NA |          NA | 0.028 |        NA |               NA |       0.052 |                NA |       NA |
| JaccardDistance \~ ec + water_content + wr                                                    | -51.372 |     0.954 |   1.310 |         NA | 0.073 |         0.058 | 0.041 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + wr                                            | -51.339 |     0.987 |   2.381 |      0.072 |    NA |         0.050 | 0.029 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + wr + fine_silt                           | -51.334 |     0.992 |   4.336 |      0.041 | 0.039 |         0.045 | 0.027 |       NA |            NA |          NA |    NA |     0.032 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + clay                                          | -51.328 |     0.998 |   1.252 |      0.102 |    NA |         0.038 |    NA |       NA |            NA |          NA | 0.029 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + water_content + oc_beregnet + shannon_bak                             | -51.320 |     1.006 |   2.671 |         NA | 0.056 |         0.050 |    NA |       NA |            NA |       0.026 |    NA |        NA |               NA |       0.039 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + oc_beregnet + coarse_silt_sand           | -51.316 |     1.010 |   3.292 |      0.048 | 0.033 |         0.045 |    NA |       NA |            NA |       0.032 |    NA |        NA |            0.026 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + oc_beregnet + finesilt_and_clay          | -51.316 |     1.010 |   3.292 |      0.048 | 0.033 |         0.045 |    NA |       NA |            NA |       0.027 |    NA |        NA |               NA |          NA |             0.026 |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + oc_beregnet + shannon_bak                | -51.315 |     1.012 |   9.240 |      0.037 | 0.043 |         0.036 |    NA |       NA |            NA |       0.033 |    NA |        NA |               NA |       0.026 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + ammonium + fine_silt + shannon_bak                       | -51.309 |     1.017 |   6.796 |      0.043 | 0.045 |            NA |    NA |    0.025 |            NA |          NA |    NA |     0.036 |               NA |       0.035 |                NA |       NA |
| JaccardDistance \~ p_h\_water + clay                                                          | -51.299 |     1.028 |   1.083 |      0.109 |    NA |            NA |    NA |       NA |            NA |          NA | 0.031 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + clay + coarse_silt_sand                       | -51.285 |     1.041 |   9.834 |      0.101 |    NA |         0.042 |    NA |       NA |            NA |          NA | 0.034 |        NA |            0.037 |          NA |                NA |       NA |
| JaccardDistance \~ ec + fine_silt + shannon_bak                                               | -51.281 |     1.045 |   2.463 |         NA | 0.056 |            NA |    NA |       NA |            NA |          NA |    NA |     0.044 |               NA |       0.051 |                NA |       NA |
| JaccardDistance \~ p_h\_water + coarse_silt_sand + shannon_bak                                | -51.279 |     1.047 |   6.685 |      0.061 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.031 |       0.034 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + ammonium + finesilt_and_clay                             | -51.273 |     1.054 |   2.491 |      0.059 | 0.038 |            NA |    NA |    0.030 |            NA |          NA |    NA |        NA |               NA |          NA |             0.040 |       NA |
| JaccardDistance \~ p_h\_water + ammonium + finesilt_and_clay                                  | -51.271 |     1.055 |   1.249 |      0.100 |    NA |            NA |    NA |    0.029 |            NA |          NA |    NA |        NA |               NA |          NA |             0.041 |       NA |
| JaccardDistance \~ p_h\_water + water_content + oc_beregnet + clay                            | -51.271 |     1.055 |   2.353 |      0.102 |    NA |         0.051 |    NA |       NA |            NA |       0.037 | 0.028 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + fine_silt + shannon_bak                  | -51.251 |     1.076 |   8.209 |      0.031 | 0.046 |         0.024 |    NA |       NA |            NA |          NA |    NA |     0.032 |               NA |       0.025 |                NA |       NA |
| JaccardDistance \~ ec + coarse_silt_sand + shannon_bak                                        | -51.241 |     1.085 |   2.606 |         NA | 0.060 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.043 |       0.051 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + wr + finesilt_and_clay                   | -51.233 |     1.093 |   4.294 |      0.041 | 0.040 |         0.047 | 0.026 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.030 |       NA |
| JaccardDistance \~ p_h\_water + wr + shannon_bak                                              | -51.230 |     1.097 |   8.286 |      0.075 |    NA |            NA | 0.030 |       NA |            NA |          NA |    NA |        NA |               NA |       0.048 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + dexter_n                                 | -51.222 |     1.105 |   2.280 |      0.056 | 0.044 |         0.049 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.020 |
| JaccardDistance \~ ec + water_content + oc_beregnet                                           | -51.221 |     1.105 |   1.944 |         NA | 0.097 |         0.056 |    NA |       NA |            NA |       0.039 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + water_content + ammonium + shannon_bak                                | -51.213 |     1.114 |   2.433 |         NA | 0.061 |         0.046 |    NA |    0.025 |            NA |          NA |    NA |        NA |               NA |       0.051 |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + wr + fine_silt                                | -51.198 |     1.128 |   2.595 |      0.070 |    NA |         0.044 | 0.030 |       NA |            NA |          NA |    NA |     0.035 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + clay                                                     | -51.198 |     1.129 |   2.180 |      0.062 | 0.036 |            NA |    NA |       NA |            NA |          NA | 0.031 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + clay + fine_silt + shannon_bak                           | -51.178 |     1.148 |   6.376 |      0.050 | 0.045 |            NA |    NA |       NA |            NA |          NA | 0.023 |     0.031 |               NA |       0.039 |                NA |       NA |
| JaccardDistance \~ ec + water_content                                                         | -51.167 |     1.159 |   1.067 |         NA | 0.101 |         0.055 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + dexter_n + shannon_bak                                   | -51.165 |     1.161 |   6.062 |      0.060 | 0.048 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.048 |                NA |    0.020 |
| JaccardDistance \~ p_h\_water + ec + nitrat_nitrit + shannon_bak                              | -51.154 |     1.173 |   6.142 |      0.057 | 0.044 |            NA |    NA |       NA |         0.020 |          NA |    NA |        NA |               NA |       0.045 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + finesilt_and_clay + shannon_bak          | -51.147 |     1.179 |   8.201 |      0.030 | 0.046 |         0.026 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.025 |             0.031 |       NA |
| JaccardDistance \~ p_h\_water + ec + ammonium                                                 | -51.146 |     1.181 |   2.092 |      0.059 | 0.039 |            NA |    NA |    0.030 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + oc_beregnet + clay                       | -51.140 |     1.186 |   3.182 |      0.049 | 0.036 |         0.051 |    NA |       NA |            NA |       0.029 | 0.024 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + water_content + wr + shannon_bak                                      | -51.135 |     1.191 |   2.941 |         NA | 0.061 |         0.056 | 0.023 |       NA |            NA |          NA |    NA |        NA |               NA |       0.034 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + ammonium + finesilt_and_clay + shannon_bak               | -51.135 |     1.192 |   6.774 |      0.043 | 0.044 |            NA |    NA |    0.026 |            NA |          NA |    NA |        NA |               NA |       0.035 |             0.033 |       NA |
| JaccardDistance \~ p_h\_water + oc_beregnet + shannon_bak                                     | -51.123 |     1.203 |   7.448 |      0.068 |    NA |            NA |    NA |       NA |            NA |       0.028 |    NA |        NA |               NA |       0.044 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + oc_beregnet + fine_silt + shannon_bak                    | -51.116 |     1.211 |   9.134 |      0.042 | 0.039 |            NA |    NA |       NA |            NA |       0.022 |    NA |     0.033 |               NA |       0.038 |                NA |       NA |
| JaccardDistance \~ p_h\_water + clay + fine_silt                                              | -51.107 |     1.219 |   3.129 |      0.105 |    NA |            NA |    NA |       NA |            NA |          NA | 0.023 |     0.035 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + oc_beregnet + fine_silt                                  | -51.104 |     1.222 |   3.277 |      0.048 | 0.039 |            NA |    NA |       NA |            NA |       0.024 |    NA |     0.043 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + clay + shannon_bak                                            | -51.099 |     1.227 |   5.841 |      0.063 |    NA |            NA |    NA |       NA |            NA |          NA | 0.028 |        NA |               NA |       0.035 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + clay + coarse_silt_sand                                  | -51.097 |     1.229 |   6.013 |      0.053 | 0.040 |            NA |    NA |       NA |            NA |          NA | 0.031 |        NA |            0.036 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + ammonium + coarse_silt_sand                              | -51.094 |     1.232 |   2.880 |      0.055 | 0.038 |            NA |    NA |    0.031 |            NA |          NA |    NA |        NA |            0.037 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ammonium                                                      | -51.089 |     1.237 |   1.046 |      0.111 |    NA |            NA |    NA |    0.028 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + water_content + wr + fine_silt                                        | -51.088 |     1.239 |   2.237 |         NA | 0.068 |         0.046 | 0.041 |       NA |            NA |          NA |    NA |     0.033 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + ammonium                                      | -51.077 |     1.249 |   2.295 |      0.112 |    NA |         0.038 |    NA |    0.025 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ammonium + coarse_silt_sand                                   | -51.071 |     1.255 |   1.287 |      0.104 |    NA |            NA |    NA |    0.030 |            NA |          NA |    NA |        NA |            0.038 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + wr + finesilt_and_clay                        | -51.067 |     1.259 |   2.515 |      0.070 |    NA |         0.045 | 0.029 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.034 |       NA |
| JaccardDistance \~ p_h\_water + ec + clay + fine_silt                                         | -51.064 |     1.262 |   3.545 |      0.062 | 0.037 |            NA |    NA |       NA |            NA |          NA | 0.023 |     0.036 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + wr + oc_beregnet                         | -51.056 |     1.270 |   4.304 |      0.044 | 0.040 |         0.052 | 0.023 |       NA |            NA |       0.028 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + oc_beregnet + fine_silt                                       | -51.049 |     1.277 |   1.491 |      0.093 |    NA |            NA |    NA |       NA |            NA |       0.022 |    NA |     0.043 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + finesilt_and_clay + shannon_bak                                       | -51.045 |     1.281 |   2.414 |         NA | 0.055 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.051 |             0.040 |       NA |
| JaccardDistance \~ p_h\_water + oc_beregnet + fine_silt + shannon_bak                         | -51.037 |     1.289 |   7.510 |      0.066 |    NA |            NA |    NA |       NA |            NA |       0.029 |    NA |     0.036 |               NA |       0.038 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + wr + ammonium + shannon_bak                              | -51.026 |     1.301 |   9.256 |      0.058 | 0.046 |            NA | 0.031 |    0.026 |            NA |          NA |    NA |        NA |               NA |       0.046 |                NA |       NA |
| JaccardDistance \~ p_h\_water + wr + fine_silt + shannon_bak                                  | -51.021 |     1.305 |   8.891 |      0.070 |    NA |            NA | 0.028 |       NA |            NA |          NA |    NA |     0.035 |               NA |       0.041 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + nitrat_nitrit                            | -51.008 |     1.318 |   2.452 |      0.055 | 0.043 |         0.042 |    NA |       NA |         0.017 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + wr + coarse_silt_sand                         | -51.006 |     1.320 |   2.393 |      0.069 |    NA |         0.045 | 0.029 |       NA |            NA |          NA |    NA |        NA |            0.033 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + clay + coarse_silt_sand + shannon_bak                    | -50.993 |     1.333 |   7.807 |      0.039 | 0.041 |            NA |    NA |       NA |            NA |          NA | 0.029 |        NA |            0.029 |       0.036 |                NA |       NA |
| JaccardDistance \~ ec + water_content + wr + finesilt_and_clay                                | -50.993 |     1.333 |   1.912 |         NA | 0.069 |         0.049 | 0.041 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.032 |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + ammonium + oc_beregnet                   | -50.984 |     1.343 |   3.146 |      0.050 | 0.039 |         0.054 |    NA |    0.021 |            NA |       0.029 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + wr + clay                                | -50.979 |     1.347 |   4.223 |      0.042 | 0.041 |         0.051 | 0.027 |       NA |            NA |          NA | 0.027 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + clay + coarse_silt_sand                                       | -50.979 |     1.347 |   4.325 |      0.108 |    NA |            NA |    NA |       NA |            NA |          NA | 0.029 |        NA |            0.033 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + wr + coarse_silt_sand                    | -50.975 |     1.351 |   4.403 |      0.041 | 0.037 |         0.045 | 0.027 |       NA |            NA |          NA |    NA |        NA |            0.027 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + clay + fine_silt                         | -50.971 |     1.355 |   4.281 |      0.055 | 0.043 |         0.036 |    NA |       NA |            NA |          NA | 0.021 |     0.027 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + oc_beregnet + shannon_bak                     | -50.961 |     1.365 |   8.086 |      0.051 |    NA |         0.035 |    NA |       NA |            NA |       0.039 |    NA |        NA |               NA |       0.023 |                NA |       NA |
| JaccardDistance \~ ec + ammonium + fine_silt + shannon_bak                                    | -50.957 |     1.370 |   2.469 |         NA | 0.057 |            NA |    NA |    0.033 |            NA |          NA |    NA |     0.042 |               NA |       0.051 |                NA |       NA |
| JaccardDistance \~ ec + shannon_bak                                                           | -50.932 |     1.394 |   2.052 |         NA | 0.052 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.052 |                NA |       NA |
| JaccardDistance \~ ec + water_content + dexter_n + shannon_bak                                | -50.921 |     1.405 |   2.175 |         NA | 0.060 |         0.056 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.052 |                NA |    0.020 |
| JaccardDistance \~ p_h\_water + ec + wr + oc_beregnet + shannon_bak                           | -50.914 |     1.412 |   9.981 |      0.058 | 0.042 |            NA | 0.030 |       NA |            NA |       0.025 |    NA |        NA |               NA |       0.050 |                NA |       NA |
| JaccardDistance \~ p_h\_water + nitrat_nitrit                                                 | -50.899 |     1.427 |   1.157 |      0.102 |    NA |            NA |    NA |       NA |         0.025 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + ammonium + coarse_silt_sand + shannon_bak                | -50.895 |     1.431 |   7.332 |      0.038 | 0.042 |            NA |    NA |    0.028 |            NA |          NA |    NA |        NA |            0.030 |       0.035 |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + wr + oc_beregnet                              | -50.888 |     1.439 |   3.225 |      0.072 |    NA |         0.050 | 0.022 |       NA |            NA |       0.031 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + oc_beregnet + coarse_silt_sand + shannon_bak             | -50.885 |     1.441 |   9.130 |      0.042 | 0.038 |            NA |    NA |       NA |            NA |       0.028 |    NA |        NA |            0.030 |       0.039 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + oc_beregnet + finesilt_and_clay + shannon_bak            | -50.885 |     1.441 |   9.130 |      0.042 | 0.038 |            NA |    NA |       NA |            NA |       0.022 |    NA |        NA |               NA |       0.039 |             0.030 |       NA |
| JaccardDistance \~ p_h\_water + oc_beregnet + coarse_silt_sand + shannon_bak                  | -50.877 |     1.449 |   7.522 |      0.066 |    NA |            NA |    NA |       NA |            NA |       0.032 |    NA |        NA |            0.034 |       0.038 |                NA |       NA |
| JaccardDistance \~ p_h\_water + oc_beregnet + finesilt_and_clay + shannon_bak                 | -50.877 |     1.449 |   7.522 |      0.066 |    NA |            NA |    NA |       NA |            NA |       0.029 |    NA |        NA |               NA |       0.038 |             0.034 |       NA |
| JaccardDistance \~ ec + clay + coarse_silt_sand + shannon_bak                                 | -50.877 |     1.449 |   4.817 |         NA | 0.063 |            NA |    NA |       NA |            NA |          NA | 0.032 |        NA |            0.044 |       0.050 |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + ammonium + oc_beregnet                        | -50.869 |     1.458 |   2.874 |      0.109 |    NA |         0.053 |    NA |    0.022 |            NA |       0.035 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + oc_beregnet + coarse_silt_sand                                | -50.867 |     1.460 |   1.685 |      0.094 |    NA |            NA |    NA |       NA |            NA |       0.027 |    NA |        NA |            0.040 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + oc_beregnet + finesilt_and_clay                               | -50.867 |     1.460 |   1.453 |      0.094 |    NA |            NA |    NA |       NA |            NA |       0.023 |    NA |        NA |               NA |          NA |             0.040 |       NA |
| JaccardDistance \~ p_h\_water + wr + finesilt_and_clay + shannon_bak                          | -50.861 |     1.465 |   8.812 |      0.071 |    NA |            NA | 0.029 |       NA |            NA |          NA |    NA |        NA |               NA |       0.042 |             0.032 |       NA |
| JaccardDistance \~ p_h\_water + water_content + shannon_bak                                   | -50.858 |     1.468 |   7.726 |      0.042 |    NA |         0.024 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.022 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + clay + shannon_bak                       | -50.858 |     1.469 |   8.234 |      0.029 | 0.047 |         0.027 |    NA |       NA |            NA |          NA | 0.027 |        NA |               NA |       0.025 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + ammonium + finesilt_and_clay             | -50.857 |     1.469 |   2.941 |      0.057 | 0.043 |         0.032 |    NA |    0.021 |            NA |          NA |    NA |        NA |               NA |          NA |             0.027 |       NA |
| JaccardDistance \~ ec + water_content + fine_silt                                             | -50.855 |     1.471 |   1.908 |         NA | 0.076 |         0.044 |    NA |       NA |            NA |          NA |    NA |     0.033 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + ammonium + coarse_silt_sand + shannon_bak                             | -50.855 |     1.471 |   2.619 |         NA | 0.059 |            NA |    NA |    0.032 |            NA |          NA |    NA |        NA |            0.040 |       0.052 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + wr + clay + shannon_bak                                  | -50.854 |     1.472 |   9.528 |      0.058 | 0.043 |            NA | 0.027 |       NA |            NA |          NA | 0.024 |        NA |               NA |       0.049 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + coarse_silt_sand + shannon_bak           | -50.842 |     1.484 |   8.345 |      0.029 | 0.042 |         0.027 |    NA |       NA |            NA |          NA |    NA |        NA |            0.027 |       0.025 |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + fine_silt + coarse_silt_sand                  | -50.840 |     1.486 |   8.749 |      0.088 |    NA |         0.041 |    NA |       NA |            NA |          NA |    NA |     0.027 |            0.024 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ammonium + fine_silt + shannon_bak                            | -50.832 |     1.494 |   6.658 |      0.056 |    NA |            NA |    NA |    0.026 |            NA |          NA |    NA |     0.038 |               NA |       0.028 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + oc_beregnet + coarse_silt_sand                           | -50.830 |     1.497 |   3.215 |      0.048 | 0.037 |            NA |    NA |       NA |            NA |       0.027 |    NA |        NA |            0.039 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + oc_beregnet + finesilt_and_clay                          | -50.830 |     1.497 |   3.215 |      0.048 | 0.037 |            NA |    NA |       NA |            NA |       0.024 |    NA |        NA |               NA |          NA |             0.039 |       NA |
| JaccardDistance \~ p_h\_water + ec + wr + fine_silt                                           | -50.827 |     1.499 |   4.327 |      0.042 | 0.038 |            NA | 0.020 |       NA |            NA |          NA |    NA |     0.039 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + oc_beregnet + fine_silt + shannon_bak                                 | -50.825 |     1.501 |   2.695 |         NA | 0.063 |            NA |    NA |       NA |            NA |       0.031 |    NA |     0.043 |               NA |       0.044 |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + oc_beregnet + fine_silt + shannon_bak         | -50.824 |     1.502 |   8.118 |      0.050 |    NA |         0.035 |    NA |       NA |            NA |       0.040 |    NA |     0.036 |               NA |       0.024 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + fine_silt + coarse_silt_sand             | -50.824 |     1.502 |   9.977 |      0.049 | 0.038 |         0.041 |    NA |       NA |            NA |          NA |    NA |     0.025 |            0.019 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ammonium + shannon_bak                                        | -50.823 |     1.503 |   6.303 |      0.057 |    NA |            NA |    NA |    0.024 |            NA |          NA |    NA |        NA |               NA |       0.034 |                NA |       NA |
| JaccardDistance \~ ec + water_content + wr + fine_silt + shannon_bak                          | -50.816 |     1.510 |   2.968 |         NA | 0.058 |         0.044 | 0.024 |       NA |            NA |          NA |    NA |     0.033 |               NA |       0.034 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + ammonium + clay                          | -50.810 |     1.517 |   2.618 |      0.056 | 0.044 |         0.039 |    NA |    0.024 |            NA |          NA | 0.027 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + ammonium + fine_silt                     | -50.808 |     1.519 |   3.233 |      0.056 | 0.042 |         0.028 |    NA |    0.019 |            NA |          NA |    NA |     0.027 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ammonium + clay                                               | -50.806 |     1.521 |   1.166 |      0.104 |    NA |            NA |    NA |    0.030 |            NA |          NA | 0.033 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + wr + fine_silt                                                | -50.797 |     1.529 |   2.324 |      0.070 |    NA |            NA | 0.018 |       NA |            NA |          NA |    NA |     0.041 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + dexter_n                                      | -50.790 |     1.536 |   1.041 |      0.112 |    NA |         0.042 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.021 |
| JaccardDistance \~ p_h\_water + clay + coarse_silt_sand + shannon_bak                         | -50.789 |     1.537 |   6.973 |      0.062 |    NA |            NA |    NA |       NA |            NA |          NA | 0.030 |        NA |            0.033 |       0.035 |                NA |       NA |
| JaccardDistance \~ ec + water_content + finesilt_and_clay                                     | -50.777 |     1.550 |   1.726 |         NA | 0.080 |         0.047 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.032 |       NA |
| JaccardDistance \~ p_h\_water + wr + coarse_silt_sand + shannon_bak                           | -50.771 |     1.555 |   9.033 |      0.070 |    NA |            NA | 0.030 |       NA |            NA |          NA |    NA |        NA |            0.031 |       0.041 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + ammonium + clay + shannon_bak                            | -50.768 |     1.558 |   6.666 |      0.045 | 0.045 |            NA |    NA |    0.025 |            NA |          NA | 0.028 |        NA |               NA |       0.038 |                NA |       NA |
| JaccardDistance \~ ec + ammonium + finesilt_and_clay + shannon_bak                            | -50.765 |     1.561 |   2.422 |         NA | 0.057 |            NA |    NA |    0.034 |            NA |          NA |    NA |        NA |               NA |       0.051 |             0.039 |       NA |
| JaccardDistance \~ p_h\_water + ec + fine_silt + dexter_n + shannon_bak                       | -50.764 |     1.563 |   6.386 |      0.049 | 0.045 |            NA |    NA |       NA |            NA |          NA |    NA |     0.032 |               NA |       0.040 |                NA |    0.017 |
| JaccardDistance \~ p_h\_water + ec + oc_beregnet                                              | -50.763 |     1.563 |   3.120 |      0.049 | 0.039 |            NA |    NA |       NA |            NA |       0.024 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + water_content + wr + coarse_silt_sand                                 | -50.757 |     1.570 |   2.367 |         NA | 0.065 |         0.045 | 0.041 |       NA |            NA |          NA |    NA |        NA |            0.028 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + wr + clay                                     | -50.748 |     1.578 |   2.420 |      0.071 |    NA |         0.048 | 0.029 |       NA |            NA |          NA | 0.029 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + wr                                                       | -50.747 |     1.579 |   4.030 |      0.048 | 0.041 |            NA | 0.024 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + ammonium + clay                                          | -50.745 |     1.581 |   2.292 |      0.059 | 0.037 |            NA |    NA |    0.031 |            NA |          NA | 0.032 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + wr + oc_beregnet + fine_silt                  | -50.744 |     1.582 |   3.349 |      0.068 |    NA |         0.047 | 0.023 |       NA |            NA |       0.031 |    NA |     0.036 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + wr + ammonium                            | -50.742 |     1.584 |   4.187 |      0.042 | 0.043 |         0.042 | 0.026 |    0.023 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + oc_beregnet + clay + fine_silt                | -50.740 |     1.587 |   4.362 |      0.092 |    NA |         0.048 |    NA |       NA |            NA |       0.039 | 0.023 |     0.030 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + nitrat_nitrit + shannon_bak                                   | -50.736 |     1.590 |   5.752 |      0.064 |    NA |            NA |    NA |       NA |         0.022 |          NA |    NA |        NA |               NA |       0.035 |                NA |       NA |
| JaccardDistance \~ p_h\_water + oc_beregnet                                                   | -50.734 |     1.592 |   1.211 |      0.103 |    NA |            NA |    NA |       NA |            NA |       0.022 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + fine_silt + dexter_n                                          | -50.729 |     1.597 |   1.219 |      0.105 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.041 |               NA |          NA |                NA |    0.017 |
| JaccardDistance \~ ec + water_content + oc_beregnet + fine_silt + shannon_bak                 | -50.729 |     1.597 |   3.225 |         NA | 0.049 |         0.036 |    NA |       NA |            NA |       0.023 |    NA |     0.029 |               NA |       0.039 |                NA |       NA |
| JaccardDistance \~ ec + water_content + nitrat_nitrit + shannon_bak                           | -50.728 |     1.599 |   2.343 |         NA | 0.056 |         0.050 |    NA |       NA |         0.017 |          NA |    NA |        NA |               NA |       0.051 |                NA |       NA |
| JaccardDistance \~ ec + water_content + wr + finesilt_and_clay + shannon_bak                  | -50.724 |     1.602 |   2.953 |         NA | 0.058 |         0.047 | 0.024 |       NA |            NA |          NA |    NA |        NA |               NA |       0.034 |             0.032 |       NA |
| JaccardDistance \~ ec + water_content + wr + clay                                             | -50.721 |     1.605 |   1.469 |         NA | 0.071 |         0.054 | 0.041 |       NA |            NA |          NA | 0.028 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + nitrat_nitrit                                            | -50.714 |     1.612 |   2.212 |      0.063 | 0.035 |            NA |    NA |       NA |         0.023 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + ammonium + shannon_bak                                                | -50.695 |     1.631 |   2.092 |         NA | 0.055 |            NA |    NA |    0.034 |            NA |          NA |    NA |        NA |               NA |       0.052 |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + nitrat_nitrit                                 | -50.691 |     1.635 |   1.424 |      0.090 |    NA |         0.035 |    NA |       NA |         0.019 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + nitrat_nitrit + oc_beregnet                   | -50.689 |     1.637 |   2.339 |      0.091 |    NA |         0.052 |    NA |       NA |         0.019 |       0.038 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ammonium + finesilt_and_clay + shannon_bak                    | -50.689 |     1.638 |   6.617 |      0.056 |    NA |            NA |    NA |    0.026 |            NA |          NA |    NA |        NA |               NA |       0.029 |             0.036 |       NA |
| JaccardDistance \~ ec + water_content + oc_beregnet + fine_silt                               | -50.685 |     1.641 |   2.921 |         NA | 0.077 |         0.042 |    NA |       NA |            NA |       0.035 |    NA |     0.030 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + clay + fine_silt + shannon_bak                                | -50.682 |     1.645 |   6.161 |      0.062 |    NA |            NA |    NA |       NA |            NA |          NA | 0.023 |     0.031 |               NA |       0.031 |                NA |       NA |
| JaccardDistance \~ p_h\_water + fine_silt + coarse_silt_sand                                  | -50.667 |     1.659 |   7.633 |      0.094 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.024 |            0.016 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + fine_silt + coarse_silt_sand + shannon_bak               | -50.663 |     1.663 |   8.656 |      0.044 | 0.040 |            NA |    NA |       NA |            NA |          NA |    NA |     0.025 |            0.016 |       0.039 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + ammonium + shannon_bak                   | -50.659 |     1.667 |   8.214 |      0.030 | 0.049 |         0.026 |    NA |    0.024 |            NA |          NA |    NA |        NA |               NA |       0.025 |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + oc_beregnet + coarse_silt_sand + shannon_bak  | -50.659 |     1.668 |   8.092 |      0.051 |    NA |         0.035 |    NA |       NA |            NA |       0.039 |    NA |        NA |            0.033 |       0.024 |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + oc_beregnet + finesilt_and_clay + shannon_bak | -50.659 |     1.668 |   8.092 |      0.051 |    NA |         0.035 |    NA |       NA |            NA |       0.039 |    NA |        NA |               NA |       0.024 |             0.033 |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + ammonium + coarse_silt_sand              | -50.657 |     1.669 |   3.356 |      0.055 | 0.039 |         0.031 |    NA |    0.023 |            NA |          NA |    NA |        NA |            0.025 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + fine_silt + shannon_bak                       | -50.652 |     1.674 |   7.779 |      0.041 |    NA |         0.023 |    NA |       NA |            NA |          NA |    NA |     0.035 |               NA |       0.022 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + nitrat_nitrit + oc_beregnet              | -50.651 |     1.675 |   3.146 |      0.050 | 0.037 |         0.051 |    NA |       NA |         0.017 |       0.033 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ water_content + shannon_bak                                                | -50.647 |     1.679 |   1.023 |         NA |    NA |         0.047 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.092 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + fine_silt + dexter_n                                     | -50.644 |     1.683 |   2.668 |      0.061 | 0.037 |            NA |    NA |       NA |            NA |          NA |    NA |     0.040 |               NA |          NA |                NA |    0.017 |
| JaccardDistance \~ p_h\_water + water_content + clay + fine_silt                              | -50.642 |     1.684 |   4.079 |      0.088 |    NA |         0.031 |    NA |       NA |            NA |          NA | 0.021 |     0.027 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + water_content + clay + fine_silt + shannon_bak                        | -50.642 |     1.684 |   4.230 |         NA | 0.057 |         0.043 |    NA |       NA |            NA |          NA | 0.022 |     0.026 |               NA |       0.051 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ammonium + coarse_silt_sand + shannon_bak                     | -50.637 |     1.690 |   6.933 |      0.056 |    NA |            NA |    NA |    0.028 |            NA |          NA |    NA |        NA |            0.035 |       0.031 |                NA |       NA |
| JaccardDistance \~ ec + water_content + oc_beregnet + coarse_silt_sand + shannon_bak          | -50.627 |     1.699 |   2.860 |         NA | 0.050 |         0.038 |    NA |       NA |            NA |       0.026 |    NA |        NA |            0.027 |       0.039 |                NA |       NA |
| JaccardDistance \~ ec + water_content + oc_beregnet + finesilt_and_clay + shannon_bak         | -50.627 |     1.699 |   2.860 |         NA | 0.050 |         0.038 |    NA |       NA |            NA |       0.023 |    NA |        NA |               NA |       0.039 |             0.027 |       NA |
| JaccardDistance \~ p_h\_water + ec + nitrat_nitrit + fine_silt + shannon_bak                  | -50.619 |     1.707 |   6.403 |      0.050 | 0.045 |            NA |    NA |       NA |         0.015 |          NA |    NA |     0.030 |               NA |       0.039 |                NA |       NA |
| JaccardDistance \~ ec + water_content + wr + oc_beregnet                                      | -50.618 |     1.708 |   3.326 |         NA | 0.068 |         0.052 | 0.029 |       NA |            NA |       0.026 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + wr + oc_beregnet + shannon_bak                                | -50.608 |     1.718 |   8.811 |      0.073 |    NA |            NA | 0.030 |       NA |            NA |       0.028 |    NA |        NA |               NA |       0.046 |                NA |       NA |
| JaccardDistance \~ p_h\_water + nitrat_nitrit + fine_silt                                     | -50.599 |     1.727 |   1.645 |      0.101 |    NA |            NA |    NA |       NA |         0.015 |          NA |    NA |     0.033 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + oc_beregnet + coarse_silt_sand + shannon_bak                          | -50.597 |     1.729 |   2.649 |         NA | 0.062 |            NA |    NA |       NA |            NA |       0.028 |    NA |        NA |            0.039 |       0.044 |                NA |       NA |
| JaccardDistance \~ ec + oc_beregnet + finesilt_and_clay + shannon_bak                         | -50.597 |     1.729 |   2.649 |         NA | 0.062 |            NA |    NA |       NA |            NA |       0.031 |    NA |        NA |               NA |       0.044 |             0.039 |       NA |
| JaccardDistance \~ ec + coarse_silt_sand + finesilt_and_clay + shannon_bak                    | -50.597 |     1.729 |   9.551 |         NA | 0.062 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.031 |       0.044 |             0.028 |       NA |
| JaccardDistance \~ p_h\_water + wr                                                            | -50.596 |     1.731 |   2.230 |      0.068 |    NA |            NA | 0.020 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + wr + finesilt_and_clay                                        | -50.595 |     1.731 |   2.323 |      0.070 |    NA |            NA | 0.019 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.038 |       NA |
| JaccardDistance \~ p_h\_water + ec + fine_silt + coarse_silt_sand                             | -50.594 |     1.732 |   8.073 |      0.050 | 0.037 |            NA |    NA |       NA |            NA |          NA |    NA |     0.023 |            0.016 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + wr + clay + coarse_silt_sand                  | -50.589 |     1.737 |   9.942 |      0.068 |    NA |         0.046 | 0.027 |       NA |            NA |          NA | 0.032 |        NA |            0.036 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + wr + finesilt_and_clay                                   | -50.589 |     1.737 |   4.266 |      0.043 | 0.038 |            NA | 0.020 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.036 |       NA |
| JaccardDistance \~ ec + water_content + oc_beregnet + coarse_silt_sand                        | -50.583 |     1.743 |   2.593 |         NA | 0.081 |         0.044 |    NA |       NA |            NA |       0.039 |    NA |        NA |            0.028 |          NA |                NA |       NA |
| JaccardDistance \~ ec + water_content + oc_beregnet + finesilt_and_clay                       | -50.583 |     1.743 |   2.593 |         NA | 0.081 |         0.044 |    NA |       NA |            NA |       0.035 |    NA |        NA |               NA |          NA |             0.028 |       NA |
| JaccardDistance \~ p_h\_water + dexter_n                                                      | -50.581 |     1.745 |   1.023 |      0.114 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.019 |
| JaccardDistance \~ p_h\_water + water_content + ammonium + coarse_silt_sand                   | -50.578 |     1.748 |   3.352 |      0.101 |    NA |         0.030 |    NA |    0.023 |            NA |          NA |    NA |        NA |            0.030 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + wr + ammonium + fine_silt                                | -50.577 |     1.749 |   4.327 |      0.042 | 0.038 |            NA | 0.024 |    0.034 |            NA |          NA |    NA |     0.039 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + wr + fine_silt                                                        | -50.577 |     1.750 |   1.615 |         NA | 0.066 |            NA | 0.040 |       NA |            NA |          NA |    NA |     0.045 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + wr + oc_beregnet + coarse_silt_sand           | -50.575 |     1.751 |   3.302 |      0.069 |    NA |         0.046 | 0.023 |       NA |            NA |       0.031 |    NA |        NA |            0.033 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + wr + oc_beregnet + finesilt_and_clay          | -50.575 |     1.751 |   3.302 |      0.069 |    NA |         0.046 | 0.023 |       NA |            NA |       0.031 |    NA |        NA |               NA |          NA |             0.033 |       NA |
| JaccardDistance \~ p_h\_water + water_content + coarse_silt_sand + shannon_bak                | -50.567 |     1.759 |   7.852 |      0.042 |    NA |         0.027 |    NA |       NA |            NA |          NA |    NA |        NA |            0.034 |       0.023 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + nitrat_nitrit + finesilt_and_clay + shannon_bak          | -50.565 |     1.762 |   6.407 |      0.051 | 0.044 |            NA |    NA |       NA |         0.017 |          NA |    NA |        NA |               NA |       0.039 |             0.029 |       NA |
| JaccardDistance \~ p_h\_water + water_content + finesilt_and_clay + shannon_bak               | -50.558 |     1.768 |   7.817 |      0.041 |    NA |         0.024 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.021 |             0.033 |       NA |
| JaccardDistance \~ p_h\_water + wr + coarse_silt_sand                                         | -50.556 |     1.770 |   2.326 |      0.070 |    NA |            NA | 0.022 |       NA |            NA |          NA |    NA |        NA |            0.038 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + nitrat_nitrit + fine_silt                                | -50.555 |     1.771 |   2.560 |      0.062 | 0.037 |            NA |    NA |       NA |         0.016 |          NA |    NA |     0.036 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + wr + ammonium + fine_silt                                     | -50.546 |     1.780 |   2.573 |      0.070 |    NA |            NA | 0.024 |    0.034 |            NA |          NA |    NA |     0.041 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + water_content + coarse_silt_sand                                      | -50.539 |     1.787 |   2.358 |         NA | 0.083 |         0.040 |    NA |       NA |            NA |          NA |    NA |        NA |            0.028 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + wr + clay + shannon_bak                                       | -50.539 |     1.787 |   8.592 |      0.072 |    NA |            NA | 0.029 |       NA |            NA |          NA | 0.027 |        NA |               NA |       0.044 |                NA |       NA |
| JaccardDistance \~ p_h\_water + nitrat_nitrit + finesilt_and_clay                             | -50.538 |     1.788 |   1.461 |      0.101 |    NA |            NA |    NA |       NA |         0.018 |          NA |    NA |        NA |               NA |          NA |             0.032 |       NA |
| JaccardDistance \~ p_h\_water + ec + wr + dexter_n + shannon_bak                              | -50.536 |     1.790 |   8.933 |      0.065 | 0.046 |            NA | 0.028 |       NA |            NA |          NA |    NA |        NA |               NA |       0.054 |                NA |    0.019 |
| JaccardDistance \~ p_h\_water + fine_silt + coarse_silt_sand + shannon_bak                    | -50.535 |     1.791 |   7.768 |      0.065 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.026 |            0.021 |       0.036 |                NA |       NA |
| JaccardDistance \~ p_h\_water + dexter_n + shannon_bak                                        | -50.535 |     1.792 |   5.332 |      0.064 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.038 |                NA |    0.019 |
| JaccardDistance \~ p_h\_water + water_content + ammonium + fine_silt                          | -50.532 |     1.794 |   3.188 |      0.091 |    NA |         0.024 |    NA |    0.020 |            NA |          NA |    NA |     0.030 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + water_content + clay                                                  | -50.530 |     1.796 |   1.357 |         NA | 0.089 |         0.052 |    NA |       NA |            NA |          NA | 0.028 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + ammonium + finesilt_and_clay                  | -50.530 |     1.796 |   2.857 |      0.095 |    NA |         0.026 |    NA |    0.021 |            NA |          NA |    NA |        NA |               NA |          NA |             0.030 |       NA |
| JaccardDistance \~ p_h\_water + ec + ammonium + clay + fine_silt                              | -50.529 |     1.798 |   3.578 |      0.059 | 0.038 |            NA |    NA |    0.030 |            NA |          NA | 0.023 |     0.035 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + oc_beregnet + shannon_bak                                             | -50.528 |     1.798 |   2.406 |         NA | 0.059 |            NA |    NA |       NA |            NA |       0.032 |    NA |        NA |               NA |       0.045 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + oc_beregnet + clay + shannon_bak                         | -50.527 |     1.800 |   8.948 |      0.045 | 0.038 |            NA |    NA |       NA |            NA |       0.022 | 0.024 |        NA |               NA |       0.042 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + dexter_n                                                 | -50.524 |     1.802 |   1.910 |      0.063 | 0.037 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.020 |
| JaccardDistance \~ p_h\_water + ec + wr + coarse_silt_sand                                    | -50.520 |     1.806 |   4.403 |      0.041 | 0.038 |            NA | 0.022 |       NA |            NA |          NA |    NA |        NA |            0.035 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + oc_beregnet + clay + shannon_bak                              | -50.514 |     1.812 |   7.519 |      0.067 |    NA |            NA |    NA |       NA |            NA |       0.029 | 0.029 |        NA |               NA |       0.040 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + nitrat_nitrit + finesilt_and_clay        | -50.513 |     1.813 |   2.660 |      0.055 | 0.042 |         0.039 |    NA |       NA |         0.016 |          NA |    NA |        NA |               NA |          NA |             0.031 |       NA |
| JaccardDistance \~ ec + clay + shannon_bak                                                    | -50.512 |     1.814 |   2.263 |         NA | 0.053 |            NA |    NA |       NA |            NA |          NA | 0.031 |        NA |               NA |       0.052 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ammonium + clay + fine_silt                                   | -50.509 |     1.818 |   3.215 |      0.098 |    NA |            NA |    NA |    0.029 |            NA |          NA | 0.023 |     0.033 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + nitrat_nitrit + fine_silt                | -50.498 |     1.828 |   2.647 |      0.055 | 0.042 |         0.037 |    NA |       NA |         0.015 |          NA |    NA |     0.030 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + water_content + wr + ammonium                                         | -50.488 |     1.838 |   2.437 |         NA | 0.074 |         0.047 | 0.040 |    0.024 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + wr + ammonium + shannon_bak                                   | -50.487 |     1.840 |   8.636 |      0.069 |    NA |            NA | 0.033 |    0.026 |            NA |          NA |    NA |        NA |               NA |       0.040 |                NA |       NA |
| JaccardDistance \~ ec + water_content + ammonium + finesilt_and_clay + shannon_bak            | -50.481 |     1.845 |   2.818 |         NA | 0.057 |         0.034 |    NA |    0.020 |            NA |          NA |    NA |        NA |               NA |       0.051 |             0.027 |       NA |
| JaccardDistance \~ p_h\_water + ec + wr + ammonium                                            | -50.481 |     1.846 |   4.050 |      0.047 | 0.040 |            NA | 0.028 |    0.034 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + fine_silt                                                             | -50.476 |     1.851 |   1.554 |         NA | 0.080 |            NA |    NA |       NA |            NA |          NA |    NA |     0.044 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + nitrat_nitrit + finesilt_and_clay                        | -50.468 |     1.859 |   2.556 |      0.062 | 0.037 |            NA |    NA |       NA |         0.018 |          NA |    NA |        NA |               NA |          NA |             0.034 |       NA |
| JaccardDistance \~ water_content + fine_silt + shannon_bak                                    | -50.467 |     1.859 |   1.691 |         NA |    NA |         0.043 |    NA |       NA |            NA |          NA |    NA |     0.035 |               NA |       0.070 |                NA |       NA |
| JaccardDistance \~ p_h\_water + wr + ammonium + fine_silt + shannon_bak                       | -50.465 |     1.861 |   9.065 |      0.066 |    NA |            NA | 0.032 |    0.030 |            NA |          NA |    NA |     0.038 |               NA |       0.037 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + fine_silt + dexter_n                     | -50.462 |     1.864 |   2.717 |      0.056 | 0.041 |         0.035 |    NA |       NA |            NA |          NA |    NA |     0.027 |               NA |          NA |                NA |    0.014 |
| JaccardDistance \~ ec + water_content + oc_beregnet + clay + shannon_bak                      | -50.459 |     1.867 |   2.673 |         NA | 0.053 |         0.044 |    NA |       NA |            NA |       0.024 | 0.025 |        NA |               NA |       0.039 |                NA |       NA |
| JaccardDistance \~ ec + water_content + ammonium + clay + shannon_bak                         | -50.457 |     1.869 |   2.530 |         NA | 0.059 |         0.041 |    NA |    0.024 |            NA |          NA | 0.027 |        NA |               NA |       0.051 |                NA |       NA |
| JaccardDistance \~ ec + water_content + wr + coarse_silt_sand + shannon_bak                   | -50.451 |     1.875 |   2.985 |         NA | 0.055 |         0.042 | 0.023 |       NA |            NA |          NA |    NA |        NA |            0.028 |       0.033 |                NA |       NA |
| JaccardDistance \~ ec + water_content + wr + clay + shannon_bak                               | -50.451 |     1.875 |   2.942 |         NA | 0.060 |         0.052 | 0.023 |       NA |            NA |          NA | 0.028 |        NA |               NA |       0.034 |                NA |       NA |
| JaccardDistance \~ ec + oc_beregnet + fine_silt                                               | -50.442 |     1.884 |   1.560 |         NA | 0.083 |            NA |    NA |       NA |            NA |       0.038 |    NA |     0.044 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + wr + nitrat_nitrit + shannon_bak                         | -50.432 |     1.894 |   9.320 |      0.060 | 0.042 |            NA | 0.027 |       NA |         0.018 |          NA |    NA |        NA |               NA |       0.050 |                NA |       NA |
| JaccardDistance \~ ec + water_content + fine_silt + coarse_silt_sand + shannon_bak            | -50.432 |     1.894 |   9.070 |         NA | 0.054 |         0.041 |    NA |       NA |            NA |          NA |    NA |     0.023 |            0.019 |       0.043 |                NA |       NA |
| JaccardDistance \~ ec + wr + fine_silt + shannon_bak                                          | -50.429 |     1.897 |   2.839 |         NA | 0.061 |            NA | 0.025 |       NA |            NA |          NA |    NA |     0.045 |               NA |       0.036 |                NA |       NA |
| JaccardDistance \~ ec + water_content + ammonium + fine_silt + shannon_bak                    | -50.426 |     1.900 |   3.077 |         NA | 0.056 |         0.030 |    NA |    0.019 |            NA |          NA |    NA |     0.026 |               NA |       0.051 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + wr + dexter_n                            | -50.424 |     1.902 |   4.166 |      0.042 | 0.043 |         0.052 | 0.026 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.019 |
| JaccardDistance \~ ec + coarse_silt_sand                                                      | -50.417 |     1.909 |   1.437 |         NA | 0.093 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.043 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + ammonium + clay                               | -50.411 |     1.915 |   2.434 |      0.103 |    NA |         0.032 |    NA |    0.024 |            NA |          NA | 0.028 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + water_content + ammonium + oc_beregnet + fine_silt            | -50.407 |     1.919 |   4.299 |      0.093 |    NA |         0.041 |    NA |    0.018 |            NA |       0.036 |    NA |     0.031 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ ec + water_content + oc_beregnet + clay                                    | -50.403 |     1.923 |   2.166 |         NA | 0.089 |         0.050 |    NA |       NA |            NA |       0.036 | 0.025 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + wr + ammonium + finesilt_and_clay                             | -50.398 |     1.928 |   2.546 |      0.070 |    NA |            NA | 0.024 |    0.035 |            NA |          NA |    NA |        NA |               NA |          NA |             0.039 |       NA |
| JaccardDistance \~ p_h\_water + oc_beregnet + clay                                            | -50.396 |     1.930 |   1.353 |      0.098 |    NA |            NA |    NA |       NA |            NA |       0.023 | 0.033 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + wr + ammonium + finesilt_and_clay                        | -50.387 |     1.939 |   4.267 |      0.043 | 0.038 |            NA | 0.025 |    0.035 |            NA |          NA |    NA |        NA |               NA |          NA |             0.037 |       NA |
| JaccardDistance \~ p_h\_water + water_content + wr + ammonium                                 | -50.372 |     1.954 |   2.480 |      0.072 |    NA |         0.039 | 0.027 |    0.023 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + wr + oc_beregnet + fine_silt + shannon_bak                    | -50.368 |     1.958 |   9.199 |      0.070 |    NA |            NA | 0.028 |       NA |            NA |       0.028 |    NA |     0.034 |               NA |       0.041 |                NA |       NA |
| JaccardDistance \~ ec + water_content + ammonium                                              | -50.358 |     1.968 |   2.417 |         NA | 0.100 |         0.046 |    NA |    0.025 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ water_content + finesilt_and_clay + shannon_bak                            | -50.357 |     1.969 |   1.550 |         NA |    NA |         0.044 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.073 |             0.034 |       NA |
| JaccardDistance \~ water_content + coarse_silt_sand + shannon_bak                             | -50.353 |     1.974 |   1.916 |         NA |    NA |         0.046 |    NA |       NA |            NA |          NA |    NA |        NA |            0.034 |       0.080 |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + nitrat_nitrit + clay + shannon_bak                       | -50.351 |     1.975 |   6.345 |      0.053 | 0.044 |            NA |    NA |       NA |         0.019 |          NA | 0.026 |        NA |               NA |       0.041 |                NA |       NA |
| JaccardDistance \~ ec + water_content + ammonium + coarse_silt_sand + shannon_bak             | -50.348 |     1.978 |   3.368 |         NA | 0.054 |         0.030 |    NA |    0.022 |            NA |          NA |    NA |        NA |            0.025 |       0.051 |                NA |       NA |
| JaccardDistance \~ ec + wr + coarse_silt_sand                                                 | -50.346 |     1.981 |   1.682 |         NA | 0.066 |            NA | 0.037 |       NA |            NA |          NA |    NA |        NA |            0.042 |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + wr + ammonium                                                 | -50.344 |     1.983 |   2.477 |      0.070 |    NA |            NA | 0.026 |    0.034 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ p_h\_water + ec + water_content + dexter_n + shannon_bak                   | -50.342 |     1.984 |   8.263 |      0.029 | 0.049 |         0.026 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.025 |                NA |    0.019 |
| JaccardDistance \~ p_h\_water + wr + ammonium + finesilt_and_clay + shannon_bak               | -50.328 |     1.998 |   9.005 |      0.067 |    NA |            NA | 0.033 |    0.030 |            NA |          NA |    NA |        NA |               NA |       0.037 |             0.036 |       NA |

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
  Vars <- Vars[Vars != "area"]
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
#> [1] "1 of 3 ready 2023-03-23 11:41:07"
#> [1] "2 of 3 ready 2023-03-23 11:42:10"
#> [1] "3 of 3 ready 2023-03-23 11:43:13"


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

This generate up to 1,536 models to evaluate, which can be downloaded as
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
  Vars <- Vars[Vars != "area"]
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
  Temp$AICc <-  try(AICc.PERMANOVA2(with(Response,adonis2(as.formula(AllFormsVegAbund$Form[x]), data = Response, by = "margin", strata = area)))$AICc, silent = T)
  
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
  
  Rs <- broom::tidy(with(Response,adonis2(as.formula(AllFormsVegAbund$Form[x]), data = Response, by = "margin"))) %>% dplyr::filter(!(term %in% c("Residual", "Total"))) %>% dplyr::select(term, R2) %>%  pivot_wider(names_from = term, values_from = R2)
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

As seen in table <a href="#tab:SummaryVegAbund">2.3</a> there are 68
models within 2 AICc of each other, you can see there how many times a
variable has been selected

| Variable          | Number_of_models | Full_Akaike_Adjusted_RSq | Subset_Akaike_Adjusted_RSq |
|:------------------|-----------------:|-------------------------:|---------------------------:|
| p_h\_water        |               68 |                    0.089 |                      0.089 |
| shannon_bak       |               62 |                    0.056 |                      0.060 |
| ec                |               36 |                    0.021 |                      0.038 |
| water_content     |               25 |                    0.015 |                      0.039 |
| ammonium          |               24 |                    0.013 |                      0.035 |
| fine_silt         |               13 |                    0.005 |                      0.026 |
| oc_beregnet       |               12 |                    0.004 |                      0.025 |
| nitrat_nitrit     |                9 |                    0.002 |                      0.021 |
| coarse_silt_sand  |                9 |                    0.002 |                      0.021 |
| wr                |                8 |                    0.003 |                      0.026 |
| finesilt_and_clay |                7 |                    0.002 |                      0.022 |
| clay              |                5 |                    0.001 |                      0.015 |
| dexter_n          |                3 |                    0.000 |                      0.011 |

Table 2.3: Number of selected models were variables are present and
their Akaike Weighted R squared for the Marginal effect of the terms

Now we can see the top models that have a delta AICc within 2 of the
best model in table <a href="#tab:BestVegAbundModels">2.4</a> if
expanded

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Show table of selected models for vegetation presence absence

</summary>

| Form                                                                                |    AICc | DeltaAICc | Max_VIF | p_h\_water |    ec | water_content |    wr | ammonium | nitrat_nitrit | oc_beregnet |  clay | fine_silt | coarse_silt_sand | shannon_bak | finesilt_and_clay | dexter_n |
|:------------------------------------------------------------------------------------|--------:|----------:|--------:|-----------:|------:|--------------:|------:|---------:|--------------:|------------:|------:|----------:|-----------------:|------------:|------------------:|---------:|
| BrayDistance \~ p_h\_water + ec + ammonium + shannon_bak                            | -60.587 |     0.000 |   6.496 |      0.099 | 0.042 |            NA |    NA |    0.038 |            NA |          NA |    NA |        NA |               NA |       0.071 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + shannon_bak                       | -60.435 |     0.151 |   8.200 |      0.066 | 0.037 |         0.036 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.053 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + shannon_bak                                       | -60.196 |     0.391 |   5.976 |      0.115 | 0.045 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.076 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + shannon_bak                            | -60.183 |     0.404 |   7.726 |      0.087 |    NA |         0.045 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.047 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + wr + shannon_bak                                  | -60.044 |     0.543 |   8.770 |      0.084 | 0.051 |            NA | 0.031 |       NA |            NA |          NA |    NA |        NA |               NA |       0.074 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + oc_beregnet + shannon_bak         | -60.010 |     0.577 |   9.240 |      0.055 | 0.044 |         0.035 |    NA |       NA |            NA |       0.028 |    NA |        NA |               NA |       0.052 |                NA |       NA |
| BrayDistance \~ p_h\_water + ammonium + shannon_bak                                 | -59.965 |     0.621 |   6.303 |      0.111 |    NA |            NA |    NA |    0.042 |            NA |          NA |    NA |        NA |               NA |       0.064 |                NA |       NA |
| BrayDistance \~ p_h\_water + ammonium + fine_silt + shannon_bak                     | -59.956 |     0.631 |   6.658 |      0.110 |    NA |            NA |    NA |    0.042 |            NA |          NA |    NA |     0.033 |               NA |       0.061 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + wr + ammonium + shannon_bak                       | -59.936 |     0.651 |   9.256 |      0.085 | 0.048 |            NA | 0.025 |    0.032 |            NA |          NA |    NA |        NA |               NA |       0.071 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + ammonium + shannon_bak            | -59.872 |     0.715 |   8.214 |      0.065 | 0.036 |         0.024 |    NA |    0.026 |            NA |          NA |    NA |        NA |               NA |       0.050 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + oc_beregnet + shannon_bak                         | -59.862 |     0.725 |   8.520 |      0.061 | 0.045 |            NA |    NA |       NA |            NA |       0.029 |    NA |        NA |               NA |       0.061 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + ammonium + fine_silt + shannon_bak                | -59.786 |     0.800 |   6.796 |      0.092 | 0.031 |            NA |    NA |    0.040 |            NA |          NA |    NA |     0.023 |               NA |       0.067 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + fine_silt + shannon_bak                | -59.785 |     0.802 |   7.779 |      0.084 |    NA |         0.040 |    NA |       NA |            NA |          NA |    NA |     0.028 |               NA |       0.048 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + ammonium + shannon_bak                 | -59.642 |     0.945 |   7.831 |      0.086 |    NA |         0.029 |    NA |    0.026 |            NA |          NA |    NA |        NA |               NA |       0.044 |                NA |       NA |
| BrayDistance \~ p_h\_water + ammonium + finesilt_and_clay + shannon_bak             | -59.622 |     0.964 |   6.617 |      0.110 |    NA |            NA |    NA |    0.042 |            NA |          NA |    NA |        NA |               NA |       0.060 |             0.029 |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + fine_silt + shannon_bak           | -59.615 |     0.971 |   8.209 |      0.067 | 0.031 |         0.038 |    NA |       NA |            NA |          NA |    NA |     0.023 |               NA |       0.053 |                NA |       NA |
| BrayDistance \~ p_h\_water + ammonium + coarse_silt_sand + shannon_bak              | -59.594 |     0.993 |   6.933 |      0.105 |    NA |            NA |    NA |    0.040 |            NA |          NA |    NA |        NA |            0.028 |       0.057 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + ammonium + oc_beregnet + shannon_bak              | -59.572 |     1.015 |   8.574 |      0.063 | 0.046 |            NA |    NA |    0.030 |            NA |       0.020 |    NA |        NA |               NA |       0.063 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + ammonium + finesilt_and_clay + shannon_bak        | -59.517 |     1.070 |   6.774 |      0.092 | 0.032 |            NA |    NA |    0.040 |            NA |          NA |    NA |        NA |               NA |       0.066 |             0.019 |       NA |
| BrayDistance \~ p_h\_water + water_content + nitrat_nitrit + shannon_bak            | -59.453 |     1.133 |   7.760 |      0.083 |    NA |         0.043 |    NA |       NA |         0.024 |          NA |    NA |        NA |               NA |       0.048 |                NA |       NA |
| BrayDistance \~ p_h\_water + ammonium + nitrat_nitrit + shannon_bak                 | -59.418 |     1.169 |   6.502 |      0.112 |    NA |            NA |    NA |    0.042 |         0.026 |          NA |    NA |        NA |               NA |       0.066 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + ammonium + coarse_silt_sand + shannon_bak         | -59.416 |     1.170 |   7.332 |      0.076 | 0.031 |            NA |    NA |    0.040 |            NA |          NA |    NA |        NA |            0.018 |       0.061 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + finesilt_and_clay + shannon_bak        | -59.401 |     1.185 |   7.817 |      0.086 |    NA |         0.039 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.048 |             0.023 |       NA |
| BrayDistance \~ p_h\_water + shannon_bak                                            | -59.392 |     1.195 |   5.305 |      0.121 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.071 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + ammonium + nitrat_nitrit + shannon_bak            | -59.344 |     1.243 |   6.704 |      0.097 | 0.033 |            NA |    NA |    0.038 |         0.017 |          NA |    NA |        NA |               NA |       0.071 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + nitrat_nitrit + shannon_bak       | -59.339 |     1.248 |   8.227 |      0.066 | 0.032 |         0.038 |    NA |       NA |         0.019 |          NA |    NA |        NA |               NA |       0.053 |                NA |       NA |
| BrayDistance \~ p_h\_water + fine_silt + shannon_bak                                | -59.335 |     1.251 |   6.160 |      0.119 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.033 |               NA |       0.065 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + fine_silt + shannon_bak                           | -59.311 |     1.276 |   6.375 |      0.105 | 0.033 |            NA |    NA |       NA |            NA |          NA |    NA |     0.021 |               NA |       0.071 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + finesilt_and_clay + shannon_bak   | -59.243 |     1.344 |   8.201 |      0.066 | 0.032 |         0.036 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.052 |             0.018 |       NA |
| BrayDistance \~ p_h\_water + water_content                                          | -59.241 |     1.345 |   1.010 |      0.158 |    NA |         0.069 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + oc_beregnet + shannon_bak              | -59.197 |     1.390 |   8.086 |      0.082 |    NA |         0.037 |    NA |       NA |            NA |       0.020 |    NA |        NA |               NA |       0.047 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + wr + oc_beregnet + shannon_bak                    | -59.140 |     1.447 |   9.981 |      0.066 | 0.047 |            NA | 0.024 |       NA |            NA |       0.022 |    NA |        NA |               NA |       0.061 |                NA |       NA |
| BrayDistance \~ p_h\_water + coarse_silt_sand + shannon_bak                         | -59.128 |     1.459 |   6.685 |      0.108 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.030 |       0.058 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + coarse_silt_sand + shannon_bak         | -59.127 |     1.460 |   7.852 |      0.087 |    NA |         0.034 |    NA |       NA |            NA |          NA |    NA |        NA |            0.019 |       0.048 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + ammonium + clay + shannon_bak                     | -59.107 |     1.480 |   6.666 |      0.093 | 0.036 |            NA |    NA |    0.039 |            NA |          NA | 0.014 |        NA |               NA |       0.067 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + finesilt_and_clay + shannon_bak                   | -59.059 |     1.528 |   6.351 |      0.104 | 0.035 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.071 |             0.018 |       NA |
| BrayDistance \~ p_h\_water + ec + water_content                                     | -59.038 |     1.549 |   2.269 |      0.062 | 0.031 |         0.060 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + finesilt_and_clay + shannon_bak                        | -59.019 |     1.568 |   6.108 |      0.118 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.064 |             0.029 |       NA |
| BrayDistance \~ p_h\_water + ec + nitrat_nitrit + shannon_bak                       | -59.011 |     1.576 |   6.142 |      0.114 | 0.037 |            NA |    NA |       NA |         0.017 |          NA |    NA |        NA |               NA |       0.077 |                NA |       NA |
| BrayDistance \~ p_h\_water + oc_beregnet + shannon_bak                              | -59.007 |     1.579 |   7.448 |      0.087 |    NA |            NA |    NA |       NA |            NA |       0.028 |    NA |        NA |               NA |       0.055 |                NA |       NA |
| BrayDistance \~ p_h\_water + ammonium + clay + shannon_bak                          | -58.954 |     1.633 |   6.417 |      0.109 |    NA |            NA |    NA |    0.042 |            NA |          NA | 0.020 |        NA |               NA |       0.060 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + coarse_silt_sand + shannon_bak                    | -58.948 |     1.639 |   6.737 |      0.078 | 0.032 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.017 |       0.062 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + coarse_silt_sand + shannon_bak    | -58.935 |     1.651 |   8.345 |      0.063 | 0.031 |         0.034 |    NA |       NA |            NA |          NA |    NA |        NA |            0.014 |       0.052 |                NA |       NA |
| BrayDistance \~ p_h\_water + oc_beregnet + fine_silt + shannon_bak                  | -58.905 |     1.681 |   7.510 |      0.087 |    NA |            NA |    NA |       NA |            NA |       0.028 |    NA |     0.033 |               NA |       0.054 |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + ammonium                               | -58.885 |     1.702 |   2.295 |      0.161 |    NA |         0.049 |    NA |    0.029 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + clay + shannon_bak                     | -58.872 |     1.714 |   7.850 |      0.087 |    NA |         0.041 |    NA |       NA |            NA |          NA | 0.016 |        NA |               NA |       0.047 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + oc_beregnet + fine_silt + shannon_bak             | -58.853 |     1.734 |   9.134 |      0.051 | 0.033 |            NA |    NA |       NA |            NA |       0.028 |    NA |     0.020 |               NA |       0.054 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + clay + shannon_bak                | -58.833 |     1.754 |   8.234 |      0.065 | 0.033 |         0.036 |    NA |       NA |            NA |          NA | 0.013 |        NA |               NA |       0.052 |                NA |       NA |
| BrayDistance \~ p_h\_water + nitrat_nitrit + shannon_bak                            | -58.832 |     1.755 |   5.752 |      0.122 |    NA |            NA |    NA |       NA |         0.026 |          NA |    NA |        NA |               NA |       0.071 |                NA |       NA |
| BrayDistance \~ p_h\_water + wr + ammonium + shannon_bak                            | -58.808 |     1.779 |   8.636 |      0.096 |    NA |            NA | 0.018 |    0.034 |            NA |          NA |    NA |        NA |               NA |       0.066 |                NA |       NA |
| BrayDistance \~ p_h\_water + wr + shannon_bak                                       | -58.805 |     1.781 |   8.286 |      0.094 |    NA |            NA | 0.026 |       NA |            NA |          NA |    NA |        NA |               NA |       0.070 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + wr + nitrat_nitrit + shannon_bak                  | -58.766 |     1.821 |   9.320 |      0.081 | 0.042 |            NA | 0.031 |       NA |         0.017 |          NA |    NA |        NA |               NA |       0.073 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + ammonium + dexter_n + shannon_bak                 | -58.759 |     1.828 |   6.635 |      0.099 | 0.042 |            NA |    NA |    0.036 |            NA |          NA |    NA |        NA |               NA |       0.070 |                NA |    0.010 |
| BrayDistance \~ p_h\_water + water_content + fine_silt                              | -58.739 |     1.848 |   1.706 |      0.127 |    NA |         0.057 |    NA |       NA |            NA |          NA |    NA |     0.027 |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + water_content + ammonium + fine_silt + shannon_bak     | -58.737 |     1.850 |   7.856 |      0.081 |    NA |         0.018 |    NA |    0.020 |            NA |          NA |    NA |     0.022 |               NA |       0.046 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + wr                                | -58.731 |     1.856 |   4.161 |      0.056 | 0.040 |         0.056 | 0.030 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + clay + shannon_bak                                | -58.703 |     1.884 |   6.217 |      0.106 | 0.039 |            NA |    NA |       NA |            NA |          NA | 0.013 |        NA |               NA |       0.071 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + water_content + dexter_n + shannon_bak            | -58.663 |     1.924 |   8.263 |      0.065 | 0.037 |         0.035 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.051 |                NA |    0.010 |
| BrayDistance \~ p_h\_water + ec + water_content + ammonium                          | -58.651 |     1.936 |   2.510 |      0.065 | 0.031 |         0.045 |    NA |    0.029 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ p_h\_water + ammonium + oc_beregnet + shannon_bak                   | -58.637 |     1.950 |   7.531 |      0.089 |    NA |            NA |    NA |    0.029 |            NA |       0.016 |    NA |        NA |               NA |       0.057 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + dexter_n + shannon_bak                            | -58.631 |     1.956 |   6.062 |      0.115 | 0.045 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.077 |                NA |    0.012 |
| BrayDistance \~ p_h\_water + water_content + ammonium + nitrat_nitrit + shannon_bak | -58.629 |     1.957 |   7.831 |      0.078 |    NA |         0.023 |    NA |    0.023 |         0.020 |          NA |    NA |        NA |               NA |       0.044 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + nitrat_nitrit + oc_beregnet + shannon_bak         | -58.608 |     1.978 |   8.620 |      0.059 | 0.037 |            NA |    NA |       NA |         0.017 |       0.029 |    NA |        NA |               NA |       0.060 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + oc_beregnet + coarse_silt_sand + shannon_bak      | -58.605 |     1.982 |   9.130 |      0.051 | 0.035 |            NA |    NA |       NA |            NA |       0.029 |    NA |        NA |            0.017 |       0.053 |                NA |       NA |
| BrayDistance \~ p_h\_water + ec + oc_beregnet + finesilt_and_clay + shannon_bak     | -58.605 |     1.982 |   9.130 |      0.051 | 0.035 |            NA |    NA |       NA |            NA |       0.028 |    NA |        NA |               NA |       0.053 |             0.017 |       NA |
| BrayDistance \~ p_h\_water + ammonium + fine_silt + coarse_silt_sand + shannon_bak  | -58.600 |     1.986 |   9.030 |      0.093 |    NA |            NA |    NA |    0.034 |            NA |          NA |    NA |     0.021 |            0.016 |       0.055 |                NA |       NA |
| BrayDistance \~ p_h\_water + fine_silt + coarse_silt_sand + shannon_bak             | -58.589 |     1.998 |   7.768 |      0.090 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.027 |            0.024 |       0.054 |                NA |       NA |
| BrayDistance \~ p_h\_water + wr + fine_silt + shannon_bak                           | -58.589 |     1.998 |   8.891 |      0.089 |    NA |            NA | 0.024 |       NA |            NA |          NA |    NA |     0.031 |               NA |       0.063 |                NA |       NA |

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
  Vars <- Vars[Vars != "area"]
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
#> [1] "1 of 3 ready 2023-03-23 11:52:11"
#> [1] "2 of 3 ready 2023-03-23 11:53:25"
#> [1] "3 of 3 ready 2023-03-23 11:54:40"


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

This generate up to 2,358 models to evaluate, which can be downloaded as
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
  Vars <- Vars[Vars != "area"]
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
  Temp$AICc <-  try(AICc.PERMANOVA2(with(Response,adonis2(as.formula(AllForms$Form[x]), data = Response, by = "margin", strata = area)))$AICc, silent = T)
  
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
  
  Rs <- broom::tidy(with(Response,adonis2(as.formula(AllForms$Form[x]), data = Response, by = "margin"))) %>% dplyr::filter(!(term %in% c("Residual", "Total"))) %>% dplyr::select(term, R2) %>%  pivot_wider(names_from = term, values_from = R2)
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
  Vars <- Vars[Vars != "area"]
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

This generate up to 2,358 models to evaluate, which can be downloaded as
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
  Vars <- Vars[Vars != "area"]
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
  Temp$AICc <-  try(AICc.PERMANOVA2(with(Response,adonis2(as.formula(AllForms$Form[x]), data = Response, by = "margin", strata = area)))$AICc, silent = T)
  
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
  
  Rs <- broom::tidy(with(Response,adonis2(as.formula(AllForms$Form[x]), data = Response, by = "margin"))) %>% dplyr::filter(!(term %in% c("Residual", "Total"))) %>% dplyr::select(term, R2) %>%  pivot_wider(names_from = term, values_from = R2)
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
