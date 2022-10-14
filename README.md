14/10, 2022

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

<!-- README.md is generated from README.Rmd. Please edit that file -->

# 1 AICcPermanova

<!-- badges: start -->
<!-- badges: end -->

The goal of AICcPermanova is to evaluate the best models for plant
communities and bacterial communities in Denmark in order to do that we
require the following packages

``` r
library(vegan)
library(dplyr)
library(readxl)
library(parallel)
library(foreach)
library(doParallel)
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
#> [1] "1 of 3 ready 2022-10-14 07:51:06"
#> [1] "2 of 3 ready 2022-10-14 07:52:14"
#> [1] "3 of 3 ready 2022-10-14 07:53:23"


AllForms <- AllForms %>% 
  purrr::reduce(bind_rows) %>% 
  dplyr::distinct(Form, AICc, .keep_all = T)

NullMod <- data.frame(Form = paste(Dataset, "~ 1", collapse = ""), AICc = NA) %>% 
  mutate(Dataset = METADATAS[1])

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

Fs <- foreach(x = 1:nrow(AllForms), .packages = c("vegan", "dplyr", "tidyr", "readxl"), .combine = bind_rows) %dopar% {
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
  Temp$AICc <-  try(AICc.PERMANOVA2(adonis2(as.formula(AllForms$Form[x]), data = Response, by = "terms"))$AICc, silent = T)
  Rs <- broom::tidy(adonis2(as.formula(AllForms$Form[x]), data = Response, by = "terms")) %>% dplyr::filter(!(term %in% c("Residual", "Total"))) %>% dplyr::select(term, R2) %>%  pivot_wider(names_from = term, values_from = R2)
  if((x %% 100) == 0){
    sink("log.txt", append = T)
    cat(paste("finished", x, "number of models", Sys.time(), "of",  nrow(Formulas)))
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

Now we can see the top models that have a delta AICc within 2 of the
best model in table <a href="#tab:BestPLantModels">2.1</a> if expanded

| Form                                                                                                 |    AICc | DeltaAICc | habitat_type |  area | p_h\_water |    ec | water_content |    wr | ammonium | nitrat_nitrit | oc_beregnet |  clay | fine_silt | coarse_silt_sand | shannon_bak | finesilt_and_clay | dexter_n |
|:-----------------------------------------------------------------------------------------------------|--------:|----------:|-------------:|------:|-----------:|------:|--------------:|------:|---------:|--------------:|------------:|------:|----------:|-----------------:|------------:|------------------:|---------:|
| JaccardDistance \~ habitat_type                                                                      | -60.591 |     0.000 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water                                                         | -60.384 |     0.207 |        0.411 |    NA |      0.029 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec                                                                 | -60.371 |     0.220 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + fine_silt                                                     | -60.298 |     0.293 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |          NA |    NA |     0.030 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + fine_silt                                                          | -60.217 |     0.374 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.027 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + fine_silt                                             | -60.092 |     0.499 |        0.411 |    NA |      0.029 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.028 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet                                                        | -60.055 |     0.536 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |       0.025 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + oc_beregnet                                                   | -60.044 |     0.547 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |       0.027 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content + oc_beregnet                                        | -60.002 |     0.589 |        0.411 |    NA |         NA |    NA |         0.016 |    NA |       NA |            NA |       0.040 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + fine_silt                                            | -59.990 |     0.601 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |       0.025 |    NA |     0.030 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + water_content + oc_beregnet                                   | -59.972 |     0.619 |        0.411 |    NA |         NA | 0.029 |         0.019 |    NA |       NA |            NA |       0.039 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + finesilt_and_clay                                                  | -59.913 |     0.678 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.024 |       NA |
| JaccardDistance \~ habitat_type + wr                                                                 | -59.843 |     0.747 |        0.411 |    NA |         NA |    NA |            NA | 0.023 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + finesilt_and_clay                                     | -59.810 |     0.780 |        0.411 |    NA |      0.029 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.025 |       NA |
| JaccardDistance \~ habitat_type + ec + oc_beregnet + fine_silt                                       | -59.803 |     0.788 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |       0.027 |    NA |     0.028 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + finesilt_and_clay                                             | -59.789 |     0.802 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.025 |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + coarse_silt_sand                                     | -59.765 |     0.826 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |       0.025 |    NA |        NA |            0.028 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + finesilt_and_clay                                    | -59.765 |     0.826 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |       0.025 |    NA |        NA |               NA |          NA |             0.028 |       NA |
| JaccardDistance \~ habitat_type + coarse_silt_sand + finesilt_and_clay                               | -59.765 |     0.826 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.018 |          NA |             0.035 |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + coarse_silt_sand + finesilt_and_clay                 | -59.765 |     0.826 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |       0.025 |    NA |        NA |            0.028 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + wr                                                            | -59.684 |     0.906 |        0.411 |    NA |         NA | 0.029 |            NA | 0.024 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + area                                                               | -59.621 |     0.970 |        0.411 | 0.312 |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + wr + fine_silt                                                | -59.597 |     0.994 |        0.411 |    NA |         NA | 0.029 |            NA | 0.024 |       NA |            NA |          NA |    NA |     0.030 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + fine_silt + coarse_silt_sand                                       | -59.573 |     1.018 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.027 |            0.024 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + clay                                                               | -59.562 |     1.029 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |          NA | 0.020 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + oc_beregnet                                           | -59.548 |     1.043 |        0.411 |    NA |      0.029 |    NA |            NA |    NA |       NA |            NA |       0.022 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + wr + fine_silt                                                     | -59.492 |     1.098 |        0.411 |    NA |         NA |    NA |            NA | 0.023 |       NA |            NA |          NA |    NA |     0.027 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + water_content + oc_beregnet                           | -59.474 |     1.117 |        0.411 |    NA |      0.029 |    NA |         0.017 |    NA |       NA |            NA |       0.035 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + oc_beregnet + fine_silt                               | -59.474 |     1.117 |        0.411 |    NA |      0.029 |    NA |            NA |    NA |       NA |            NA |       0.022 |    NA |     0.030 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + oc_beregnet + coarse_silt_sand                                | -59.436 |     1.155 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |       0.027 |    NA |        NA |            0.025 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + oc_beregnet + finesilt_and_clay                               | -59.436 |     1.155 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |       0.027 |    NA |        NA |               NA |          NA |             0.025 |       NA |
| JaccardDistance \~ habitat_type + ec + coarse_silt_sand + finesilt_and_clay                          | -59.436 |     1.155 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.018 |          NA |             0.034 |       NA |
| JaccardDistance \~ habitat_type + ec + oc_beregnet + coarse_silt_sand + finesilt_and_clay            | -59.436 |     1.155 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |       0.027 |    NA |        NA |            0.025 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + coarse_silt_sand                                                   | -59.425 |     1.165 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.018 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + dexter_n                                                           | -59.409 |     1.182 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.018 |
| JaccardDistance \~ habitat_type + p_h\_water + wr                                                    | -59.396 |     1.194 |        0.411 |    NA |      0.029 |    NA |            NA | 0.021 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + nitrat_nitrit                                                      | -59.345 |     1.246 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |         0.017 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + ec                                                    | -59.335 |     1.256 |        0.411 |    NA |      0.029 | 0.020 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + clay                                                 | -59.332 |     1.259 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |       0.025 | 0.023 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + clay + fine_silt                                              | -59.314 |     1.277 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |          NA | 0.019 |     0.032 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + clay                                                  | -59.309 |     1.281 |        0.411 |    NA |      0.029 |    NA |            NA |    NA |       NA |            NA |          NA | 0.020 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content + oc_beregnet + fine_silt                            | -59.293 |     1.298 |        0.411 |    NA |         NA |    NA |         0.016 |    NA |       NA |            NA |       0.040 |    NA |     0.024 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + fine_silt + coarse_silt_sand                                  | -59.285 |     1.306 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |          NA |    NA |     0.030 |            0.021 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content + oc_beregnet + coarse_silt_sand                     | -59.266 |     1.325 |        0.411 |    NA |         NA |    NA |         0.016 |    NA |       NA |            NA |       0.040 |    NA |        NA |            0.024 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content + oc_beregnet + finesilt_and_clay                    | -59.266 |     1.325 |        0.411 |    NA |         NA |    NA |         0.016 |    NA |       NA |            NA |       0.040 |    NA |        NA |               NA |          NA |             0.024 |       NA |
| JaccardDistance \~ habitat_type + water_content + coarse_silt_sand + finesilt_and_clay               | -59.266 |     1.325 |        0.411 |    NA |         NA |    NA |         0.016 |    NA |       NA |            NA |          NA |    NA |        NA |            0.020 |          NA |             0.043 |       NA |
| JaccardDistance \~ habitat_type + water_content + oc_beregnet + coarse_silt_sand + finesilt_and_clay | -59.266 |     1.325 |        0.411 |    NA |         NA |    NA |         0.016 |    NA |       NA |            NA |       0.040 |    NA |        NA |            0.024 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + clay + fine_silt                                                   | -59.260 |     1.331 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |          NA | 0.020 |     0.028 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + clay                                                          | -59.258 |     1.333 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |          NA | 0.019 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + dexter_n                                                      | -59.223 |     1.368 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.019 |
| JaccardDistance \~ habitat_type + ec + water_content                                                 | -59.205 |     1.386 |        0.411 |    NA |         NA | 0.029 |         0.019 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + wr + oc_beregnet                                                   | -59.192 |     1.399 |        0.411 |    NA |         NA |    NA |            NA | 0.023 |       NA |            NA |       0.024 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + coarse_silt_sand                                      | -59.190 |     1.401 |        0.411 |    NA |      0.029 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.018 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content                                                      | -59.183 |     1.408 |        0.411 |    NA |         NA |    NA |         0.016 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + wr + oc_beregnet                                              | -59.180 |     1.411 |        0.411 |    NA |         NA | 0.029 |            NA | 0.024 |       NA |            NA |       0.026 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + oc_beregnet + coarse_silt_sand                        | -59.174 |     1.417 |        0.411 |    NA |      0.029 |    NA |            NA |    NA |       NA |            NA |       0.022 |    NA |        NA |            0.027 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + oc_beregnet + finesilt_and_clay                       | -59.174 |     1.417 |        0.411 |    NA |      0.029 |    NA |            NA |    NA |       NA |            NA |       0.022 |    NA |        NA |               NA |          NA |             0.027 |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + coarse_silt_sand + finesilt_and_clay                  | -59.174 |     1.417 |        0.411 |    NA |      0.029 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.018 |          NA |             0.031 |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + oc_beregnet + coarse_silt_sand + finesilt_and_clay    | -59.174 |     1.417 |        0.411 |    NA |      0.029 |    NA |            NA |    NA |       NA |            NA |       0.022 |    NA |        NA |            0.027 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + dexter_n                                              | -59.161 |     1.429 |        0.411 |    NA |      0.029 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.018 |
| JaccardDistance \~ habitat_type + water_content + oc_beregnet + clay                                 | -59.150 |     1.441 |        0.411 |    NA |         NA |    NA |         0.016 |    NA |       NA |            NA |       0.040 | 0.022 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + nitrat_nitrit                                                 | -59.141 |     1.450 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |         0.018 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + coarse_silt_sand                                              | -59.136 |     1.454 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.018 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content + wr + oc_beregnet                                   | -59.120 |     1.470 |        0.411 |    NA |         NA |    NA |         0.016 | 0.025 |       NA |            NA |       0.037 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + wr + finesilt_and_clay                                             | -59.120 |     1.471 |        0.411 |    NA |         NA |    NA |            NA | 0.023 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.024 |       NA |
| JaccardDistance \~ habitat_type + ec + water_content + wr + oc_beregnet                              | -59.105 |     1.486 |        0.411 |    NA |         NA | 0.029 |         0.019 | 0.024 |       NA |            NA |       0.037 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + nitrat_nitrit                                         | -59.100 |     1.491 |        0.411 |    NA |      0.029 |    NA |            NA |    NA |       NA |         0.018 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + water_content + oc_beregnet + fine_silt                       | -59.096 |     1.494 |        0.411 |    NA |         NA | 0.029 |         0.019 |    NA |       NA |            NA |       0.039 |    NA |     0.022 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + water_content                                         | -59.067 |     1.524 |        0.411 |    NA |      0.029 |    NA |         0.017 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + oc_beregnet + clay                                            | -59.059 |     1.532 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |       0.027 | 0.021 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + wr + fine_silt                                        | -59.050 |     1.541 |        0.411 |    NA |      0.029 |    NA |            NA | 0.021 |       NA |            NA |          NA |    NA |     0.028 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + wr + finesilt_and_clay                                        | -59.048 |     1.543 |        0.411 |    NA |         NA | 0.029 |            NA | 0.024 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.025 |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + clay + fine_silt                                     | -59.044 |     1.547 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |       0.025 | 0.023 |     0.028 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + clay + coarse_silt_sand                              | -59.044 |     1.547 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |       0.025 | 0.023 |        NA |            0.028 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + fine_silt + coarse_silt_sand                         | -59.044 |     1.547 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |       0.025 |    NA |     0.030 |            0.021 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + clay + fine_silt + coarse_silt_sand                                | -59.044 |     1.547 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |          NA | 0.020 |     0.028 |            0.029 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + clay + fine_silt + coarse_silt_sand                  | -59.044 |     1.547 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |       0.025 | 0.023 |     0.028 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + wr + oc_beregnet + fine_silt                                       | -59.030 |     1.560 |        0.411 |    NA |         NA |    NA |            NA | 0.023 |       NA |            NA |       0.024 |    NA |     0.030 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + fine_silt + dexter_n                                          | -59.030 |     1.561 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |          NA |    NA |     0.030 |               NA |          NA |                NA |    0.018 |
| JaccardDistance \~ habitat_type + ec + water_content + oc_beregnet + coarse_silt_sand                | -58.975 |     1.616 |        0.411 |    NA |         NA | 0.029 |         0.019 |    NA |       NA |            NA |       0.039 |    NA |        NA |            0.021 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + water_content + oc_beregnet + finesilt_and_clay               | -58.975 |     1.616 |        0.411 |    NA |         NA | 0.029 |         0.019 |    NA |       NA |            NA |       0.039 |    NA |        NA |               NA |          NA |             0.021 |       NA |
| JaccardDistance \~ habitat_type + ec + water_content + coarse_silt_sand + finesilt_and_clay          | -58.975 |     1.616 |        0.411 |    NA |         NA | 0.029 |         0.019 |    NA |       NA |            NA |          NA |    NA |        NA |            0.018 |          NA |             0.041 |       NA |
| JaccardDistance \~ habitat_type + ec + oc_beregnet + clay + fine_silt                                | -58.923 |     1.668 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |       0.027 | 0.021 |     0.030 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + oc_beregnet + clay + coarse_silt_sand                         | -58.923 |     1.668 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |       0.027 | 0.021 |        NA |            0.030 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + oc_beregnet + fine_silt + coarse_silt_sand                    | -58.923 |     1.668 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |       0.027 |    NA |     0.028 |            0.022 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + clay + fine_silt + coarse_silt_sand                           | -58.923 |     1.668 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |          NA | 0.019 |     0.032 |            0.027 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + wr + oc_beregnet + fine_silt                                  | -58.919 |     1.672 |        0.411 |    NA |         NA | 0.029 |            NA | 0.024 |       NA |            NA |       0.026 |    NA |     0.028 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + shannon_bak                                                        | -58.918 |     1.672 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.013 |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + water_content + fine_silt                                     | -58.906 |     1.685 |        0.411 |    NA |         NA | 0.029 |         0.019 |    NA |       NA |            NA |          NA |    NA |     0.028 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + water_content + oc_beregnet + clay                            | -58.900 |     1.690 |        0.411 |    NA |         NA | 0.029 |         0.019 |    NA |       NA |            NA |       0.039 | 0.020 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + fine_silt + dexter_n                                               | -58.878 |     1.713 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.027 |               NA |          NA |                NA |    0.017 |
| JaccardDistance \~ habitat_type + p_h\_water + ec + fine_silt                                        | -58.856 |     1.734 |        0.411 |    NA |      0.029 | 0.020 |            NA |    NA |       NA |            NA |          NA |    NA |     0.026 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + clay + fine_silt                                      | -58.808 |     1.783 |        0.411 |    NA |      0.029 |    NA |            NA |    NA |       NA |            NA |          NA | 0.020 |     0.026 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + ec + oc_beregnet                                      | -58.806 |     1.784 |        0.411 |    NA |      0.029 | 0.020 |            NA |    NA |       NA |            NA |       0.026 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + nitrat_nitrit + oc_beregnet                                        | -58.800 |     1.791 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |         0.017 |       0.026 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + nitrat_nitrit + fine_silt                                     | -58.782 |     1.809 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |         0.018 |          NA |    NA |     0.028 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + nitrat_nitrit + oc_beregnet                                   | -58.762 |     1.829 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |         0.018 |       0.027 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + wr + finesilt_and_clay                                | -58.722 |     1.868 |        0.411 |    NA |      0.029 |    NA |            NA | 0.021 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.024 |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + fine_silt + dexter_n                                  | -58.706 |     1.885 |        0.411 |    NA |      0.029 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.028 |               NA |          NA |                NA |    0.017 |
| JaccardDistance \~ habitat_type + nitrat_nitrit + fine_silt                                          | -58.698 |     1.893 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |         0.017 |          NA |    NA |     0.025 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + shannon_bak                                                   | -58.693 |     1.898 |        0.411 |    NA |         NA | 0.029 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.013 |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content + fine_silt                                          | -58.690 |     1.900 |        0.411 |    NA |         NA |    NA |         0.016 |    NA |       NA |            NA |          NA |    NA |     0.026 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + clay + coarse_silt_sand                                            | -58.689 |     1.902 |        0.411 |    NA |         NA |    NA |            NA |    NA |       NA |            NA |          NA | 0.020 |        NA |            0.022 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + area + fine_silt                                                   | -58.665 |     1.926 |        0.411 | 0.312 |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.021 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + wr + oc_beregnet + coarse_silt_sand                                | -58.664 |     1.927 |        0.411 |    NA |         NA |    NA |            NA | 0.023 |       NA |            NA |       0.024 |    NA |        NA |            0.026 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + wr + oc_beregnet + finesilt_and_clay                               | -58.664 |     1.927 |        0.411 |    NA |         NA |    NA |            NA | 0.023 |       NA |            NA |       0.024 |    NA |        NA |               NA |          NA |             0.026 |       NA |
| JaccardDistance \~ habitat_type + wr + coarse_silt_sand + finesilt_and_clay                          | -58.664 |     1.927 |        0.411 |    NA |         NA |    NA |            NA | 0.023 |       NA |            NA |          NA |    NA |        NA |            0.017 |          NA |             0.033 |       NA |
| JaccardDistance \~ habitat_type + wr + oc_beregnet + coarse_silt_sand + finesilt_and_clay            | -58.664 |     1.927 |        0.411 |    NA |         NA |    NA |            NA | 0.023 |       NA |            NA |       0.024 |    NA |        NA |            0.026 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + ec + finesilt_and_clay                                | -58.656 |     1.935 |        0.411 |    NA |      0.029 | 0.020 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.024 |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + water_content + oc_beregnet + fine_silt               | -58.650 |     1.941 |        0.411 |    NA |      0.029 |    NA |         0.017 |    NA |       NA |            NA |       0.035 |    NA |     0.023 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + fine_silt + coarse_silt_sand                          | -58.640 |     1.951 |        0.411 |    NA |      0.029 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.028 |            0.016 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + wr + dexter_n                                                      | -58.602 |     1.989 |        0.411 |    NA |         NA |    NA |            NA | 0.023 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.018 |
| JaccardDistance \~ habitat_type + ammonium                                                           | -58.594 |     1.997 |        0.411 |    NA |         NA |    NA |            NA |    NA |    0.009 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |

Table 2.1: Best models
