17/10, 2022

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
#> [1] "1 of 3 ready 2022-10-17 03:56:59"
#> [1] "2 of 3 ready 2022-10-17 03:58:10"
#> [1] "3 of 3 ready 2022-10-17 03:59:10"


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

As seen in table <a href="#tab:SummaryPlantPA">2.1</a> there are 113
models within 2 AICc of each other, you can see there how many times a
variable has been selected

| Variable          | Number_of_models |
|:------------------|-----------------:|
| habitat_type      |              113 |
| oc_beregnet       |               49 |
| ec                |               42 |
| fine_silt         |               36 |
| coarse_silt_sand  |               30 |
| p_h\_water        |               26 |
| water_content     |               22 |
| wr                |               20 |
| finesilt_and_clay |               19 |
| clay              |               18 |
| nitrat_nitrit     |                7 |
| dexter_n          |                7 |
| area              |                2 |
| shannon_bak       |                2 |
| ammonium          |                1 |

Table 2.1: Number of selected models were variables are present

Now we can see the top models that have a delta AICc within 2 of the
best model in table <a href="#tab:BestPLantModels">2.2</a> if expanded

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Show table of selected models for vegetation presence absence

</summary>

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

Table 2.2: Best models

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
#> [1] "1 of 3 ready 2022-10-17 04:04:36"
#> [1] "2 of 3 ready 2022-10-17 04:05:58"
#> [1] "3 of 3 ready 2022-10-17 04:07:20"


AllForms <- AllForms %>% 
  purrr::reduce(bind_rows) %>% 
  dplyr::distinct(Form, AICc, .keep_all = T)

NullMod <- data.frame(Form = paste(Dataset, "~ 1", collapse = ""), AICc = NA) %>% 
  mutate(Dataset = METADATAS[1])

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

Fs <- foreach(x = 1:nrow(AllForms), .packages = c("vegan", "dplyr", "tidyr", "readxl", "ampvis2"), .combine = bind_rows) %dopar% {
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
  Temp$AICc <-  try(AICc.PERMANOVA2(adonis2(as.formula(AllForms$Form[x]), data = Response, by = "terms"))$AICc, silent = T)
  Rs <- broom::tidy(adonis2(as.formula(AllForms$Form[x]), data = Response, by = "terms")) %>% dplyr::filter(!(term %in% c("Residual", "Total"))) %>% dplyr::select(term, R2) %>%  pivot_wider(names_from = term, values_from = R2)
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
69 models within 2 AICc of each other, you can see there how many times
a variable has been selected

| Variable          | Number_of_models |
|:------------------|-----------------:|
| habitat_type      |               69 |
| oc_beregnet       |               29 |
| water_content     |               26 |
| shannon_veg       |               20 |
| p_h\_water        |               19 |
| coarse_silt_sand  |               14 |
| wr                |               12 |
| finesilt_and_clay |                8 |
| fine_silt         |                7 |
| ec                |                6 |
| ammonium          |                6 |
| clay              |                2 |
| dexter_n          |                2 |

Table 3.1: Number of selected models were variables are present

Now we can see the top models that have a delta AICc within 2 of the
best model in table <a href="#tab:BestBacterialPAModels">3.4</a> if
expanded

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Show table of selected models for bacterial abundance

</summary>

| Form                                                                                              |    AICc | DeltaAICc | habitat_type | p_h\_water |    ec | water_content |    wr | ammonium | oc_beregnet |  clay | fine_silt | coarse_silt_sand | shannon_veg | finesilt_and_clay | dexter_n |
|:--------------------------------------------------------------------------------------------------|--------:|----------:|-------------:|-----------:|------:|--------------:|------:|---------:|------------:|------:|----------:|-----------------:|------------:|------------------:|---------:|
| BrayDistance \~ habitat_type + oc_beregnet                                                        | -86.086 |     0.000 |        0.529 |         NA |    NA |            NA |    NA |       NA |       0.028 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + oc_beregnet + shannon_veg                                          | -85.789 |     0.297 |        0.529 |         NA |    NA |            NA |    NA |       NA |       0.028 |    NA |        NA |               NA |       0.022 |                NA |       NA |
| BrayDistance \~ habitat_type                                                                      | -85.730 |     0.355 |        0.529 |         NA |    NA |            NA |    NA |       NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + p_h\_water                                                         | -85.661 |     0.425 |        0.529 |      0.024 |    NA |            NA |    NA |       NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + water_content                                                      | -85.591 |     0.495 |        0.529 |         NA |    NA |         0.024 |    NA |       NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + water_content + oc_beregnet                                        | -85.582 |     0.504 |        0.529 |         NA |    NA |         0.024 |    NA |       NA |       0.025 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + wr                                                                 | -85.547 |     0.539 |        0.529 |         NA |    NA |            NA | 0.023 |       NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + p_h\_water + oc_beregnet                                           | -85.541 |     0.545 |        0.529 |      0.024 |    NA |            NA |    NA |       NA |       0.024 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + p_h\_water + water_content                                         | -85.533 |     0.552 |        0.529 |      0.024 |    NA |         0.024 |    NA |       NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + shannon_veg                                                        | -85.331 |     0.755 |        0.529 |         NA |    NA |            NA |    NA |       NA |          NA |    NA |        NA |               NA |       0.021 |                NA |       NA |
| BrayDistance \~ habitat_type + water_content + wr                                                 | -85.283 |     0.803 |        0.529 |         NA |    NA |         0.024 | 0.022 |       NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + p_h\_water + shannon_veg                                           | -85.249 |     0.837 |        0.529 |      0.024 |    NA |            NA |    NA |       NA |          NA |    NA |        NA |               NA |       0.021 |                NA |       NA |
| BrayDistance \~ habitat_type + p_h\_water + oc_beregnet + shannon_veg                             | -85.182 |     0.903 |        0.529 |      0.024 |    NA |            NA |    NA |       NA |       0.024 |    NA |        NA |               NA |       0.022 |                NA |       NA |
| BrayDistance \~ habitat_type + wr + oc_beregnet                                                   | -85.149 |     0.937 |        0.529 |         NA |    NA |            NA | 0.023 |       NA |       0.021 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + water_content + oc_beregnet + shannon_veg                          | -85.130 |     0.956 |        0.529 |         NA |    NA |         0.024 |    NA |       NA |       0.025 |    NA |        NA |               NA |       0.021 |                NA |       NA |
| BrayDistance \~ habitat_type + wr + shannon_veg                                                   | -85.121 |     0.965 |        0.529 |         NA |    NA |            NA | 0.023 |       NA |          NA |    NA |        NA |               NA |       0.021 |                NA |       NA |
| BrayDistance \~ habitat_type + water_content + shannon_veg                                        | -85.047 |     1.039 |        0.529 |         NA |    NA |         0.024 |    NA |       NA |          NA |    NA |        NA |               NA |       0.020 |                NA |       NA |
| BrayDistance \~ habitat_type + p_h\_water + water_content + shannon_veg                           | -84.938 |     1.148 |        0.529 |      0.024 |    NA |         0.024 |    NA |       NA |          NA |    NA |        NA |               NA |       0.020 |                NA |       NA |
| BrayDistance \~ habitat_type + coarse_silt_sand                                                   | -84.869 |     1.216 |        0.529 |         NA |    NA |            NA |    NA |       NA |          NA |    NA |        NA |            0.017 |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + p_h\_water + water_content + oc_beregnet                           | -84.797 |     1.289 |        0.529 |      0.024 |    NA |         0.024 |    NA |       NA |       0.019 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + ec + oc_beregnet                                                   | -84.719 |     1.367 |        0.529 |         NA | 0.015 |            NA |    NA |       NA |       0.026 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + ammonium                                                           | -84.683 |     1.403 |        0.529 |         NA |    NA |            NA |    NA |    0.016 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + water_content + coarse_silt_sand                                   | -84.656 |     1.430 |        0.529 |         NA |    NA |         0.024 |    NA |       NA |          NA |    NA |        NA |            0.017 |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + ec + water_content                                                 | -84.649 |     1.437 |        0.529 |         NA | 0.015 |         0.026 |    NA |       NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + p_h\_water + wr                                                    | -84.635 |     1.451 |        0.529 |      0.024 |    NA |            NA | 0.016 |       NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + water_content + wr + oc_beregnet                                   | -84.589 |     1.497 |        0.529 |         NA |    NA |         0.024 | 0.022 |       NA |       0.019 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + ec                                                                 | -84.556 |     1.530 |        0.529 |         NA | 0.015 |            NA |    NA |       NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + oc_beregnet + fine_silt                                            | -84.527 |     1.559 |        0.529 |         NA |    NA |            NA |    NA |       NA |       0.028 |    NA |     0.012 |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + oc_beregnet + coarse_silt_sand                                     | -84.491 |     1.594 |        0.529 |         NA |    NA |            NA |    NA |       NA |       0.028 |    NA |        NA |            0.011 |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + oc_beregnet + finesilt_and_clay                                    | -84.491 |     1.594 |        0.529 |         NA |    NA |            NA |    NA |       NA |       0.028 |    NA |        NA |               NA |          NA |             0.011 |       NA |
| BrayDistance \~ habitat_type + coarse_silt_sand + finesilt_and_clay                               | -84.491 |     1.594 |        0.529 |         NA |    NA |            NA |    NA |       NA |          NA |    NA |        NA |            0.017 |          NA |             0.022 |       NA |
| BrayDistance \~ habitat_type + oc_beregnet + coarse_silt_sand + finesilt_and_clay                 | -84.491 |     1.594 |        0.529 |         NA |    NA |            NA |    NA |       NA |       0.028 |    NA |        NA |            0.011 |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + water_content + oc_beregnet + fine_silt                            | -84.464 |     1.622 |        0.529 |         NA |    NA |         0.024 |    NA |       NA |       0.025 |    NA |     0.016 |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + oc_beregnet + clay                                                 | -84.448 |     1.638 |        0.529 |         NA |    NA |            NA |    NA |       NA |       0.028 | 0.011 |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + coarse_silt_sand + shannon_veg                                     | -84.442 |     1.644 |        0.529 |         NA |    NA |            NA |    NA |       NA |          NA |    NA |        NA |            0.017 |       0.021 |                NA |       NA |
| BrayDistance \~ habitat_type + ec + water_content + oc_beregnet                                   | -84.405 |     1.681 |        0.529 |         NA | 0.015 |         0.026 |    NA |       NA |       0.023 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + water_content + wr + shannon_veg                                   | -84.403 |     1.683 |        0.529 |         NA |    NA |         0.024 | 0.022 |       NA |          NA |    NA |        NA |               NA |       0.018 |                NA |       NA |
| BrayDistance \~ habitat_type + wr + oc_beregnet + shannon_veg                                     | -84.385 |     1.700 |        0.529 |         NA |    NA |            NA | 0.023 |       NA |       0.021 |    NA |        NA |               NA |       0.019 |                NA |       NA |
| BrayDistance \~ habitat_type + ec + oc_beregnet + shannon_veg                                     | -84.370 |     1.716 |        0.529 |         NA | 0.015 |            NA |    NA |       NA |       0.026 |    NA |        NA |               NA |       0.022 |                NA |       NA |
| BrayDistance \~ habitat_type + water_content + fine_silt                                          | -84.367 |     1.719 |        0.529 |         NA |    NA |         0.024 |    NA |       NA |          NA |    NA |     0.015 |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + p_h\_water + water_content + wr                                    | -84.352 |     1.734 |        0.529 |      0.024 |    NA |         0.024 | 0.015 |       NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + p_h\_water + coarse_silt_sand                                      | -84.348 |     1.738 |        0.529 |      0.024 |    NA |            NA |    NA |       NA |          NA |    NA |        NA |            0.014 |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + finesilt_and_clay                                                  | -84.339 |     1.746 |        0.529 |         NA |    NA |            NA |    NA |       NA |          NA |    NA |        NA |               NA |          NA |             0.013 |       NA |
| BrayDistance \~ habitat_type + clay                                                               | -84.328 |     1.757 |        0.529 |         NA |    NA |            NA |    NA |       NA |          NA | 0.013 |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + fine_silt                                                          | -84.314 |     1.771 |        0.529 |         NA |    NA |            NA |    NA |       NA |          NA |    NA |     0.012 |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + p_h\_water + ammonium                                              | -84.307 |     1.779 |        0.529 |      0.024 |    NA |            NA |    NA |    0.013 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + p_h\_water + water_content + oc_beregnet + shannon_veg             | -84.306 |     1.780 |        0.529 |      0.024 |    NA |         0.024 |    NA |       NA |       0.019 |    NA |        NA |               NA |       0.021 |                NA |       NA |
| BrayDistance \~ habitat_type + dexter_n                                                           | -84.266 |     1.820 |        0.529 |         NA |    NA |            NA |    NA |       NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.012 |
| BrayDistance \~ habitat_type + water_content + ammonium                                           | -84.265 |     1.820 |        0.529 |         NA |    NA |         0.024 |    NA |    0.014 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + p_h\_water + water_content + fine_silt                             | -84.244 |     1.842 |        0.529 |      0.024 |    NA |         0.024 |    NA |       NA |          NA |    NA |     0.014 |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + p_h\_water + wr + oc_beregnet                                      | -84.236 |     1.850 |        0.529 |      0.024 |    NA |            NA | 0.016 |       NA |       0.022 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + water_content + finesilt_and_clay                                  | -84.235 |     1.850 |        0.529 |         NA |    NA |         0.024 |    NA |       NA |          NA |    NA |        NA |               NA |          NA |             0.013 |       NA |
| BrayDistance \~ habitat_type + wr + ammonium                                                      | -84.208 |     1.878 |        0.529 |         NA |    NA |            NA | 0.023 |    0.014 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + water_content + oc_beregnet + coarse_silt_sand                     | -84.206 |     1.880 |        0.529 |         NA |    NA |         0.024 |    NA |       NA |       0.025 |    NA |        NA |            0.014 |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + water_content + oc_beregnet + finesilt_and_clay                    | -84.206 |     1.880 |        0.529 |         NA |    NA |         0.024 |    NA |       NA |       0.025 |    NA |        NA |               NA |          NA |             0.014 |       NA |
| BrayDistance \~ habitat_type + water_content + coarse_silt_sand + finesilt_and_clay               | -84.206 |     1.880 |        0.529 |         NA |    NA |         0.024 |    NA |       NA |          NA |    NA |        NA |            0.017 |          NA |             0.021 |       NA |
| BrayDistance \~ habitat_type + water_content + oc_beregnet + coarse_silt_sand + finesilt_and_clay | -84.206 |     1.880 |        0.529 |         NA |    NA |         0.024 |    NA |       NA |       0.025 |    NA |        NA |            0.014 |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + ammonium + shannon_veg                                             | -84.175 |     1.911 |        0.529 |         NA |    NA |            NA |    NA |    0.016 |          NA |    NA |        NA |               NA |       0.021 |                NA |       NA |
| BrayDistance \~ habitat_type + oc_beregnet + fine_silt + shannon_veg                              | -84.170 |     1.916 |        0.529 |         NA |    NA |            NA |    NA |       NA |       0.028 |    NA |     0.012 |               NA |       0.022 |                NA |       NA |
| BrayDistance \~ habitat_type + p_h\_water + ec + water_content                                    | -84.151 |     1.934 |        0.529 |      0.024 | 0.011 |         0.026 |    NA |       NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + p_h\_water + dexter_n                                              | -84.135 |     1.951 |        0.529 |      0.024 |    NA |            NA |    NA |       NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.012 |
| BrayDistance \~ habitat_type + p_h\_water + fine_silt                                             | -84.124 |     1.962 |        0.529 |      0.024 |    NA |            NA |    NA |       NA |          NA |    NA |     0.012 |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + p_h\_water + water_content + coarse_silt_sand                      | -84.114 |     1.972 |        0.529 |      0.024 |    NA |         0.024 |    NA |       NA |          NA |    NA |        NA |            0.013 |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + oc_beregnet + coarse_silt_sand + shannon_veg                       | -84.104 |     1.981 |        0.529 |         NA |    NA |            NA |    NA |       NA |       0.028 |    NA |        NA |            0.011 |       0.022 |                NA |       NA |
| BrayDistance \~ habitat_type + oc_beregnet + finesilt_and_clay + shannon_veg                      | -84.104 |     1.981 |        0.529 |         NA |    NA |            NA |    NA |       NA |       0.028 |    NA |        NA |               NA |       0.022 |             0.011 |       NA |
| BrayDistance \~ habitat_type + coarse_silt_sand + finesilt_and_clay + shannon_veg                 | -84.104 |     1.981 |        0.529 |         NA |    NA |            NA |    NA |       NA |          NA |    NA |        NA |            0.017 |       0.022 |             0.022 |       NA |
| BrayDistance \~ habitat_type + oc_beregnet + coarse_silt_sand + finesilt_and_clay + shannon_veg   | -84.104 |     1.981 |        0.529 |         NA |    NA |            NA |    NA |       NA |       0.028 |    NA |        NA |            0.011 |       0.022 |                NA |       NA |
| BrayDistance \~ habitat_type + ammonium + oc_beregnet                                             | -84.102 |     1.983 |        0.529 |         NA |    NA |            NA |    NA |    0.016 |       0.020 |    NA |        NA |               NA |          NA |                NA |       NA |
| BrayDistance \~ habitat_type + p_h\_water + wr + shannon_veg                                      | -84.086 |     2.000 |        0.529 |      0.024 |    NA |            NA | 0.016 |       NA |          NA |    NA |        NA |               NA |       0.020 |                NA |       NA |

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
#> [1] "1 of 3 ready 2022-10-17 05:03:21"
#> [1] "2 of 3 ready 2022-10-17 05:04:43"
#> [1] "3 of 3 ready 2022-10-17 05:06:05"


AllForms <- AllForms %>% 
  purrr::reduce(bind_rows) %>% 
  dplyr::distinct(Form, AICc, .keep_all = T)

NullMod <- data.frame(Form = paste(Dataset, "~ 1", collapse = ""), AICc = NA) %>% 
  mutate(Dataset = METADATAS[1])

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

Fs <- foreach(x = 1:nrow(AllForms), .packages = c("vegan", "dplyr", "tidyr", "readxl", "ampvis2"), .combine = bind_rows) %dopar% {
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
  Temp$AICc <-  try(AICc.PERMANOVA2(adonis2(as.formula(AllForms$Form[x]), data = Response, by = "terms"))$AICc, silent = T)
  Rs <- broom::tidy(adonis2(as.formula(AllForms$Form[x]), data = Response, by = "terms")) %>% dplyr::filter(!(term %in% c("Residual", "Total"))) %>% dplyr::select(term, R2) %>%  pivot_wider(names_from = term, values_from = R2)
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

As seen in table <a href="#tab:SummaryBacterialPA">3.3</a> there are 43
models within 2 AICc of each other, you can see there how many times a
variable has been selected

| Variable          | Number_of_models |
|:------------------|-----------------:|
| habitat_type      |               43 |
| water_content     |               14 |
| oc_beregnet       |               13 |
| coarse_silt_sand  |                8 |
| shannon_veg       |                7 |
| p_h\_water        |                5 |
| wr                |                4 |
| fine_silt         |                4 |
| finesilt_and_clay |                4 |
| ec                |                3 |
| nitrat_nitrit     |                3 |
| clay              |                3 |
| ammonium          |                2 |
| dexter_n          |                2 |

Table 3.3: Number of selected models were variables are present

Now we can see the top models that have a delta AICc within 2 of the
best model in table <a href="#tab:BestBacterialPAModels">3.4</a> if
expanded

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Show table of selected models for bacterial presence

</summary>

| Form                                                                                 |    AICc | DeltaAICc | habitat_type | p_h\_water |    ec | water_content |    wr | ammonium | nitrat_nitrit | oc_beregnet |  clay | fine_silt | coarse_silt_sand | shannon_veg | finesilt_and_clay | dexter_n |
|:-------------------------------------------------------------------------------------|--------:|----------:|-------------:|-----------:|------:|--------------:|------:|---------:|--------------:|------------:|------:|----------:|-----------------:|------------:|------------------:|---------:|
| JaccardDistance \~ habitat_type + water_content                                      | -99.586 |     0.000 |        0.577 |         NA |    NA |         0.023 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type                                                      | -99.502 |     0.084 |        0.577 |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet                                        | -99.386 |     0.200 |        0.577 |         NA |    NA |            NA |    NA |       NA |            NA |       0.021 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + water_content                         | -98.820 |     0.766 |        0.577 |      0.016 |    NA |         0.023 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content + oc_beregnet                        | -98.786 |     0.800 |        0.577 |         NA |    NA |         0.023 |    NA |       NA |            NA |       0.016 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + wr                                                 | -98.751 |     0.835 |        0.577 |         NA |    NA |            NA | 0.016 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content + wr                                 | -98.719 |     0.867 |        0.577 |         NA |    NA |         0.023 | 0.016 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water                                         | -98.718 |     0.868 |        0.577 |      0.016 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + coarse_silt_sand                                   | -98.586 |     1.000 |        0.577 |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.015 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + shannon_veg                                        | -98.582 |     1.004 |        0.577 |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.015 |                NA |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + shannon_veg                          | -98.532 |     1.054 |        0.577 |         NA |    NA |            NA |    NA |       NA |            NA |       0.021 |    NA |        NA |               NA |       0.016 |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content + shannon_veg                        | -98.445 |     1.141 |        0.577 |         NA |    NA |         0.023 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.014 |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + oc_beregnet                           | -98.393 |     1.193 |        0.577 |      0.016 |    NA |            NA |    NA |       NA |            NA |       0.020 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content + fine_silt                          | -98.365 |     1.221 |        0.577 |         NA |    NA |         0.023 |    NA |       NA |            NA |          NA |    NA |     0.013 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + fine_silt                                          | -98.333 |     1.253 |        0.577 |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.013 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ammonium                                           | -98.313 |     1.274 |        0.577 |         NA |    NA |            NA |    NA |    0.013 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content + coarse_silt_sand                   | -98.261 |     1.325 |        0.577 |         NA |    NA |         0.023 |    NA |       NA |            NA |          NA |    NA |        NA |            0.012 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + water_content                                 | -98.242 |     1.344 |        0.577 |         NA | 0.011 |         0.024 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + nitrat_nitrit                                      | -98.166 |     1.420 |        0.577 |         NA |    NA |            NA |    NA |       NA |         0.012 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content + ammonium                           | -98.162 |     1.424 |        0.577 |         NA |    NA |         0.023 |    NA |    0.011 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + finesilt_and_clay                                  | -98.153 |     1.433 |        0.577 |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.012 |       NA |
| JaccardDistance \~ habitat_type + wr + oc_beregnet                                   | -98.112 |     1.474 |        0.577 |         NA |    NA |            NA | 0.016 |       NA |            NA |       0.018 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content + nitrat_nitrit                      | -98.085 |     1.501 |        0.577 |         NA |    NA |         0.023 |    NA |       NA |         0.011 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content + finesilt_and_clay                  | -98.080 |     1.506 |        0.577 |         NA |    NA |         0.023 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.011 |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + fine_silt                            | -98.078 |     1.508 |        0.577 |         NA |    NA |            NA |    NA |       NA |            NA |       0.021 |    NA |     0.012 |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + ec                                                 | -98.071 |     1.515 |        0.577 |         NA | 0.011 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + clay                                               | -98.016 |     1.570 |        0.577 |         NA |    NA |            NA |    NA |       NA |            NA |          NA | 0.011 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + nitrat_nitrit + oc_beregnet                        | -97.967 |     1.619 |        0.577 |         NA |    NA |            NA |    NA |       NA |         0.012 |       0.021 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content + clay                               | -97.939 |     1.647 |        0.577 |         NA |    NA |         0.023 |    NA |       NA |            NA |          NA | 0.010 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + wr + shannon_veg                                   | -97.938 |     1.648 |        0.577 |         NA |    NA |            NA | 0.016 |       NA |            NA |          NA |    NA |        NA |               NA |       0.016 |                NA |       NA |
| JaccardDistance \~ habitat_type + ec + oc_beregnet                                   | -97.849 |     1.737 |        0.577 |         NA | 0.011 |            NA |    NA |       NA |            NA |       0.021 |    NA |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + coarse_silt_sand                     | -97.823 |     1.764 |        0.577 |         NA |    NA |            NA |    NA |       NA |            NA |       0.021 |    NA |        NA |            0.010 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + finesilt_and_clay                    | -97.823 |     1.764 |        0.577 |         NA |    NA |            NA |    NA |       NA |            NA |       0.021 |    NA |        NA |               NA |          NA |             0.010 |       NA |
| JaccardDistance \~ habitat_type + coarse_silt_sand + finesilt_and_clay               | -97.823 |     1.764 |        0.577 |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.015 |          NA |             0.017 |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + coarse_silt_sand + finesilt_and_clay | -97.823 |     1.764 |        0.577 |         NA |    NA |            NA |    NA |       NA |            NA |       0.021 |    NA |        NA |            0.010 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + fine_silt + coarse_silt_sand                       | -97.817 |     1.769 |        0.577 |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.013 |            0.019 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + shannon_veg                           | -97.734 |     1.852 |        0.577 |      0.016 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.015 |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content + dexter_n                           | -97.662 |     1.924 |        0.577 |         NA |    NA |         0.023 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.008 |
| JaccardDistance \~ habitat_type + dexter_n                                           | -97.658 |     1.928 |        0.577 |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.008 |
| JaccardDistance \~ habitat_type + coarse_silt_sand + shannon_veg                     | -97.651 |     1.935 |        0.577 |         NA |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.015 |       0.015 |                NA |       NA |
| JaccardDistance \~ habitat_type + oc_beregnet + clay                                 | -97.643 |     1.943 |        0.577 |         NA |    NA |            NA |    NA |       NA |            NA |       0.021 | 0.009 |        NA |               NA |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + p_h\_water + coarse_silt_sand                      | -97.640 |     1.947 |        0.577 |      0.016 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.014 |          NA |                NA |       NA |
| JaccardDistance \~ habitat_type + water_content + oc_beregnet + shannon_veg          | -97.633 |     1.953 |        0.577 |         NA |    NA |         0.023 |    NA |       NA |            NA |       0.016 |    NA |        NA |               NA |       0.014 |                NA |       NA |

Table 3.4: Best models for bacterial presence absence datasets

</details>
