04/10, 2022

- <a href="#1-aiccpermanova" id="toc-1-aiccpermanova">1 AICcPermanova</a>
- <a href="#2-vegetation-models" id="toc-2-vegetation-models">2 Vegetation
  models</a>
  - <a href="#21-generation-of-all-possible-models"
    id="toc-21-generation-of-all-possible-models">2.1 Generation of all
    possible models</a>
  - <a href="#22-model-fitting" id="toc-22-model-fitting">2.2 Model
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

## 2.1 Generation of all possible models

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
  Dataset <- "vegetation_data_no_ID"
  Response = env.data
  
  Forms <- list()
  
  Models <- for(i in 1:floor(nrow(env.data)/10)){
    Test <- combn(Vars, i, simplify = F)
    cl <- makeCluster(21)
    registerDoParallel(cl)
    Formulas <- foreach(j = 1:length(Test), .combine = "rbind", .packages = c("dplyr")) %dopar% {
      Dataset <- "vegetation_data_no_ID"
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
  
  Dataset <- "vegetation_data_no_ID"
}
#> [1] "1 of 3 ready 2022-10-03 13:59:32"
#> [1] "2 of 3 ready 2022-10-03 14:00:42"
#> [1] "3 of 3 ready 2022-10-03 14:01:52"


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

## 2.2 Model fitting

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
  
  
  # Removing columns from env.data that is not used in the analysis 
  
  env.data = subset(meta.data, select = -c(order))
  
  
  
  env.data <- env.data %>% tidyr::drop_na()  
  
  Vars <- colnames(env.data)
  Dataset <- "vegetation_data_no_ID"
  
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

| Form                                                                        |    AICc | DeltaAICc | habitat_type |  area | water_content | oc_beregnet | fine_silt | coarse_silt_sand |
|:----------------------------------------------------------------------------|--------:|----------:|-------------:|------:|--------------:|------------:|----------:|-----------------:|
| vegetation_data_no_ID \~ habitat_type + area + water_content + oc_beregnet  | -89.995 |     0.000 |         0.53 | 0.285 |         0.005 |       0.032 |        NA |               NA |
| vegetation_data_no_ID \~ habitat_type + area + fine_silt + coarse_silt_sand | -89.691 |     0.305 |         0.53 | 0.285 |            NA |          NA |     0.019 |            0.017 |
| vegetation_data_no_ID \~ habitat_type + area + fine_silt                    | -89.291 |     0.704 |         0.53 | 0.285 |            NA |          NA |     0.019 |               NA |
| vegetation_data_no_ID \~ habitat_type + area + oc_beregnet                  | -88.839 |     1.156 |         0.53 | 0.285 |            NA |       0.018 |        NA |               NA |
| vegetation_data_no_ID \~ habitat_type + area                                | -88.676 |     1.319 |         0.53 | 0.285 |            NA |          NA |        NA |               NA |
| vegetation_data_no_ID \~ habitat_type + area + oc_beregnet + fine_silt      | -88.418 |     1.577 |         0.53 | 0.285 |            NA |       0.018 |     0.015 |               NA |

Table 2.1: Best models
