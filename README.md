23/04, 2023

- <a href="#1-aiccpermanova" id="toc-1-aiccpermanova">1 AICcPermanova</a>
- <a href="#2-vegetation-models" id="toc-2-vegetation-models">2 Vegetation
  models</a>
  - <a href="#21-presence-absence-data"
    id="toc-21-presence-absence-data">2.1 Presence-absence data</a>
    - <a href="#211-generation-of-all-possible-models"
      id="toc-211-generation-of-all-possible-models">2.1.1 Generation of all
      possible models</a>
  - <a href="#22-avoiding-colinearity" id="toc-22-avoiding-colinearity">2.2
    Avoiding colinearity</a>
    - <a href="#221-model-fitting" id="toc-221-model-fitting">2.2.1 Model
      fitting</a>
  - <a href="#23-abundance-data" id="toc-23-abundance-data">2.3 Abundance
    data</a>
    - <a href="#231-model-fitting" id="toc-231-model-fitting">2.3.1 Model
      fitting</a>
- <a href="#3-baterial-models" id="toc-3-baterial-models">3 Baterial
  models</a>
  - <a href="#31-abundance-data" id="toc-31-abundance-data">3.1 Abundance
    data</a>
    - <a href="#311-generation-of-all-possible-models"
      id="toc-311-generation-of-all-possible-models">3.1.1 Generation of all
      possible models</a>
  - <a href="#32-avoiding-colinearity" id="toc-32-avoiding-colinearity">3.2
    Avoiding colinearity</a>
    - <a href="#321-model-fitting" id="toc-321-model-fitting">3.2.1 Model
      fitting</a>

<!-- README.md is generated from README.Rmd. Please edit that file -->

# 1 AICcPermanova

<!-- badges: start -->
<!-- badges: end -->

The goal of AICcPermanova is to evaluate the best models for plant
communities and bacterial communities in Denmark in order to do that we
require the following packages

``` r
library(ampvis2)
library(readxl)
library(tidyverse)
library(AICcPermanova)
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
METADATAS <- list.files(pattern = "PERMANOVA_VEGETATION_", full.names = T) |> 
  purrr::map(read_excel) |> 
  purrr::reduce(full_join) |> 
  janitor::clean_names()

AllForms <- AICcPermanova::make_models(vars = c("p_h_water", "ec", "water_content", "wr", "ammonium", "nitrat_nitrit", 
"oc_beregnet", "clay", "fine_silt", "coarse_silt_sand", "shannon_bak", 
"finesilt_and_clay", "dexter_n")
, ncores = 21, k = floor(nrow(METADATAS)/10))

saveRDS(AllForms, "AllForms.rds")
openxlsx::write.xlsx(AllForms, "AllForms.xlsx")
```

</details>

This generate up to 2,380 models to evaluate, which can be downloaded as
an excel file
[here](https://github.com/Sustainscapes/AICcPermanova/raw/master/AllForms.xlsx)
an rds
[here](https://github.com/Sustainscapes/AICcPermanova/blob/master/AllForms.rds).

## 2.2 Avoiding colinearity

Now before the model fitting, we will filter out models with high degree
of multicolinearity using VIF with 5 as a cut off

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

VIF filter

</summary>

``` r
Filtered <- AICcPermanova::filter_vif(all_forms = AllForms , env_data = METADATAS, ncores = 21)
```

</details>

This is then narrowed down to 1,351 models to evaluate

### 2.2.1 Model fitting

Then in the following code each model is fitted and AICc is calculated
to order it

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Model fitting code

</summary>

``` r
vegetation_data = read_excel("Presence_absence_vegetation_AC_Danielsen.xlsx")%>% 
    janitor::clean_names()

vegetation_data_no_ID = subset(vegetation_data, select = -plot)

FS <- AICcPermanova::fit_models(all_forms = Filtered, 
                 env_data = METADATAS, 
                 veg_data = vegetation_data_no_ID, 
                 ncores = 21, 
                 method = "jaccard",
                logfile = "new.txt",
                strata = "area"
                 )


saveRDS(FS, "FS.rds")
openxlsx::write.xlsx(FS, "FS.xlsx")
```

</details>

As seen in table <a href="#tab:SummaryPlantPA">2.1</a> there are 190
models within 2 AICc where the max VIF is lower or equal than 6 of each
other, you can see there how many times a variable has been selected

| Variable          | Full_Akaike_Adjusted_RSq | Number_of_models |
|:------------------|-------------------------:|-----------------:|
| p_h\_water        |                    0.050 |              125 |
| ec                |                    0.033 |              122 |
| water_content     |                    0.027 |              110 |
| fine_silt         |                    0.009 |               49 |
| oc_beregnet       |                    0.008 |               48 |
| wr                |                    0.006 |               47 |
| ammonium          |                    0.006 |               44 |
| shannon_bak       |                    0.010 |               42 |
| finesilt_and_clay |                    0.006 |               35 |
| coarse_silt_sand  |                    0.005 |               32 |
| clay              |                    0.004 |               31 |
| nitrat_nitrit     |                    0.001 |               13 |
| dexter_n          |                    0.001 |               13 |

Table 2.1: Number of selected models were variables are present and
their Akaike Weighted R squared for the Marginal effect of the terms

Now we can see the top models that have a delta AICc within 2 of the
best model in table <a href="#tab:BestPLantModels">2.2</a> if expanded

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Show table of selected models for vegetation presence absence

</summary>

| form                                                                            | max_vif |    AICc |   k |   N | p_h\_water |    ec | water_content |    wr | ammonium | nitrat_nitrit | oc_beregnet |  clay | fine_silt | coarse_silt_sand | shannon_bak | finesilt_and_clay | dexter_n | DeltaAICc | AICWeight |
|:--------------------------------------------------------------------------------|--------:|--------:|----:|----:|-----------:|------:|--------------:|------:|---------:|--------------:|------------:|------:|----------:|-----------------:|------------:|------------------:|---------:|----------:|----------:|
| Distance \~ p_h\_water + ec + water_content                                     |   2.269 | -52.326 |   4 |  52 |      0.056 | 0.045 |         0.048 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     0.000 |     0.010 |
| Distance \~ p_h\_water + ec + water_content + oc_beregnet                       |   3.146 | -52.069 |   5 |  52 |      0.050 | 0.040 |         0.058 |    NA |       NA |            NA |       0.033 |    NA |        NA |               NA |          NA |                NA |       NA |     0.258 |     0.009 |
| Distance \~ p_h\_water + ec + water_content + fine_silt                         |   2.583 | -52.058 |   5 |  52 |      0.056 | 0.043 |         0.038 |    NA |       NA |            NA |          NA |    NA |     0.033 |               NA |          NA |                NA |       NA |     0.269 |     0.009 |
| Distance \~ ec + water_content + shannon_bak                                    |   2.143 | -52.055 |   4 |  52 |         NA | 0.060 |         0.056 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.052 |                NA |       NA |     0.272 |     0.009 |
| Distance \~ p_h\_water + fine_silt                                              |   1.161 | -52.004 |   3 |  52 |      0.105 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.043 |               NA |          NA |                NA |       NA |     0.322 |     0.009 |
| Distance \~ p_h\_water + ec + fine_silt                                         |   2.470 | -51.996 |   4 |  52 |      0.062 | 0.037 |            NA |    NA |       NA |            NA |          NA |    NA |     0.043 |               NA |          NA |                NA |       NA |     0.330 |     0.009 |
| Distance \~ p_h\_water + ec + water_content + finesilt_and_clay                 |   2.564 | -51.984 |   5 |  52 |      0.056 | 0.043 |         0.041 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.032 |       NA |     0.342 |     0.009 |
| Distance \~ p_h\_water + water_content + oc_beregnet                            |   2.325 | -51.888 |   4 |  52 |      0.108 |    NA |         0.056 |    NA |       NA |            NA |       0.038 |    NA |        NA |               NA |          NA |                NA |       NA |     0.438 |     0.008 |
| Distance \~ p_h\_water + water_content                                          |   1.010 | -51.863 |   3 |  52 |      0.112 |    NA |         0.041 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     0.464 |     0.008 |
| Distance \~ p_h\_water + finesilt_and_clay                                      |   1.140 | -51.784 |   3 |  52 |      0.106 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.039 |       NA |     0.543 |     0.008 |
| Distance \~ p_h\_water + ec + finesilt_and_clay                                 |   2.405 | -51.741 |   4 |  52 |      0.062 | 0.037 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.039 |       NA |     0.585 |     0.008 |
| Distance \~ p_h\_water + water_content + oc_beregnet + fine_silt                |   2.857 | -51.733 |   5 |  52 |      0.093 |    NA |         0.048 |    NA |       NA |            NA |       0.038 |    NA |     0.035 |               NA |          NA |                NA |       NA |     0.594 |     0.008 |
| Distance \~ p_h\_water + ec + water_content + clay                              |   2.462 | -51.717 |   5 |  52 |      0.056 | 0.043 |         0.045 |    NA |       NA |            NA |          NA | 0.028 |        NA |               NA |          NA |                NA |       NA |     0.609 |     0.008 |
| Distance \~ p_h\_water + ec + water_content + wr                                |   4.161 | -51.715 |   5 |  52 |      0.043 | 0.043 |         0.053 | 0.028 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     0.611 |     0.008 |
| Distance \~ ec + water_content + fine_silt + shannon_bak                        |   2.463 | -51.709 |   5 |  52 |         NA | 0.057 |         0.044 |    NA |       NA |            NA |          NA |    NA |     0.032 |               NA |       0.051 |                NA |       NA |     0.617 |     0.008 |
| Distance \~ p_h\_water + water_content + fine_silt                              |   1.706 | -51.702 |   4 |  52 |      0.089 |    NA |         0.033 |    NA |       NA |            NA |          NA |    NA |     0.035 |               NA |          NA |                NA |       NA |     0.624 |     0.008 |
| Distance \~ p_h\_water                                                          |   0.000 | -51.678 |   2 |  52 |      0.115 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     0.648 |     0.007 |
| Distance \~ p_h\_water + ec + water_content + coarse_silt_sand                  |   2.883 | -51.662 |   5 |  52 |      0.055 | 0.039 |         0.039 |    NA |       NA |            NA |          NA |    NA |        NA |            0.027 |          NA |                NA |       NA |     0.665 |     0.007 |
| Distance \~ ec + water_content + shannon_bak + finesilt_and_clay                |   2.420 | -51.649 |   5 |  52 |         NA | 0.058 |         0.047 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.051 |             0.031 |       NA |     0.677 |     0.007 |
| Distance \~ p_h\_water + ec                                                     |   1.841 | -51.634 |   3 |  52 |      0.063 | 0.037 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     0.692 |     0.007 |
| Distance \~ p_h\_water + water_content + finesilt_and_clay                      |   1.542 | -51.615 |   4 |  52 |      0.093 |    NA |         0.035 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.033 |       NA |     0.711 |     0.007 |
| Distance \~ p_h\_water + water_content + oc_beregnet + coarse_silt_sand         |   2.588 | -51.581 |   5 |  52 |      0.096 |    NA |         0.049 |    NA |       NA |            NA |       0.038 |    NA |        NA |            0.033 |          NA |                NA |       NA |     0.745 |     0.007 |
| Distance \~ p_h\_water + water_content + oc_beregnet + finesilt_and_clay        |   2.588 | -51.581 |   5 |  52 |      0.096 |    NA |         0.049 |    NA |       NA |            NA |       0.037 |    NA |        NA |               NA |          NA |             0.033 |       NA |     0.745 |     0.007 |
| Distance \~ p_h\_water + ec + coarse_silt_sand                                  |   2.875 | -51.552 |   4 |  52 |      0.056 | 0.038 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.036 |          NA |                NA |       NA |     0.774 |     0.007 |
| Distance \~ p_h\_water + ec + water_content + ammonium                          |   2.510 | -51.548 |   5 |  52 |      0.056 | 0.045 |         0.044 |    NA |    0.025 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     0.778 |     0.007 |
| Distance \~ p_h\_water + water_content + coarse_silt_sand                       |   1.885 | -51.540 |   4 |  52 |      0.099 |    NA |         0.038 |    NA |       NA |            NA |          NA |    NA |        NA |            0.032 |          NA |                NA |       NA |     0.787 |     0.007 |
| Distance \~ p_h\_water + coarse_silt_sand                                       |   1.030 | -51.528 |   3 |  52 |      0.112 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.035 |          NA |                NA |       NA |     0.798 |     0.007 |
| Distance \~ p_h\_water + ec + ammonium + fine_silt                              |   2.551 | -51.492 |   5 |  52 |      0.059 | 0.038 |            NA |    NA |    0.030 |            NA |          NA |    NA |     0.043 |               NA |          NA |                NA |       NA |     0.834 |     0.007 |
| Distance \~ p_h\_water + ec + water_content + oc_beregnet + fine_silt           |   3.436 | -51.448 |   6 |  52 |      0.049 | 0.033 |         0.042 |    NA |       NA |            NA |       0.028 |    NA |     0.028 |               NA |          NA |                NA |       NA |     0.878 |     0.007 |
| Distance \~ p_h\_water + ammonium + fine_silt                                   |   1.277 | -51.442 |   4 |  52 |      0.098 |    NA |            NA |    NA |    0.028 |            NA |          NA |    NA |     0.043 |               NA |          NA |                NA |       NA |     0.885 |     0.007 |
| Distance \~ ec + water_content + coarse_silt_sand + shannon_bak                 |   2.650 | -51.417 |   5 |  52 |         NA | 0.054 |         0.040 |    NA |       NA |            NA |          NA |    NA |        NA |            0.028 |       0.051 |                NA |       NA |     0.909 |     0.007 |
| Distance \~ ec + water_content + clay + shannon_bak                             |   2.301 | -51.417 |   5 |  52 |         NA | 0.059 |         0.052 |    NA |       NA |            NA |          NA | 0.028 |        NA |               NA |       0.052 |                NA |       NA |     0.910 |     0.007 |
| Distance \~ ec + water_content + wr                                             |   1.310 | -51.372 |   4 |  52 |         NA | 0.073 |         0.058 | 0.041 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     0.954 |     0.006 |
| Distance \~ p_h\_water + water_content + wr                                     |   2.381 | -51.339 |   4 |  52 |      0.072 |    NA |         0.050 | 0.029 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     0.987 |     0.006 |
| Distance \~ p_h\_water + ec + water_content + wr + fine_silt                    |   4.336 | -51.334 |   6 |  52 |      0.041 | 0.039 |         0.045 | 0.027 |       NA |            NA |          NA |    NA |     0.032 |               NA |          NA |                NA |       NA |     0.992 |     0.006 |
| Distance \~ p_h\_water + water_content + clay                                   |   1.252 | -51.328 |   4 |  52 |      0.102 |    NA |         0.038 |    NA |       NA |            NA |          NA | 0.029 |        NA |               NA |          NA |                NA |       NA |     0.998 |     0.006 |
| Distance \~ ec + water_content + oc_beregnet + shannon_bak                      |   2.671 | -51.320 |   5 |  52 |         NA | 0.056 |         0.050 |    NA |       NA |            NA |       0.026 |    NA |        NA |               NA |       0.039 |                NA |       NA |     1.006 |     0.006 |
| Distance \~ p_h\_water + ec + water_content + oc_beregnet + coarse_silt_sand    |   3.292 | -51.316 |   6 |  52 |      0.048 | 0.033 |         0.045 |    NA |       NA |            NA |       0.032 |    NA |        NA |            0.026 |          NA |                NA |       NA |     1.010 |     0.006 |
| Distance \~ p_h\_water + ec + water_content + oc_beregnet + finesilt_and_clay   |   3.292 | -51.316 |   6 |  52 |      0.048 | 0.033 |         0.045 |    NA |       NA |            NA |       0.027 |    NA |        NA |               NA |          NA |             0.026 |       NA |     1.010 |     0.006 |
| Distance \~ p_h\_water + clay                                                   |   1.083 | -51.299 |   3 |  52 |      0.109 |    NA |            NA |    NA |       NA |            NA |          NA | 0.031 |        NA |               NA |          NA |                NA |       NA |     1.028 |     0.006 |
| Distance \~ ec + fine_silt + shannon_bak                                        |   2.463 | -51.281 |   4 |  52 |         NA | 0.056 |            NA |    NA |       NA |            NA |          NA |    NA |     0.044 |               NA |       0.051 |                NA |       NA |     1.045 |     0.006 |
| Distance \~ p_h\_water + ec + ammonium + finesilt_and_clay                      |   2.491 | -51.273 |   5 |  52 |      0.059 | 0.038 |            NA |    NA |    0.030 |            NA |          NA |    NA |        NA |               NA |          NA |             0.040 |       NA |     1.054 |     0.006 |
| Distance \~ p_h\_water + ammonium + finesilt_and_clay                           |   1.249 | -51.271 |   4 |  52 |      0.100 |    NA |            NA |    NA |    0.029 |            NA |          NA |    NA |        NA |               NA |          NA |             0.041 |       NA |     1.055 |     0.006 |
| Distance \~ p_h\_water + water_content + oc_beregnet + clay                     |   2.353 | -51.271 |   5 |  52 |      0.102 |    NA |         0.051 |    NA |       NA |            NA |       0.037 | 0.028 |        NA |               NA |          NA |                NA |       NA |     1.055 |     0.006 |
| Distance \~ ec + coarse_silt_sand + shannon_bak                                 |   2.606 | -51.241 |   4 |  52 |         NA | 0.060 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.043 |       0.051 |                NA |       NA |     1.085 |     0.006 |
| Distance \~ p_h\_water + ec + water_content + wr + finesilt_and_clay            |   4.294 | -51.233 |   6 |  52 |      0.041 | 0.040 |         0.047 | 0.026 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.030 |       NA |     1.093 |     0.006 |
| Distance \~ p_h\_water + ec + water_content + dexter_n                          |   2.280 | -51.222 |   5 |  52 |      0.056 | 0.044 |         0.049 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.020 |     1.105 |     0.006 |
| Distance \~ ec + water_content + oc_beregnet                                    |   1.944 | -51.221 |   4 |  52 |         NA | 0.097 |         0.056 |    NA |       NA |            NA |       0.039 |    NA |        NA |               NA |          NA |                NA |       NA |     1.105 |     0.006 |
| Distance \~ ec + water_content + ammonium + shannon_bak                         |   2.433 | -51.213 |   5 |  52 |         NA | 0.061 |         0.046 |    NA |    0.025 |            NA |          NA |    NA |        NA |               NA |       0.051 |                NA |       NA |     1.114 |     0.006 |
| Distance \~ p_h\_water + water_content + wr + fine_silt                         |   2.595 | -51.198 |   5 |  52 |      0.070 |    NA |         0.044 | 0.030 |       NA |            NA |          NA |    NA |     0.035 |               NA |          NA |                NA |       NA |     1.128 |     0.006 |
| Distance \~ p_h\_water + ec + clay                                              |   2.180 | -51.198 |   4 |  52 |      0.062 | 0.036 |            NA |    NA |       NA |            NA |          NA | 0.031 |        NA |               NA |          NA |                NA |       NA |     1.129 |     0.006 |
| Distance \~ ec + water_content                                                  |   1.067 | -51.167 |   3 |  52 |         NA | 0.101 |         0.055 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.159 |     0.006 |
| Distance \~ p_h\_water + ec + ammonium                                          |   2.092 | -51.146 |   4 |  52 |      0.059 | 0.039 |            NA |    NA |    0.030 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.181 |     0.006 |
| Distance \~ p_h\_water + ec + water_content + oc_beregnet + clay                |   3.182 | -51.140 |   6 |  52 |      0.049 | 0.036 |         0.051 |    NA |       NA |            NA |       0.029 | 0.024 |        NA |               NA |          NA |                NA |       NA |     1.186 |     0.006 |
| Distance \~ ec + water_content + wr + shannon_bak                               |   2.941 | -51.135 |   5 |  52 |         NA | 0.061 |         0.056 | 0.023 |       NA |            NA |          NA |    NA |        NA |               NA |       0.034 |                NA |       NA |     1.191 |     0.006 |
| Distance \~ p_h\_water + clay + fine_silt                                       |   3.129 | -51.107 |   4 |  52 |      0.105 |    NA |            NA |    NA |       NA |            NA |          NA | 0.023 |     0.035 |               NA |          NA |                NA |       NA |     1.219 |     0.006 |
| Distance \~ p_h\_water + ec + oc_beregnet + fine_silt                           |   3.277 | -51.104 |   5 |  52 |      0.048 | 0.039 |            NA |    NA |       NA |            NA |       0.024 |    NA |     0.043 |               NA |          NA |                NA |       NA |     1.222 |     0.006 |
| Distance \~ p_h\_water + ec + ammonium + coarse_silt_sand                       |   2.880 | -51.094 |   5 |  52 |      0.055 | 0.038 |            NA |    NA |    0.031 |            NA |          NA |    NA |        NA |            0.037 |          NA |                NA |       NA |     1.232 |     0.006 |
| Distance \~ p_h\_water + ammonium                                               |   1.046 | -51.089 |   3 |  52 |      0.111 |    NA |            NA |    NA |    0.028 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.237 |     0.006 |
| Distance \~ ec + water_content + wr + fine_silt                                 |   2.237 | -51.088 |   5 |  52 |         NA | 0.068 |         0.046 | 0.041 |       NA |            NA |          NA |    NA |     0.033 |               NA |          NA |                NA |       NA |     1.239 |     0.006 |
| Distance \~ p_h\_water + water_content + ammonium                               |   2.295 | -51.077 |   4 |  52 |      0.112 |    NA |         0.038 |    NA |    0.025 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.249 |     0.006 |
| Distance \~ p_h\_water + ammonium + coarse_silt_sand                            |   1.287 | -51.071 |   4 |  52 |      0.104 |    NA |            NA |    NA |    0.030 |            NA |          NA |    NA |        NA |            0.038 |          NA |                NA |       NA |     1.255 |     0.006 |
| Distance \~ p_h\_water + water_content + wr + finesilt_and_clay                 |   2.515 | -51.067 |   5 |  52 |      0.070 |    NA |         0.045 | 0.029 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.034 |       NA |     1.259 |     0.006 |
| Distance \~ p_h\_water + ec + clay + fine_silt                                  |   3.545 | -51.064 |   5 |  52 |      0.062 | 0.037 |            NA |    NA |       NA |            NA |          NA | 0.023 |     0.036 |               NA |          NA |                NA |       NA |     1.262 |     0.005 |
| Distance \~ p_h\_water + ec + water_content + wr + oc_beregnet                  |   4.304 | -51.056 |   6 |  52 |      0.044 | 0.040 |         0.052 | 0.023 |       NA |            NA |       0.028 |    NA |        NA |               NA |          NA |                NA |       NA |     1.270 |     0.005 |
| Distance \~ p_h\_water + oc_beregnet + fine_silt                                |   1.491 | -51.049 |   4 |  52 |      0.093 |    NA |            NA |    NA |       NA |            NA |       0.022 |    NA |     0.043 |               NA |          NA |                NA |       NA |     1.277 |     0.005 |
| Distance \~ ec + shannon_bak + finesilt_and_clay                                |   2.414 | -51.045 |   4 |  52 |         NA | 0.055 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.051 |             0.040 |       NA |     1.281 |     0.005 |
| Distance \~ p_h\_water + ec + water_content + nitrat_nitrit                     |   2.452 | -51.008 |   5 |  52 |      0.055 | 0.043 |         0.042 |    NA |       NA |         0.017 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.318 |     0.005 |
| Distance \~ p_h\_water + water_content + wr + coarse_silt_sand                  |   2.393 | -51.006 |   5 |  52 |      0.069 |    NA |         0.045 | 0.029 |       NA |            NA |          NA |    NA |        NA |            0.033 |          NA |                NA |       NA |     1.320 |     0.005 |
| Distance \~ ec + water_content + wr + finesilt_and_clay                         |   1.912 | -50.993 |   5 |  52 |         NA | 0.069 |         0.049 | 0.041 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.032 |       NA |     1.333 |     0.005 |
| Distance \~ p_h\_water + ec + water_content + ammonium + oc_beregnet            |   3.146 | -50.984 |   6 |  52 |      0.050 | 0.039 |         0.054 |    NA |    0.021 |            NA |       0.029 |    NA |        NA |               NA |          NA |                NA |       NA |     1.343 |     0.005 |
| Distance \~ p_h\_water + ec + water_content + wr + clay                         |   4.223 | -50.979 |   6 |  52 |      0.042 | 0.041 |         0.051 | 0.027 |       NA |            NA |          NA | 0.027 |        NA |               NA |          NA |                NA |       NA |     1.347 |     0.005 |
| Distance \~ p_h\_water + clay + coarse_silt_sand                                |   4.325 | -50.979 |   4 |  52 |      0.108 |    NA |            NA |    NA |       NA |            NA |          NA | 0.029 |        NA |            0.033 |          NA |                NA |       NA |     1.347 |     0.005 |
| Distance \~ p_h\_water + ec + water_content + wr + coarse_silt_sand             |   4.403 | -50.975 |   6 |  52 |      0.041 | 0.037 |         0.045 | 0.027 |       NA |            NA |          NA |    NA |        NA |            0.027 |          NA |                NA |       NA |     1.351 |     0.005 |
| Distance \~ p_h\_water + ec + water_content + clay + fine_silt                  |   4.281 | -50.971 |   6 |  52 |      0.055 | 0.043 |         0.036 |    NA |       NA |            NA |          NA | 0.021 |     0.027 |               NA |          NA |                NA |       NA |     1.355 |     0.005 |
| Distance \~ ec + ammonium + fine_silt + shannon_bak                             |   2.469 | -50.957 |   5 |  52 |         NA | 0.057 |            NA |    NA |    0.033 |            NA |          NA |    NA |     0.042 |               NA |       0.051 |                NA |       NA |     1.370 |     0.005 |
| Distance \~ ec + shannon_bak                                                    |   2.052 | -50.932 |   3 |  52 |         NA | 0.052 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.052 |                NA |       NA |     1.394 |     0.005 |
| Distance \~ ec + water_content + shannon_bak + dexter_n                         |   2.175 | -50.921 |   5 |  52 |         NA | 0.060 |         0.056 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.052 |                NA |    0.020 |     1.405 |     0.005 |
| Distance \~ p_h\_water + nitrat_nitrit                                          |   1.157 | -50.899 |   3 |  52 |      0.102 |    NA |            NA |    NA |       NA |         0.025 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.427 |     0.005 |
| Distance \~ p_h\_water + water_content + wr + oc_beregnet                       |   3.225 | -50.888 |   5 |  52 |      0.072 |    NA |         0.050 | 0.022 |       NA |            NA |       0.031 |    NA |        NA |               NA |          NA |                NA |       NA |     1.439 |     0.005 |
| Distance \~ ec + clay + coarse_silt_sand + shannon_bak                          |   4.817 | -50.877 |   5 |  52 |         NA | 0.063 |            NA |    NA |       NA |            NA |          NA | 0.032 |        NA |            0.044 |       0.050 |                NA |       NA |     1.449 |     0.005 |
| Distance \~ p_h\_water + water_content + ammonium + oc_beregnet                 |   2.874 | -50.869 |   5 |  52 |      0.109 |    NA |         0.053 |    NA |    0.022 |            NA |       0.035 |    NA |        NA |               NA |          NA |                NA |       NA |     1.458 |     0.005 |
| Distance \~ p_h\_water + oc_beregnet + coarse_silt_sand                         |   1.685 | -50.867 |   4 |  52 |      0.094 |    NA |            NA |    NA |       NA |            NA |       0.027 |    NA |        NA |            0.040 |          NA |                NA |       NA |     1.460 |     0.005 |
| Distance \~ p_h\_water + oc_beregnet + finesilt_and_clay                        |   1.453 | -50.867 |   4 |  52 |      0.094 |    NA |            NA |    NA |       NA |            NA |       0.023 |    NA |        NA |               NA |          NA |             0.040 |       NA |     1.460 |     0.005 |
| Distance \~ p_h\_water + ec + water_content + ammonium + finesilt_and_clay      |   2.941 | -50.857 |   6 |  52 |      0.057 | 0.043 |         0.032 |    NA |    0.021 |            NA |          NA |    NA |        NA |               NA |          NA |             0.027 |       NA |     1.469 |     0.005 |
| Distance \~ ec + water_content + fine_silt                                      |   1.908 | -50.855 |   4 |  52 |         NA | 0.076 |         0.044 |    NA |       NA |            NA |          NA |    NA |     0.033 |               NA |          NA |                NA |       NA |     1.471 |     0.005 |
| Distance \~ ec + ammonium + coarse_silt_sand + shannon_bak                      |   2.619 | -50.855 |   5 |  52 |         NA | 0.059 |            NA |    NA |    0.032 |            NA |          NA |    NA |        NA |            0.040 |       0.052 |                NA |       NA |     1.471 |     0.005 |
| Distance \~ p_h\_water + ec + oc_beregnet + coarse_silt_sand                    |   3.215 | -50.830 |   5 |  52 |      0.048 | 0.037 |            NA |    NA |       NA |            NA |       0.027 |    NA |        NA |            0.039 |          NA |                NA |       NA |     1.497 |     0.005 |
| Distance \~ p_h\_water + ec + oc_beregnet + finesilt_and_clay                   |   3.215 | -50.830 |   5 |  52 |      0.048 | 0.037 |            NA |    NA |       NA |            NA |       0.024 |    NA |        NA |               NA |          NA |             0.039 |       NA |     1.497 |     0.005 |
| Distance \~ p_h\_water + ec + wr + fine_silt                                    |   4.327 | -50.827 |   5 |  52 |      0.042 | 0.038 |            NA | 0.020 |       NA |            NA |          NA |    NA |     0.039 |               NA |          NA |                NA |       NA |     1.499 |     0.005 |
| Distance \~ ec + oc_beregnet + fine_silt + shannon_bak                          |   2.695 | -50.825 |   5 |  52 |         NA | 0.063 |            NA |    NA |       NA |            NA |       0.031 |    NA |     0.043 |               NA |       0.044 |                NA |       NA |     1.501 |     0.005 |
| Distance \~ ec + water_content + wr + fine_silt + shannon_bak                   |   2.968 | -50.816 |   6 |  52 |         NA | 0.058 |         0.044 | 0.024 |       NA |            NA |          NA |    NA |     0.033 |               NA |       0.034 |                NA |       NA |     1.510 |     0.005 |
| Distance \~ p_h\_water + ec + water_content + ammonium + clay                   |   2.618 | -50.810 |   6 |  52 |      0.056 | 0.044 |         0.039 |    NA |    0.024 |            NA |          NA | 0.027 |        NA |               NA |          NA |                NA |       NA |     1.517 |     0.005 |
| Distance \~ p_h\_water + ec + water_content + ammonium + fine_silt              |   3.233 | -50.808 |   6 |  52 |      0.056 | 0.042 |         0.028 |    NA |    0.019 |            NA |          NA |    NA |     0.027 |               NA |          NA |                NA |       NA |     1.519 |     0.005 |
| Distance \~ p_h\_water + ammonium + clay                                        |   1.166 | -50.806 |   4 |  52 |      0.104 |    NA |            NA |    NA |    0.030 |            NA |          NA | 0.033 |        NA |               NA |          NA |                NA |       NA |     1.521 |     0.005 |
| Distance \~ p_h\_water + wr + fine_silt                                         |   2.324 | -50.797 |   4 |  52 |      0.070 |    NA |            NA | 0.018 |       NA |            NA |          NA |    NA |     0.041 |               NA |          NA |                NA |       NA |     1.529 |     0.005 |
| Distance \~ p_h\_water + water_content + dexter_n                               |   1.041 | -50.790 |   4 |  52 |      0.112 |    NA |         0.042 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.021 |     1.536 |     0.005 |
| Distance \~ ec + water_content + finesilt_and_clay                              |   1.726 | -50.777 |   4 |  52 |         NA | 0.080 |         0.047 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.032 |       NA |     1.550 |     0.005 |
| Distance \~ ec + ammonium + shannon_bak + finesilt_and_clay                     |   2.422 | -50.765 |   5 |  52 |         NA | 0.057 |            NA |    NA |    0.034 |            NA |          NA |    NA |        NA |               NA |       0.051 |             0.039 |       NA |     1.561 |     0.005 |
| Distance \~ p_h\_water + ec + oc_beregnet                                       |   3.120 | -50.763 |   4 |  52 |      0.049 | 0.039 |            NA |    NA |       NA |            NA |       0.024 |    NA |        NA |               NA |          NA |                NA |       NA |     1.563 |     0.005 |
| Distance \~ ec + water_content + wr + coarse_silt_sand                          |   2.367 | -50.757 |   5 |  52 |         NA | 0.065 |         0.045 | 0.041 |       NA |            NA |          NA |    NA |        NA |            0.028 |          NA |                NA |       NA |     1.570 |     0.005 |
| Distance \~ p_h\_water + water_content + wr + clay                              |   2.420 | -50.748 |   5 |  52 |      0.071 |    NA |         0.048 | 0.029 |       NA |            NA |          NA | 0.029 |        NA |               NA |          NA |                NA |       NA |     1.578 |     0.005 |
| Distance \~ p_h\_water + ec + wr                                                |   4.030 | -50.747 |   4 |  52 |      0.048 | 0.041 |            NA | 0.024 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.579 |     0.005 |
| Distance \~ p_h\_water + ec + ammonium + clay                                   |   2.292 | -50.745 |   5 |  52 |      0.059 | 0.037 |            NA |    NA |    0.031 |            NA |          NA | 0.032 |        NA |               NA |          NA |                NA |       NA |     1.581 |     0.005 |
| Distance \~ p_h\_water + water_content + wr + oc_beregnet + fine_silt           |   3.349 | -50.744 |   6 |  52 |      0.068 |    NA |         0.047 | 0.023 |       NA |            NA |       0.031 |    NA |     0.036 |               NA |          NA |                NA |       NA |     1.582 |     0.005 |
| Distance \~ p_h\_water + ec + water_content + wr + ammonium                     |   4.187 | -50.742 |   6 |  52 |      0.042 | 0.043 |         0.042 | 0.026 |    0.023 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.584 |     0.005 |
| Distance \~ p_h\_water + water_content + oc_beregnet + clay + fine_silt         |   4.362 | -50.740 |   6 |  52 |      0.092 |    NA |         0.048 |    NA |       NA |            NA |       0.039 | 0.023 |     0.030 |               NA |          NA |                NA |       NA |     1.587 |     0.005 |
| Distance \~ p_h\_water + oc_beregnet                                            |   1.211 | -50.734 |   3 |  52 |      0.103 |    NA |            NA |    NA |       NA |            NA |       0.022 |    NA |        NA |               NA |          NA |                NA |       NA |     1.592 |     0.005 |
| Distance \~ p_h\_water + fine_silt + dexter_n                                   |   1.219 | -50.729 |   4 |  52 |      0.105 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.041 |               NA |          NA |                NA |    0.017 |     1.597 |     0.005 |
| Distance \~ ec + water_content + oc_beregnet + fine_silt + shannon_bak          |   3.225 | -50.729 |   6 |  52 |         NA | 0.049 |         0.036 |    NA |       NA |            NA |       0.023 |    NA |     0.029 |               NA |       0.039 |                NA |       NA |     1.597 |     0.005 |
| Distance \~ ec + water_content + nitrat_nitrit + shannon_bak                    |   2.343 | -50.728 |   5 |  52 |         NA | 0.056 |         0.050 |    NA |       NA |         0.017 |          NA |    NA |        NA |               NA |       0.051 |                NA |       NA |     1.599 |     0.005 |
| Distance \~ ec + water_content + wr + shannon_bak + finesilt_and_clay           |   2.953 | -50.724 |   6 |  52 |         NA | 0.058 |         0.047 | 0.024 |       NA |            NA |          NA |    NA |        NA |               NA |       0.034 |             0.032 |       NA |     1.602 |     0.005 |
| Distance \~ ec + water_content + wr + clay                                      |   1.469 | -50.721 |   5 |  52 |         NA | 0.071 |         0.054 | 0.041 |       NA |            NA |          NA | 0.028 |        NA |               NA |          NA |                NA |       NA |     1.605 |     0.005 |
| Distance \~ p_h\_water + ec + nitrat_nitrit                                     |   2.212 | -50.714 |   4 |  52 |      0.063 | 0.035 |            NA |    NA |       NA |         0.023 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.612 |     0.005 |
| Distance \~ ec + ammonium + shannon_bak                                         |   2.092 | -50.695 |   4 |  52 |         NA | 0.055 |            NA |    NA |    0.034 |            NA |          NA |    NA |        NA |               NA |       0.052 |                NA |       NA |     1.631 |     0.005 |
| Distance \~ p_h\_water + water_content + nitrat_nitrit                          |   1.424 | -50.691 |   4 |  52 |      0.090 |    NA |         0.035 |    NA |       NA |         0.019 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.635 |     0.005 |
| Distance \~ p_h\_water + water_content + nitrat_nitrit + oc_beregnet            |   2.339 | -50.689 |   5 |  52 |      0.091 |    NA |         0.052 |    NA |       NA |         0.019 |       0.038 |    NA |        NA |               NA |          NA |                NA |       NA |     1.637 |     0.005 |
| Distance \~ ec + water_content + oc_beregnet + fine_silt                        |   2.921 | -50.685 |   5 |  52 |         NA | 0.077 |         0.042 |    NA |       NA |            NA |       0.035 |    NA |     0.030 |               NA |          NA |                NA |       NA |     1.641 |     0.005 |
| Distance \~ p_h\_water + ec + water_content + oc_beregnet + dexter_n            |   3.146 | -50.672 |   6 |  52 |      0.050 | 0.040 |         0.057 |    NA |       NA |            NA |       0.030 |    NA |        NA |               NA |          NA |                NA |    0.017 |     1.654 |     0.005 |
| Distance \~ p_h\_water + ec + water_content + ammonium + coarse_silt_sand       |   3.356 | -50.657 |   6 |  52 |      0.055 | 0.039 |         0.031 |    NA |    0.023 |            NA |          NA |    NA |        NA |            0.025 |          NA |                NA |       NA |     1.669 |     0.004 |
| Distance \~ p_h\_water + ec + water_content + nitrat_nitrit + oc_beregnet       |   3.146 | -50.651 |   6 |  52 |      0.050 | 0.037 |         0.051 |    NA |       NA |         0.017 |       0.033 |    NA |        NA |               NA |          NA |                NA |       NA |     1.675 |     0.004 |
| Distance \~ water_content + shannon_bak                                         |   1.023 | -50.647 |   3 |  52 |         NA |    NA |         0.047 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.092 |                NA |       NA |     1.679 |     0.004 |
| Distance \~ p_h\_water + ec + fine_silt + dexter_n                              |   2.668 | -50.644 |   5 |  52 |      0.061 | 0.037 |            NA |    NA |       NA |            NA |          NA |    NA |     0.040 |               NA |          NA |                NA |    0.017 |     1.683 |     0.004 |
| Distance \~ p_h\_water + water_content + clay + fine_silt                       |   4.079 | -50.642 |   5 |  52 |      0.088 |    NA |         0.031 |    NA |       NA |            NA |          NA | 0.021 |     0.027 |               NA |          NA |                NA |       NA |     1.684 |     0.004 |
| Distance \~ ec + water_content + clay + fine_silt + shannon_bak                 |   4.230 | -50.642 |   6 |  52 |         NA | 0.057 |         0.043 |    NA |       NA |            NA |          NA | 0.022 |     0.026 |               NA |       0.051 |                NA |       NA |     1.684 |     0.004 |
| Distance \~ ec + water_content + oc_beregnet + coarse_silt_sand + shannon_bak   |   2.860 | -50.627 |   6 |  52 |         NA | 0.050 |         0.038 |    NA |       NA |            NA |       0.026 |    NA |        NA |            0.027 |       0.039 |                NA |       NA |     1.699 |     0.004 |
| Distance \~ ec + water_content + oc_beregnet + shannon_bak + finesilt_and_clay  |   2.860 | -50.627 |   6 |  52 |         NA | 0.050 |         0.038 |    NA |       NA |            NA |       0.023 |    NA |        NA |               NA |       0.039 |             0.027 |       NA |     1.699 |     0.004 |
| Distance \~ ec + water_content + wr + oc_beregnet                               |   3.326 | -50.618 |   5 |  52 |         NA | 0.068 |         0.052 | 0.029 |       NA |            NA |       0.026 |    NA |        NA |               NA |          NA |                NA |       NA |     1.708 |     0.004 |
| Distance \~ p_h\_water + nitrat_nitrit + fine_silt                              |   1.645 | -50.599 |   4 |  52 |      0.101 |    NA |            NA |    NA |       NA |         0.015 |          NA |    NA |     0.033 |               NA |          NA |                NA |       NA |     1.727 |     0.004 |
| Distance \~ ec + oc_beregnet + coarse_silt_sand + shannon_bak                   |   2.649 | -50.597 |   5 |  52 |         NA | 0.062 |            NA |    NA |       NA |            NA |       0.028 |    NA |        NA |            0.039 |       0.044 |                NA |       NA |     1.729 |     0.004 |
| Distance \~ ec + oc_beregnet + shannon_bak + finesilt_and_clay                  |   2.649 | -50.597 |   5 |  52 |         NA | 0.062 |            NA |    NA |       NA |            NA |       0.031 |    NA |        NA |               NA |       0.044 |             0.039 |       NA |     1.729 |     0.004 |
| Distance \~ p_h\_water + wr                                                     |   2.230 | -50.596 |   3 |  52 |      0.068 |    NA |            NA | 0.020 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.731 |     0.004 |
| Distance \~ p_h\_water + wr + finesilt_and_clay                                 |   2.323 | -50.595 |   4 |  52 |      0.070 |    NA |            NA | 0.019 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.038 |       NA |     1.731 |     0.004 |
| Distance \~ p_h\_water + ec + wr + finesilt_and_clay                            |   4.266 | -50.589 |   5 |  52 |      0.043 | 0.038 |            NA | 0.020 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.036 |       NA |     1.737 |     0.004 |
| Distance \~ ec + water_content + oc_beregnet + coarse_silt_sand                 |   2.593 | -50.583 |   5 |  52 |         NA | 0.081 |         0.044 |    NA |       NA |            NA |       0.039 |    NA |        NA |            0.028 |          NA |                NA |       NA |     1.743 |     0.004 |
| Distance \~ ec + water_content + oc_beregnet + finesilt_and_clay                |   2.593 | -50.583 |   5 |  52 |         NA | 0.081 |         0.044 |    NA |       NA |            NA |       0.035 |    NA |        NA |               NA |          NA |             0.028 |       NA |     1.743 |     0.004 |
| Distance \~ p_h\_water + dexter_n                                               |   1.023 | -50.581 |   3 |  52 |      0.114 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.019 |     1.745 |     0.004 |
| Distance \~ p_h\_water + water_content + ammonium + coarse_silt_sand            |   3.352 | -50.578 |   5 |  52 |      0.101 |    NA |         0.030 |    NA |    0.023 |            NA |          NA |    NA |        NA |            0.030 |          NA |                NA |       NA |     1.748 |     0.004 |
| Distance \~ p_h\_water + ec + wr + ammonium + fine_silt                         |   4.327 | -50.577 |   6 |  52 |      0.042 | 0.038 |            NA | 0.024 |    0.034 |            NA |          NA |    NA |     0.039 |               NA |          NA |                NA |       NA |     1.749 |     0.004 |
| Distance \~ ec + wr + fine_silt                                                 |   1.615 | -50.577 |   4 |  52 |         NA | 0.066 |            NA | 0.040 |       NA |            NA |          NA |    NA |     0.045 |               NA |          NA |                NA |       NA |     1.750 |     0.004 |
| Distance \~ p_h\_water + water_content + wr + oc_beregnet + coarse_silt_sand    |   3.302 | -50.575 |   6 |  52 |      0.069 |    NA |         0.046 | 0.023 |       NA |            NA |       0.031 |    NA |        NA |            0.033 |          NA |                NA |       NA |     1.751 |     0.004 |
| Distance \~ p_h\_water + water_content + wr + oc_beregnet + finesilt_and_clay   |   3.302 | -50.575 |   6 |  52 |      0.069 |    NA |         0.046 | 0.023 |       NA |            NA |       0.031 |    NA |        NA |               NA |          NA |             0.033 |       NA |     1.751 |     0.004 |
| Distance \~ p_h\_water + wr + coarse_silt_sand                                  |   2.326 | -50.556 |   4 |  52 |      0.070 |    NA |            NA | 0.022 |       NA |            NA |          NA |    NA |        NA |            0.038 |          NA |                NA |       NA |     1.770 |     0.004 |
| Distance \~ p_h\_water + ec + nitrat_nitrit + fine_silt                         |   2.560 | -50.555 |   5 |  52 |      0.062 | 0.037 |            NA |    NA |       NA |         0.016 |          NA |    NA |     0.036 |               NA |          NA |                NA |       NA |     1.771 |     0.004 |
| Distance \~ p_h\_water + water_content + oc_beregnet + dexter_n                 |   2.476 | -50.554 |   5 |  52 |      0.108 |    NA |         0.055 |    NA |       NA |            NA |       0.034 |    NA |        NA |               NA |          NA |                NA |    0.017 |     1.772 |     0.004 |
| Distance \~ p_h\_water + wr + ammonium + fine_silt                              |   2.573 | -50.546 |   5 |  52 |      0.070 |    NA |            NA | 0.024 |    0.034 |            NA |          NA |    NA |     0.041 |               NA |          NA |                NA |       NA |     1.780 |     0.004 |
| Distance \~ ec + water_content + coarse_silt_sand                               |   2.358 | -50.539 |   4 |  52 |         NA | 0.083 |         0.040 |    NA |       NA |            NA |          NA |    NA |        NA |            0.028 |          NA |                NA |       NA |     1.787 |     0.004 |
| Distance \~ p_h\_water + nitrat_nitrit + finesilt_and_clay                      |   1.461 | -50.538 |   4 |  52 |      0.101 |    NA |            NA |    NA |       NA |         0.018 |          NA |    NA |        NA |               NA |          NA |             0.032 |       NA |     1.788 |     0.004 |
| Distance \~ p_h\_water + water_content + ammonium + fine_silt                   |   3.188 | -50.532 |   5 |  52 |      0.091 |    NA |         0.024 |    NA |    0.020 |            NA |          NA |    NA |     0.030 |               NA |          NA |                NA |       NA |     1.794 |     0.004 |
| Distance \~ ec + water_content + clay                                           |   1.357 | -50.530 |   4 |  52 |         NA | 0.089 |         0.052 |    NA |       NA |            NA |          NA | 0.028 |        NA |               NA |          NA |                NA |       NA |     1.796 |     0.004 |
| Distance \~ p_h\_water + water_content + ammonium + finesilt_and_clay           |   2.857 | -50.530 |   5 |  52 |      0.095 |    NA |         0.026 |    NA |    0.021 |            NA |          NA |    NA |        NA |               NA |          NA |             0.030 |       NA |     1.796 |     0.004 |
| Distance \~ p_h\_water + ec + ammonium + clay + fine_silt                       |   3.578 | -50.529 |   6 |  52 |      0.059 | 0.038 |            NA |    NA |    0.030 |            NA |          NA | 0.023 |     0.035 |               NA |          NA |                NA |       NA |     1.798 |     0.004 |
| Distance \~ ec + oc_beregnet + shannon_bak                                      |   2.406 | -50.528 |   4 |  52 |         NA | 0.059 |            NA |    NA |       NA |            NA |       0.032 |    NA |        NA |               NA |       0.045 |                NA |       NA |     1.798 |     0.004 |
| Distance \~ p_h\_water + ec + dexter_n                                          |   1.910 | -50.524 |   4 |  52 |      0.063 | 0.037 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.020 |     1.802 |     0.004 |
| Distance \~ p_h\_water + ec + wr + coarse_silt_sand                             |   4.403 | -50.520 |   5 |  52 |      0.041 | 0.038 |            NA | 0.022 |       NA |            NA |          NA |    NA |        NA |            0.035 |          NA |                NA |       NA |     1.806 |     0.004 |
| Distance \~ p_h\_water + ec + water_content + nitrat_nitrit + finesilt_and_clay |   2.660 | -50.513 |   6 |  52 |      0.055 | 0.042 |         0.039 |    NA |       NA |         0.016 |          NA |    NA |        NA |               NA |          NA |             0.031 |       NA |     1.813 |     0.004 |
| Distance \~ ec + clay + shannon_bak                                             |   2.263 | -50.512 |   4 |  52 |         NA | 0.053 |            NA |    NA |       NA |            NA |          NA | 0.031 |        NA |               NA |       0.052 |                NA |       NA |     1.814 |     0.004 |
| Distance \~ p_h\_water + ammonium + clay + fine_silt                            |   3.215 | -50.509 |   5 |  52 |      0.098 |    NA |            NA |    NA |    0.029 |            NA |          NA | 0.023 |     0.033 |               NA |          NA |                NA |       NA |     1.818 |     0.004 |
| Distance \~ p_h\_water + ec + water_content + nitrat_nitrit + fine_silt         |   2.647 | -50.498 |   6 |  52 |      0.055 | 0.042 |         0.037 |    NA |       NA |         0.015 |          NA |    NA |     0.030 |               NA |          NA |                NA |       NA |     1.828 |     0.004 |
| Distance \~ ec + water_content + wr + ammonium                                  |   2.437 | -50.488 |   5 |  52 |         NA | 0.074 |         0.047 | 0.040 |    0.024 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.838 |     0.004 |
| Distance \~ ec + water_content + ammonium + shannon_bak + finesilt_and_clay     |   2.818 | -50.481 |   6 |  52 |         NA | 0.057 |         0.034 |    NA |    0.020 |            NA |          NA |    NA |        NA |               NA |       0.051 |             0.027 |       NA |     1.845 |     0.004 |
| Distance \~ p_h\_water + ec + wr + ammonium                                     |   4.050 | -50.481 |   5 |  52 |      0.047 | 0.040 |            NA | 0.028 |    0.034 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.846 |     0.004 |
| Distance \~ ec + fine_silt                                                      |   1.554 | -50.476 |   3 |  52 |         NA | 0.080 |            NA |    NA |       NA |            NA |          NA |    NA |     0.044 |               NA |          NA |                NA |       NA |     1.851 |     0.004 |
| Distance \~ p_h\_water + ec + nitrat_nitrit + finesilt_and_clay                 |   2.556 | -50.468 |   5 |  52 |      0.062 | 0.037 |            NA |    NA |       NA |         0.018 |          NA |    NA |        NA |               NA |          NA |             0.034 |       NA |     1.859 |     0.004 |
| Distance \~ water_content + fine_silt + shannon_bak                             |   1.691 | -50.467 |   4 |  52 |         NA |    NA |         0.043 |    NA |       NA |            NA |          NA |    NA |     0.035 |               NA |       0.070 |                NA |       NA |     1.859 |     0.004 |
| Distance \~ p_h\_water + ec + water_content + fine_silt + dexter_n              |   2.717 | -50.462 |   6 |  52 |      0.056 | 0.041 |         0.035 |    NA |       NA |            NA |          NA |    NA |     0.027 |               NA |          NA |                NA |    0.014 |     1.864 |     0.004 |
| Distance \~ ec + water_content + oc_beregnet + clay + shannon_bak               |   2.673 | -50.459 |   6 |  52 |         NA | 0.053 |         0.044 |    NA |       NA |            NA |       0.024 | 0.025 |        NA |               NA |       0.039 |                NA |       NA |     1.867 |     0.004 |
| Distance \~ ec + water_content + ammonium + clay + shannon_bak                  |   2.530 | -50.457 |   6 |  52 |         NA | 0.059 |         0.041 |    NA |    0.024 |            NA |          NA | 0.027 |        NA |               NA |       0.051 |                NA |       NA |     1.869 |     0.004 |
| Distance \~ p_h\_water + finesilt_and_clay + dexter_n                           |   1.261 | -50.455 |   4 |  52 |      0.106 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.036 |    0.016 |     1.871 |     0.004 |
| Distance \~ ec + water_content + wr + coarse_silt_sand + shannon_bak            |   2.985 | -50.451 |   6 |  52 |         NA | 0.055 |         0.042 | 0.023 |       NA |            NA |          NA |    NA |        NA |            0.028 |       0.033 |                NA |       NA |     1.875 |     0.004 |
| Distance \~ ec + water_content + wr + clay + shannon_bak                        |   2.942 | -50.451 |   6 |  52 |         NA | 0.060 |         0.052 | 0.023 |       NA |            NA |          NA | 0.028 |        NA |               NA |       0.034 |                NA |       NA |     1.875 |     0.004 |
| Distance \~ ec + oc_beregnet + fine_silt                                        |   1.560 | -50.442 |   4 |  52 |         NA | 0.083 |            NA |    NA |       NA |            NA |       0.038 |    NA |     0.044 |               NA |          NA |                NA |       NA |     1.884 |     0.004 |
| Distance \~ ec + wr + fine_silt + shannon_bak                                   |   2.839 | -50.429 |   5 |  52 |         NA | 0.061 |            NA | 0.025 |       NA |            NA |          NA |    NA |     0.045 |               NA |       0.036 |                NA |       NA |     1.897 |     0.004 |
| Distance \~ ec + water_content + ammonium + fine_silt + shannon_bak             |   3.077 | -50.426 |   6 |  52 |         NA | 0.056 |         0.030 |    NA |    0.019 |            NA |          NA |    NA |     0.026 |               NA |       0.051 |                NA |       NA |     1.900 |     0.004 |
| Distance \~ p_h\_water + ec + water_content + wr + dexter_n                     |   4.166 | -50.424 |   6 |  52 |      0.042 | 0.043 |         0.052 | 0.026 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.019 |     1.902 |     0.004 |
| Distance \~ ec + coarse_silt_sand                                               |   1.437 | -50.417 |   3 |  52 |         NA | 0.093 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.043 |          NA |                NA |       NA |     1.909 |     0.004 |
| Distance \~ p_h\_water + water_content + ammonium + clay                        |   2.434 | -50.411 |   5 |  52 |      0.103 |    NA |         0.032 |    NA |    0.024 |            NA |          NA | 0.028 |        NA |               NA |          NA |                NA |       NA |     1.915 |     0.004 |
| Distance \~ p_h\_water + water_content + ammonium + oc_beregnet + fine_silt     |   4.299 | -50.407 |   6 |  52 |      0.093 |    NA |         0.041 |    NA |    0.018 |            NA |       0.036 |    NA |     0.031 |               NA |          NA |                NA |       NA |     1.919 |     0.004 |
| Distance \~ ec + water_content + oc_beregnet + clay                             |   2.166 | -50.403 |   5 |  52 |         NA | 0.089 |         0.050 |    NA |       NA |            NA |       0.036 | 0.025 |        NA |               NA |          NA |                NA |       NA |     1.923 |     0.004 |
| Distance \~ p_h\_water + wr + ammonium + finesilt_and_clay                      |   2.546 | -50.398 |   5 |  52 |      0.070 |    NA |            NA | 0.024 |    0.035 |            NA |          NA |    NA |        NA |               NA |          NA |             0.039 |       NA |     1.928 |     0.004 |
| Distance \~ p_h\_water + oc_beregnet + clay                                     |   1.353 | -50.396 |   4 |  52 |      0.098 |    NA |            NA |    NA |       NA |            NA |       0.023 | 0.033 |        NA |               NA |          NA |                NA |       NA |     1.930 |     0.004 |
| Distance \~ p_h\_water + ec + wr + ammonium + finesilt_and_clay                 |   4.267 | -50.387 |   6 |  52 |      0.043 | 0.038 |            NA | 0.025 |    0.035 |            NA |          NA |    NA |        NA |               NA |          NA |             0.037 |       NA |     1.939 |     0.004 |
| Distance \~ p_h\_water + water_content + wr + ammonium                          |   2.480 | -50.372 |   5 |  52 |      0.072 |    NA |         0.039 | 0.027 |    0.023 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.954 |     0.004 |
| Distance \~ ec + water_content + ammonium                                       |   2.417 | -50.358 |   4 |  52 |         NA | 0.100 |         0.046 |    NA |    0.025 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.968 |     0.004 |
| Distance \~ water_content + shannon_bak + finesilt_and_clay                     |   1.550 | -50.357 |   4 |  52 |         NA |    NA |         0.044 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.073 |             0.034 |       NA |     1.969 |     0.004 |
| Distance \~ water_content + coarse_silt_sand + shannon_bak                      |   1.916 | -50.353 |   4 |  52 |         NA |    NA |         0.046 |    NA |       NA |            NA |          NA |    NA |        NA |            0.034 |       0.080 |                NA |       NA |     1.974 |     0.004 |
| Distance \~ ec + water_content + ammonium + coarse_silt_sand + shannon_bak      |   3.368 | -50.348 |   6 |  52 |         NA | 0.054 |         0.030 |    NA |    0.022 |            NA |          NA |    NA |        NA |            0.025 |       0.051 |                NA |       NA |     1.978 |     0.004 |
| Distance \~ ec + wr + coarse_silt_sand                                          |   1.682 | -50.346 |   4 |  52 |         NA | 0.066 |            NA | 0.037 |       NA |            NA |          NA |    NA |        NA |            0.042 |          NA |                NA |       NA |     1.981 |     0.004 |
| Distance \~ p_h\_water + wr + ammonium                                          |   2.477 | -50.344 |   4 |  52 |      0.070 |    NA |            NA | 0.026 |    0.034 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.983 |     0.004 |
| Distance \~ p_h\_water + ec + finesilt_and_clay + dexter_n                      |   2.689 | -50.332 |   5 |  52 |      0.060 | 0.036 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.035 |    0.016 |     1.994 |     0.004 |

Table 2.2: Best models for vegetation presence absence

</details>

## 2.3 Abundance data

Since the environmental variables are the same, we don´t need to do any
of the filtering or model generation again, only the model fitting,
based on the filtered data.frame

### 2.3.1 Model fitting

Then in the following code each model is fitted and AICc is calculated
to order it

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Model fitting code for vegetation abundance

</summary>

``` r
vegetation_data = read_excel("Pinpoint-data til ordination_minus_Lønstrup.xlsx")%>% 
    janitor::clean_names()
  # Remocing the first column (ID) in vegetation-dataset and cheking there are no NAs
  
vegetation_data_no_ID = vegetation_data
  
Fs <- AICcPermanova::fit_models(all_forms = Filtered, 
                 env_data = METADATAS, 
                 veg_data = vegetation_data_no_ID, 
                 ncores = 21, 
                 method = "bray",
                logfile = "PlantAbund.txt",
                strata = "area"
                 )

saveRDS(Fs, "FSVegAbund.rds")

Fs <- Fs %>% arrange(AICc)

saveRDS(Fs, "FSVegAbund.rds")
openxlsx::write.xlsx(Fs, "FSVegAbund.xlsx")
```

</details>

As seen in table <a href="#tab:SummaryVegAbund">2.3</a> there are 70
models within 2 AICc of each other, you can see there how many times a
variable has been selected

| Variable          | Full_Akaike_Adjusted_RSq | Number_of_models |
|:------------------|-------------------------:|-----------------:|
| p_h\_water        |                    0.082 |               55 |
| water_content     |                    0.037 |               48 |
| ec                |                    0.025 |               39 |
| oc_beregnet       |                    0.013 |               26 |
| ammonium          |                    0.011 |               25 |
| fine_silt         |                    0.005 |               14 |
| shannon_bak       |                    0.008 |               13 |
| wr                |                    0.004 |               11 |
| coarse_silt_sand  |                    0.003 |                8 |
| finesilt_and_clay |                    0.002 |                6 |
| nitrat_nitrit     |                    0.002 |                6 |
| clay              |                    0.001 |                4 |
| dexter_n          |                    0.000 |                2 |

Table 2.3: Number of selected models were variables are present and
their Akaike Weighted R squared for the Marginal effect of the terms

Now we can see the top models that have a delta AICc within 2 of the
best model in table <a href="#tab:BestVegAbundModels">2.4</a> if
expanded

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Show table of selected models for vegetation presence absence

</summary>

| form                                                                   | max_vif |    AICc |   k |   N | p_h\_water |    ec | water_content |    wr | ammonium | nitrat_nitrit | oc_beregnet |  clay | fine_silt | coarse_silt_sand | shannon_bak | finesilt_and_clay | dexter_n | DeltaAICc | AICWeight |
|:-----------------------------------------------------------------------|--------:|--------:|----:|----:|-----------:|------:|--------------:|------:|---------:|--------------:|------------:|------:|----------:|-----------------:|------------:|------------------:|---------:|----------:|----------:|
| Distance \~ p_h\_water + water_content                                 |   1.010 | -59.241 |   3 |  52 |      0.158 |    NA |         0.069 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     0.000 |     0.029 |
| Distance \~ p_h\_water + ec + water_content                            |   2.269 | -59.038 |   4 |  52 |      0.062 | 0.031 |         0.060 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     0.204 |     0.026 |
| Distance \~ p_h\_water + water_content + ammonium                      |   2.295 | -58.885 |   4 |  52 |      0.161 |    NA |         0.049 |    NA |    0.029 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     0.356 |     0.024 |
| Distance \~ p_h\_water + water_content + fine_silt                     |   1.706 | -58.739 |   4 |  52 |      0.127 |    NA |         0.057 |    NA |       NA |            NA |          NA |    NA |     0.027 |               NA |          NA |                NA |       NA |     0.503 |     0.022 |
| Distance \~ p_h\_water + ec + water_content + wr                       |   4.161 | -58.731 |   5 |  52 |      0.056 | 0.040 |         0.056 | 0.030 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     0.511 |     0.022 |
| Distance \~ p_h\_water + ec + water_content + ammonium                 |   2.510 | -58.651 |   5 |  52 |      0.065 | 0.031 |         0.045 |    NA |    0.029 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     0.591 |     0.021 |
| Distance \~ p_h\_water + ec + water_content + oc_beregnet              |   3.146 | -58.581 |   5 |  52 |      0.049 | 0.039 |         0.043 |    NA |       NA |            NA |       0.028 |    NA |        NA |               NA |          NA |                NA |       NA |     0.660 |     0.020 |
| Distance \~ p_h\_water + water_content + nitrat_nitrit                 |   1.424 | -58.433 |   4 |  52 |      0.128 |    NA |         0.066 |    NA |       NA |         0.022 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     0.808 |     0.019 |
| Distance \~ ec + water_content + oc_beregnet + shannon_bak             |   2.671 | -58.408 |   5 |  52 |         NA | 0.072 |         0.041 |    NA |       NA |            NA |       0.039 |    NA |        NA |               NA |       0.047 |                NA |       NA |     0.834 |     0.019 |
| Distance \~ p_h\_water + water_content + finesilt_and_clay             |   1.542 | -58.390 |   4 |  52 |      0.133 |    NA |         0.055 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.022 |       NA |     0.851 |     0.019 |
| Distance \~ p_h\_water + water_content + wr                            |   2.381 | -58.313 |   4 |  52 |      0.082 |    NA |         0.063 | 0.021 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     0.929 |     0.018 |
| Distance \~ p_h\_water + water_content + oc_beregnet                   |   2.325 | -58.281 |   4 |  52 |      0.124 |    NA |         0.045 |    NA |       NA |            NA |       0.020 |    NA |        NA |               NA |          NA |                NA |       NA |     0.960 |     0.018 |
| Distance \~ p_h\_water + ec + water_content + fine_silt                |   2.583 | -58.202 |   5 |  52 |      0.064 | 0.027 |         0.056 |    NA |       NA |            NA |          NA |    NA |     0.023 |               NA |          NA |                NA |       NA |     1.039 |     0.017 |
| Distance \~ p_h\_water + water_content + coarse_silt_sand              |   1.885 | -58.136 |   4 |  52 |      0.147 |    NA |         0.044 |    NA |       NA |            NA |          NA |    NA |        NA |            0.018 |          NA |                NA |       NA |     1.105 |     0.016 |
| Distance \~ ec + water_content + shannon_bak                           |   2.143 | -58.083 |   4 |  52 |         NA | 0.057 |         0.085 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.048 |                NA |       NA |     1.158 |     0.016 |
| Distance \~ p_h\_water + ec + water_content + ammonium + oc_beregnet   |   3.146 | -58.012 |   6 |  52 |      0.048 | 0.038 |         0.042 |    NA |    0.027 |            NA |       0.026 |    NA |        NA |               NA |          NA |                NA |       NA |     1.230 |     0.015 |
| Distance \~ p_h\_water + ammonium + fine_silt                          |   1.277 | -58.006 |   4 |  52 |      0.143 |    NA |            NA |    NA |    0.046 |            NA |          NA |    NA |     0.036 |               NA |          NA |                NA |       NA |     1.235 |     0.015 |
| Distance \~ p_h\_water + ammonium + coarse_silt_sand                   |   1.287 | -57.958 |   4 |  52 |      0.157 |    NA |            NA |    NA |    0.041 |            NA |          NA |    NA |        NA |            0.036 |          NA |                NA |       NA |     1.284 |     0.015 |
| Distance \~ p_h\_water + ec + water_content + wr + ammonium            |   4.187 | -57.952 |   6 |  52 |      0.053 | 0.040 |         0.045 | 0.025 |    0.024 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.290 |     0.015 |
| Distance \~ ec + water_content + wr + shannon_bak                      |   2.941 | -57.950 |   5 |  52 |         NA | 0.061 |         0.056 | 0.033 |       NA |            NA |          NA |    NA |        NA |               NA |       0.045 |                NA |       NA |     1.291 |     0.015 |
| Distance \~ ec + oc_beregnet + shannon_bak                             |   2.406 | -57.949 |   4 |  52 |         NA | 0.072 |            NA |    NA |       NA |            NA |       0.083 |    NA |        NA |               NA |       0.048 |                NA |       NA |     1.292 |     0.015 |
| Distance \~ p_h\_water + water_content + clay                          |   1.252 | -57.948 |   4 |  52 |      0.146 |    NA |         0.058 |    NA |       NA |            NA |          NA | 0.015 |        NA |               NA |          NA |                NA |       NA |     1.293 |     0.015 |
| Distance \~ p_h\_water + ec + oc_beregnet                              |   3.120 | -57.948 |   4 |  52 |      0.048 | 0.040 |            NA |    NA |       NA |            NA |       0.045 |    NA |        NA |               NA |          NA |                NA |       NA |     1.294 |     0.015 |
| Distance \~ p_h\_water + ec + water_content + nitrat_nitrit            |   2.452 | -57.931 |   5 |  52 |      0.062 | 0.027 |         0.062 |    NA |       NA |         0.019 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.310 |     0.015 |
| Distance \~ ec + water_content + ammonium + oc_beregnet + shannon_bak  |   3.193 | -57.917 |   6 |  52 |         NA | 0.073 |         0.041 |    NA |    0.028 |            NA |       0.039 |    NA |        NA |               NA |       0.047 |                NA |       NA |     1.325 |     0.015 |
| Distance \~ p_h\_water + ammonium                                      |   1.046 | -57.904 |   3 |  52 |      0.158 |    NA |            NA |    NA |    0.049 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.338 |     0.015 |
| Distance \~ p_h\_water + ec + ammonium                                 |   2.092 | -57.903 |   4 |  52 |      0.076 | 0.035 |            NA |    NA |    0.044 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.338 |     0.015 |
| Distance \~ p_h\_water + water_content + ammonium + fine_silt          |   3.188 | -57.896 |   5 |  52 |      0.130 |    NA |         0.033 |    NA |    0.023 |            NA |          NA |    NA |     0.021 |               NA |          NA |                NA |       NA |     1.345 |     0.015 |
| Distance \~ p_h\_water + ec + water_content + wr + oc_beregnet         |   4.304 | -57.895 |   6 |  52 |      0.054 | 0.041 |         0.045 | 0.025 |       NA |            NA |       0.023 |    NA |        NA |               NA |          NA |                NA |       NA |     1.347 |     0.015 |
| Distance \~ p_h\_water + ec + water_content + finesilt_and_clay        |   2.564 | -57.881 |   5 |  52 |      0.064 | 0.027 |         0.054 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.018 |       NA |     1.360 |     0.014 |
| Distance \~ p_h\_water + water_content + ammonium + nitrat_nitrit      |   3.186 | -57.877 |   5 |  52 |      0.129 |    NA |         0.045 |    NA |    0.027 |         0.020 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.365 |     0.014 |
| Distance \~ p_h\_water + ammonium + finesilt_and_clay                  |   1.249 | -57.773 |   4 |  52 |      0.146 |    NA |            NA |    NA |    0.046 |            NA |          NA |    NA |        NA |               NA |          NA |             0.033 |       NA |     1.468 |     0.014 |
| Distance \~ p_h\_water + water_content + ammonium + oc_beregnet        |   2.874 | -57.758 |   5 |  52 |      0.128 |    NA |         0.044 |    NA |    0.027 |            NA |       0.019 |    NA |        NA |               NA |          NA |                NA |       NA |     1.483 |     0.014 |
| Distance \~ p_h\_water + water_content + dexter_n                      |   1.041 | -57.706 |   4 |  52 |      0.160 |    NA |         0.068 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.012 |     1.536 |     0.013 |
| Distance \~ p_h\_water + water_content + ammonium + finesilt_and_clay  |   2.857 | -57.650 |   5 |  52 |      0.137 |    NA |         0.033 |    NA |    0.024 |            NA |          NA |    NA |        NA |               NA |          NA |             0.017 |       NA |     1.592 |     0.013 |
| Distance \~ p_h\_water + ec + water_content + coarse_silt_sand         |   2.883 | -57.649 |   5 |  52 |      0.058 | 0.028 |         0.044 |    NA |       NA |            NA |          NA |    NA |        NA |            0.015 |          NA |                NA |       NA |     1.592 |     0.013 |
| Distance \~ p_h\_water + oc_beregnet                                   |   1.211 | -57.609 |   3 |  52 |      0.131 |    NA |            NA |    NA |       NA |            NA |       0.045 |    NA |        NA |               NA |          NA |                NA |       NA |     1.632 |     0.013 |
| Distance \~ ec + oc_beregnet + fine_silt + shannon_bak                 |   2.695 | -57.602 |   5 |  52 |         NA | 0.069 |            NA |    NA |       NA |            NA |       0.081 |    NA |     0.030 |               NA |       0.049 |                NA |       NA |     1.640 |     0.013 |
| Distance \~ ec + water_content + ammonium + shannon_bak                |   2.433 | -57.563 |   5 |  52 |         NA | 0.058 |         0.058 |    NA |    0.027 |            NA |          NA |    NA |        NA |               NA |       0.050 |                NA |       NA |     1.678 |     0.012 |
| Distance \~ ec + water_content + oc_beregnet                           |   1.944 | -57.563 |   4 |  52 |         NA | 0.114 |         0.042 |    NA |       NA |            NA |       0.041 |    NA |        NA |               NA |          NA |                NA |       NA |     1.679 |     0.012 |
| Distance \~ p_h\_water + water_content + wr + ammonium                 |   2.480 | -57.561 |   5 |  52 |      0.081 |    NA |         0.049 | 0.016 |    0.024 |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.680 |     0.012 |
| Distance \~ p_h\_water + ec + water_content + clay                     |   2.462 | -57.545 |   5 |  52 |      0.062 | 0.029 |         0.056 |    NA |       NA |            NA |          NA | 0.014 |        NA |               NA |          NA |                NA |       NA |     1.696 |     0.012 |
| Distance \~ p_h\_water + oc_beregnet + fine_silt                       |   1.491 | -57.536 |   4 |  52 |      0.114 |    NA |            NA |    NA |       NA |            NA |       0.040 |    NA |     0.034 |               NA |          NA |                NA |       NA |     1.705 |     0.012 |
| Distance \~ p_h\_water + coarse_silt_sand                              |   1.030 | -57.534 |   3 |  52 |      0.163 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |        NA |            0.043 |          NA |                NA |       NA |     1.708 |     0.012 |
| Distance \~ p_h\_water + water_content + oc_beregnet + fine_silt       |   2.857 | -57.517 |   5 |  52 |      0.114 |    NA |         0.035 |    NA |       NA |            NA |       0.017 |    NA |     0.024 |               NA |          NA |                NA |       NA |     1.725 |     0.012 |
| Distance \~ p_h\_water + water_content + ammonium + coarse_silt_sand   |   3.352 | -57.504 |   5 |  52 |      0.151 |    NA |         0.028 |    NA |    0.026 |            NA |          NA |    NA |        NA |            0.015 |          NA |                NA |       NA |     1.737 |     0.012 |
| Distance \~ p_h\_water + water_content + wr + oc_beregnet              |   3.225 | -57.449 |   5 |  52 |      0.082 |    NA |         0.046 | 0.023 |       NA |            NA |       0.023 |    NA |        NA |               NA |          NA |                NA |       NA |     1.792 |     0.012 |
| Distance \~ p_h\_water + ec + ammonium + oc_beregnet                   |   3.133 | -57.449 |   5 |  52 |      0.046 | 0.040 |            NA |    NA |    0.028 |            NA |       0.028 |    NA |        NA |               NA |          NA |                NA |       NA |     1.792 |     0.012 |
| Distance \~ ec + ammonium + oc_beregnet + shannon_bak                  |   2.454 | -57.440 |   5 |  52 |         NA | 0.073 |            NA |    NA |    0.028 |            NA |       0.057 |    NA |        NA |               NA |       0.046 |                NA |       NA |     1.801 |     0.012 |
| Distance \~ p_h\_water + ec + ammonium + coarse_silt_sand              |   2.880 | -57.432 |   5 |  52 |      0.063 | 0.027 |            NA |    NA |    0.041 |            NA |          NA |    NA |        NA |            0.028 |          NA |                NA |       NA |     1.809 |     0.012 |
| Distance \~ p_h\_water + water_content + ammonium + clay               |   2.434 | -57.425 |   5 |  52 |      0.151 |    NA |         0.039 |    NA |    0.027 |            NA |          NA | 0.014 |        NA |               NA |          NA |                NA |       NA |     1.817 |     0.011 |
| Distance \~ p_h\_water + ec + water_content + dexter_n                 |   2.280 | -57.425 |   5 |  52 |      0.064 | 0.031 |         0.060 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.012 |     1.817 |     0.011 |
| Distance \~ p_h\_water + water_content + clay + fine_silt              |   4.079 | -57.424 |   5 |  52 |      0.125 |    NA |         0.060 |    NA |       NA |            NA |          NA | 0.016 |     0.027 |               NA |          NA |                NA |       NA |     1.817 |     0.011 |
| Distance \~ ec + oc_beregnet + coarse_silt_sand + shannon_bak          |   2.649 | -57.408 |   5 |  52 |         NA | 0.070 |            NA |    NA |       NA |            NA |       0.057 |    NA |        NA |            0.027 |       0.049 |                NA |       NA |     1.833 |     0.011 |
| Distance \~ ec + oc_beregnet + shannon_bak + finesilt_and_clay         |   2.649 | -57.408 |   5 |  52 |         NA | 0.070 |            NA |    NA |       NA |            NA |       0.081 |    NA |        NA |               NA |       0.049 |             0.027 |       NA |     1.833 |     0.011 |
| Distance \~ p_h\_water + ec + oc_beregnet + fine_silt                  |   3.277 | -57.401 |   5 |  52 |      0.046 | 0.033 |            NA |    NA |       NA |            NA |       0.045 |    NA |     0.027 |               NA |          NA |                NA |       NA |     1.841 |     0.011 |
| Distance \~ ec + water_content + wr + ammonium + shannon_bak           |   2.943 | -57.399 |   6 |  52 |         NA | 0.062 |         0.052 | 0.033 |    0.027 |            NA |          NA |    NA |        NA |               NA |       0.045 |                NA |       NA |     1.842 |     0.011 |
| Distance \~ ec + water_content + oc_beregnet + fine_silt + shannon_bak |   3.225 | -57.392 |   6 |  52 |         NA | 0.066 |         0.032 |    NA |       NA |            NA |       0.038 |    NA |     0.021 |               NA |       0.045 |                NA |       NA |     1.850 |     0.011 |
| Distance \~ p_h\_water + ec + water_content + wr + nitrat_nitrit       |   4.186 | -57.390 |   6 |  52 |      0.055 | 0.037 |         0.055 | 0.027 |       NA |         0.016 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.852 |     0.011 |
| Distance \~ p_h\_water + ec + ammonium + fine_silt                     |   2.551 | -57.357 |   5 |  52 |      0.077 | 0.026 |            NA |    NA |    0.044 |            NA |          NA |    NA |     0.027 |               NA |          NA |                NA |       NA |     1.885 |     0.011 |
| Distance \~ ec + water_content + wr + oc_beregnet + shannon_bak        |   3.438 | -57.357 |   6 |  52 |         NA | 0.066 |         0.043 | 0.020 |       NA |            NA |       0.027 |    NA |        NA |               NA |       0.047 |                NA |       NA |     1.885 |     0.011 |
| Distance \~ p_h\_water + ec + water_content + oc_beregnet + fine_silt  |   3.436 | -57.354 |   6 |  52 |      0.045 | 0.033 |         0.034 |    NA |       NA |            NA |       0.023 |    NA |     0.018 |               NA |          NA |                NA |       NA |     1.888 |     0.011 |
| Distance \~ p_h\_water + water_content + nitrat_nitrit + oc_beregnet   |   2.339 | -57.346 |   5 |  52 |      0.109 |    NA |         0.042 |    NA |       NA |         0.022 |       0.019 |    NA |        NA |               NA |          NA |                NA |       NA |     1.895 |     0.011 |
| Distance \~ p_h\_water + ec + water_content + ammonium + nitrat_nitrit |   3.265 | -57.326 |   6 |  52 |      0.064 | 0.027 |         0.044 |    NA |    0.027 |         0.017 |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.915 |     0.011 |
| Distance \~ p_h\_water + ec                                            |   1.841 | -57.312 |   3 |  52 |      0.088 | 0.040 |            NA |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.929 |     0.011 |
| Distance \~ p_h\_water + oc_beregnet + coarse_silt_sand                |   1.685 | -57.305 |   4 |  52 |      0.116 |    NA |            NA |    NA |       NA |            NA |       0.032 |    NA |        NA |            0.031 |          NA |                NA |       NA |     1.936 |     0.011 |
| Distance \~ p_h\_water + oc_beregnet + finesilt_and_clay               |   1.453 | -57.305 |   4 |  52 |      0.116 |    NA |            NA |    NA |       NA |            NA |       0.039 |    NA |        NA |               NA |          NA |             0.031 |       NA |     1.936 |     0.011 |
| Distance \~ p_h\_water + ec + water_content + ammonium + fine_silt     |   3.233 | -57.288 |   6 |  52 |      0.066 | 0.027 |         0.034 |    NA |    0.022 |            NA |          NA |    NA |     0.016 |               NA |          NA |                NA |       NA |     1.953 |     0.011 |
| Distance \~ p_h\_water + fine_silt                                     |   1.161 | -57.248 |   3 |  52 |      0.151 |    NA |            NA |    NA |       NA |            NA |          NA |    NA |     0.039 |               NA |          NA |                NA |       NA |     1.994 |     0.011 |
| Distance \~ ec + water_content + wr                                    |   1.310 | -57.245 |   4 |  52 |         NA | 0.066 |         0.059 | 0.036 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.996 |     0.011 |

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
METADATAS <- list.files(pattern = "PERMANOVA_BACTERIA_", full.names = T) |> 
  purrr::map(read_excel) |> 
  purrr::reduce(full_join) |> 
  janitor::clean_names()

AllForms <- AICcPermanova::make_models(vars = c("p_h_water", "ec", 
"water_content", "wr", "ammonium", "nitrat_nitrit", "oc_beregnet", 
"clay", "fine_silt", "coarse_silt_sand", "shannon_veg", "finesilt_and_clay", 
"dexter_n")
, ncores = 21, k = floor(nrow(METADATAS)/10))
  

saveRDS(AllForms, "AllFormsBacterialAbund.rds")
openxlsx::write.xlsx(AllForms, "AllFormsBacterialAbund.xlsx")
```

</details>

This generate up to 2,380 models to evaluate, which can be downloaded as
an excel file
[here](https://github.com/Sustainscapes/AICcPermanova/raw/master/AllFormsBacterialAbund.xlsx)
an rds
[here](https://github.com/Sustainscapes/AICcPermanova/blob/master/AllFormsBacterialAbund.rds).

## 3.2 Avoiding colinearity

Now before the model fitting, we will filter out models with high degree
of multicolinearity using VIF with 5 as a cut off

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

VIF filter

</summary>

``` r
Filtered <- AICcPermanova::filter_vif(all_forms = AllForms , env_data = METADATAS, ncores = 21)
```

</details>

This is then narrowed down to 1,543 models to evaluate

### 3.2.1 Model fitting

Then in the following code each model is fitted and AICc is calculated
to order it

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Model fitting code bacterial abund

</summary>

``` r

d <- amp_load(
otutable = "AC_otutale_new.txt",
metadata = "Metadata_mix-samples_AC_Danielsen_final.xlsx")
bac.data.subset = amp_subset_taxa(d, "d__Bacteria")
bacterial_data = amp_subset_samples(d, Investigator == "AC")
bacterial_data = as.data.frame(t(bacterial_data$abund)) %>% 
    janitor::clean_names()
METADATAS <- METADATAS[match(rownames(bacterial_data), METADATAS$seq_id),]

Fs <- AICcPermanova::fit_models(all_forms = Filtered, 
                 env_data = METADATAS, 
                 veg_data = bacterial_data, 
                 ncores = 21, 
                 method = "bray",
                logfile = "BactLog.txt",
                strata = "area"
                 )



saveRDS(Fs, "FSBacterialAbund.rds")
openxlsx::write.xlsx(Fs, "FSBacterialAbund.xlsx")
```

</details>

As seen in table <a href="#tab:SummaryBacterialAbund">3.1</a> there are
111 models within 2 AICc of each other, you can see there how many times
a variable has been selected

| Variable          | Full_Akaike_Adjusted_RSq | Number_of_models |
|:------------------|-------------------------:|-----------------:|
| p_h\_water        |                    0.216 |              111 |
| water_content     |                    0.031 |               75 |
| oc_beregnet       |                    0.025 |               66 |
| shannon_veg       |                    0.020 |               63 |
| wr                |                    0.005 |               27 |
| fine_silt         |                    0.006 |               25 |
| dexter_n          |                    0.003 |               21 |
| coarse_silt_sand  |                    0.005 |               20 |
| ec                |                    0.003 |               20 |
| finesilt_and_clay |                    0.005 |               19 |
| clay              |                    0.003 |               17 |
| ammonium          |                    0.001 |               10 |
| nitrat_nitrit     |                    0.000 |                3 |

Table 3.1: Number of selected models were variables are present and
their Akaike Weighted R squared for the Marginal effect of the terms

Now we can see the top models that have a delta AICc within 2 of the
best model in table <a href="#tab:BestBacterialAbundModels">3.2</a> if
expanded

<details style="\&quot;margin-bottom:10px;\&quot;">
<summary>

Show table of selected models for vegetation presence absence

</summary>

| form                                                                                   | max_vif |    AICc |   k |   N | p_h\_water |    ec | water_content |    wr | ammonium | nitrat_nitrit | oc_beregnet |  clay | fine_silt | coarse_silt_sand | shannon_veg | finesilt_and_clay | dexter_n | DeltaAICc | AICWeight |
|:---------------------------------------------------------------------------------------|--------:|--------:|----:|----:|-----------:|------:|--------------:|------:|---------:|--------------:|------------:|------:|----------:|-----------------:|------------:|------------------:|---------:|----------:|----------:|
| Distance \~ p_h\_water + water_content + oc_beregnet + fine_silt + shannon_veg         |   2.877 | -82.167 |   6 |  52 |      0.241 |    NA |         0.034 |    NA |       NA |            NA |       0.032 |    NA |     0.027 |               NA |       0.034 |                NA |       NA |     0.000 |     0.017 |
| Distance \~ p_h\_water + water_content + oc_beregnet + shannon_veg                     |   2.341 | -82.119 |   5 |  52 |      0.270 |    NA |         0.034 |    NA |       NA |            NA |       0.033 |    NA |        NA |               NA |       0.035 |                NA |       NA |     0.048 |     0.017 |
| Distance \~ p_h\_water + water_content + oc_beregnet + coarse_silt_sand + shannon_veg  |   2.601 | -82.057 |   6 |  52 |      0.243 |    NA |         0.032 |    NA |       NA |            NA |       0.033 |    NA |        NA |            0.026 |       0.034 |                NA |       NA |     0.111 |     0.016 |
| Distance \~ p_h\_water + water_content + oc_beregnet + shannon_veg + finesilt_and_clay |   2.601 | -82.057 |   6 |  52 |      0.243 |    NA |         0.032 |    NA |       NA |            NA |       0.033 |    NA |        NA |               NA |       0.034 |             0.026 |       NA |     0.111 |     0.016 |
| Distance \~ p_h\_water + water_content + oc_beregnet + clay + shannon_veg              |   2.370 | -81.743 |   6 |  52 |      0.253 |    NA |         0.031 |    NA |       NA |            NA |       0.033 | 0.023 |        NA |               NA |       0.035 |                NA |       NA |     0.424 |     0.014 |
| Distance \~ p_h\_water + water_content + wr + shannon_veg                              |   2.658 | -81.664 |   5 |  52 |      0.164 |    NA |         0.059 | 0.028 |       NA |            NA |          NA |    NA |        NA |               NA |       0.036 |                NA |       NA |     0.503 |     0.014 |
| Distance \~ p_h\_water + water_content + fine_silt + shannon_veg                       |   1.717 | -81.649 |   5 |  52 |      0.246 |    NA |         0.061 |    NA |       NA |            NA |          NA |    NA |     0.028 |               NA |       0.033 |                NA |       NA |     0.518 |     0.013 |
| Distance \~ p_h\_water + oc_beregnet + coarse_silt_sand + shannon_veg                  |   1.722 | -81.593 |   5 |  52 |      0.243 |    NA |            NA |    NA |       NA |            NA |       0.047 |    NA |        NA |            0.029 |       0.034 |                NA |       NA |     0.575 |     0.013 |
| Distance \~ p_h\_water + oc_beregnet + shannon_veg + finesilt_and_clay                 |   1.514 | -81.593 |   5 |  52 |      0.243 |    NA |            NA |    NA |       NA |            NA |       0.060 |    NA |        NA |               NA |       0.034 |             0.029 |       NA |     0.575 |     0.013 |
| Distance \~ p_h\_water + water_content + wr + coarse_silt_sand + shannon_veg           |   2.669 | -81.591 |   6 |  52 |      0.159 |    NA |         0.041 | 0.028 |       NA |            NA |          NA |    NA |        NA |            0.026 |       0.036 |                NA |       NA |     0.577 |     0.013 |
| Distance \~ p_h\_water + water_content + shannon_veg                                   |   1.055 | -81.568 |   4 |  52 |      0.301 |    NA |         0.067 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.033 |                NA |       NA |     0.599 |     0.013 |
| Distance \~ p_h\_water + oc_beregnet + fine_silt + shannon_veg                         |   1.555 | -81.534 |   5 |  52 |      0.241 |    NA |            NA |    NA |       NA |            NA |       0.060 |    NA |     0.028 |               NA |       0.034 |                NA |       NA |     0.633 |     0.013 |
| Distance \~ p_h\_water + water_content + wr + fine_silt + shannon_veg                  |   2.880 | -81.512 |   6 |  52 |      0.163 |    NA |         0.048 | 0.026 |       NA |            NA |          NA |    NA |     0.026 |               NA |       0.036 |                NA |       NA |     0.655 |     0.013 |
| Distance \~ p_h\_water + water_content + oc_beregnet + fine_silt                       |   2.857 | -81.508 |   5 |  52 |      0.249 |    NA |         0.033 |    NA |       NA |            NA |       0.031 |    NA |     0.028 |               NA |          NA |                NA |       NA |     0.659 |     0.013 |
| Distance \~ p_h\_water + water_content + shannon_veg + finesilt_and_clay               |   1.546 | -81.502 |   5 |  52 |      0.252 |    NA |         0.059 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.033 |             0.027 |       NA |     0.665 |     0.012 |
| Distance \~ p_h\_water + water_content + wr + shannon_veg + finesilt_and_clay          |   2.802 | -81.495 |   6 |  52 |      0.162 |    NA |         0.048 | 0.027 |       NA |            NA |          NA |    NA |        NA |               NA |       0.036 |             0.025 |       NA |     0.673 |     0.012 |
| Distance \~ p_h\_water + water_content + coarse_silt_sand + shannon_veg                |   1.886 | -81.477 |   5 |  52 |      0.272 |    NA |         0.046 |    NA |       NA |            NA |          NA |    NA |        NA |            0.026 |       0.033 |                NA |       NA |     0.690 |     0.012 |
| Distance \~ p_h\_water + oc_beregnet + shannon_veg                                     |   1.271 | -81.468 |   4 |  52 |      0.282 |    NA |            NA |    NA |       NA |            NA |       0.066 |    NA |        NA |               NA |       0.034 |                NA |       NA |     0.700 |     0.012 |
| Distance \~ p_h\_water + water_content + oc_beregnet                                   |   2.325 | -81.451 |   4 |  52 |      0.278 |    NA |         0.034 |    NA |       NA |            NA |       0.032 |    NA |        NA |               NA |          NA |                NA |       NA |     0.716 |     0.012 |
| Distance \~ p_h\_water + water_content + oc_beregnet + coarse_silt_sand                |   2.588 | -81.391 |   5 |  52 |      0.252 |    NA |         0.032 |    NA |       NA |            NA |       0.032 |    NA |        NA |            0.027 |          NA |                NA |       NA |     0.776 |     0.012 |
| Distance \~ p_h\_water + water_content + oc_beregnet + finesilt_and_clay               |   2.588 | -81.391 |   5 |  52 |      0.252 |    NA |         0.032 |    NA |       NA |            NA |       0.031 |    NA |        NA |               NA |          NA |             0.027 |       NA |     0.776 |     0.012 |
| Distance \~ p_h\_water + oc_beregnet + clay + shannon_veg                              |   1.409 | -81.374 |   5 |  52 |      0.254 |    NA |            NA |    NA |       NA |            NA |       0.061 | 0.026 |        NA |               NA |       0.034 |                NA |       NA |     0.793 |     0.012 |
| Distance \~ p_h\_water + ec + water_content + oc_beregnet + shannon_veg                |   3.208 | -81.312 |   6 |  52 |      0.110 | 0.019 |         0.035 |    NA |       NA |            NA |       0.032 |    NA |        NA |               NA |       0.035 |                NA |       NA |     0.855 |     0.011 |
| Distance \~ p_h\_water + water_content + wr + clay + shannon_veg                       |   2.704 | -81.268 |   6 |  52 |      0.161 |    NA |         0.050 | 0.028 |       NA |            NA |          NA | 0.023 |        NA |               NA |       0.036 |                NA |       NA |     0.900 |     0.011 |
| Distance \~ p_h\_water + water_content + wr + oc_beregnet + shannon_veg                |   3.610 | -81.228 |   6 |  52 |      0.161 |    NA |         0.034 | 0.018 |       NA |            NA |       0.023 |    NA |        NA |               NA |       0.034 |                NA |       NA |     0.940 |     0.011 |
| Distance \~ p_h\_water + water_content + clay + shannon_veg                            |   1.252 | -81.171 |   5 |  52 |      0.270 |    NA |         0.059 |    NA |       NA |            NA |          NA | 0.023 |        NA |               NA |       0.033 |                NA |       NA |     0.996 |     0.011 |
| Distance \~ p_h\_water + water_content + fine_silt                                     |   1.706 | -81.169 |   4 |  52 |      0.256 |    NA |         0.058 |    NA |       NA |            NA |          NA |    NA |     0.029 |               NA |          NA |                NA |       NA |     0.998 |     0.011 |
| Distance \~ p_h\_water + water_content                                                 |   1.010 | -81.083 |   3 |  52 |      0.312 |    NA |         0.063 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.084 |     0.010 |
| Distance \~ p_h\_water + water_content + oc_beregnet + shannon_veg + dexter_n          |   2.478 | -81.057 |   6 |  52 |      0.269 |    NA |         0.034 |    NA |       NA |            NA |       0.034 |    NA |        NA |               NA |       0.030 |                NA |    0.016 |     1.111 |     0.010 |
| Distance \~ p_h\_water + water_content + oc_beregnet + clay                            |   2.353 | -81.055 |   5 |  52 |      0.261 |    NA |         0.031 |    NA |       NA |            NA |       0.032 | 0.023 |        NA |               NA |          NA |                NA |       NA |     1.112 |     0.010 |
| Distance \~ p_h\_water + water_content + finesilt_and_clay                             |   1.542 | -81.033 |   4 |  52 |      0.262 |    NA |         0.056 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.027 |       NA |     1.134 |     0.010 |
| Distance \~ p_h\_water + oc_beregnet + coarse_silt_sand                                |   1.685 | -81.015 |   4 |  52 |      0.251 |    NA |            NA |    NA |       NA |            NA |       0.045 |    NA |        NA |            0.029 |          NA |                NA |       NA |     1.153 |     0.010 |
| Distance \~ p_h\_water + oc_beregnet + finesilt_and_clay                               |   1.453 | -81.015 |   4 |  52 |      0.251 |    NA |            NA |    NA |       NA |            NA |       0.056 |    NA |        NA |               NA |          NA |             0.029 |       NA |     1.153 |     0.010 |
| Distance \~ p_h\_water + water_content + coarse_silt_sand                              |   1.885 | -80.997 |   4 |  52 |      0.282 |    NA |         0.045 |    NA |       NA |            NA |          NA |    NA |        NA |            0.027 |          NA |                NA |       NA |     1.170 |     0.010 |
| Distance \~ p_h\_water + ec + water_content + wr + shannon_veg                         |   4.606 | -80.984 |   6 |  52 |      0.087 | 0.020 |         0.055 | 0.028 |       NA |            NA |          NA |    NA |        NA |               NA |       0.036 |                NA |       NA |     1.183 |     0.010 |
| Distance \~ p_h\_water + oc_beregnet + fine_silt                                       |   1.491 | -80.978 |   4 |  52 |      0.249 |    NA |            NA |    NA |       NA |            NA |       0.056 |    NA |     0.028 |               NA |          NA |                NA |       NA |     1.189 |     0.010 |
| Distance \~ p_h\_water + oc_beregnet                                                   |   1.211 | -80.931 |   3 |  52 |      0.293 |    NA |            NA |    NA |       NA |            NA |       0.061 |    NA |        NA |               NA |          NA |                NA |       NA |     1.237 |     0.009 |
| Distance \~ p_h\_water + water_content + wr                                            |   2.381 | -80.921 |   4 |  52 |      0.179 |    NA |         0.056 | 0.026 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.246 |     0.009 |
| Distance \~ p_h\_water + ec + water_content + shannon_veg                              |   2.269 | -80.913 |   5 |  52 |      0.126 | 0.020 |         0.058 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.033 |                NA |       NA |     1.254 |     0.009 |
| Distance \~ p_h\_water + ec + water_content + fine_silt + shannon_veg                  |   2.588 | -80.898 |   6 |  52 |      0.129 | 0.019 |         0.058 |    NA |       NA |            NA |          NA |    NA |     0.027 |               NA |       0.033 |                NA |       NA |     1.270 |     0.009 |
| Distance \~ p_h\_water + water_content + clay + fine_silt + shannon_veg                |   4.127 | -80.832 |   6 |  52 |      0.247 |    NA |         0.064 |    NA |       NA |            NA |          NA | 0.019 |     0.024 |               NA |       0.034 |                NA |       NA |     1.335 |     0.009 |
| Distance \~ p_h\_water + water_content + oc_beregnet + dexter_n                        |   2.476 | -80.814 |   5 |  52 |      0.278 |    NA |         0.033 |    NA |       NA |            NA |       0.033 |    NA |        NA |               NA |          NA |                NA |    0.020 |     1.353 |     0.009 |
| Distance \~ p_h\_water + water_content + wr + coarse_silt_sand                         |   2.393 | -80.808 |   5 |  52 |      0.174 |    NA |         0.040 | 0.025 |       NA |            NA |          NA |    NA |        NA |            0.026 |          NA |                NA |       NA |     1.359 |     0.009 |
| Distance \~ p_h\_water + oc_beregnet + clay                                            |   1.353 | -80.779 |   4 |  52 |      0.263 |    NA |            NA |    NA |       NA |            NA |       0.057 | 0.026 |        NA |               NA |          NA |                NA |       NA |     1.388 |     0.009 |
| Distance \~ p_h\_water + water_content + ammonium + oc_beregnet + shannon_veg          |   2.894 | -80.752 |   6 |  52 |      0.264 |    NA |         0.033 |    NA |    0.013 |            NA |       0.031 |    NA |        NA |               NA |       0.034 |                NA |       NA |     1.415 |     0.009 |
| Distance \~ p_h\_water + water_content + oc_beregnet + fine_silt + dexter_n            |   2.899 | -80.730 |   6 |  52 |      0.248 |    NA |         0.034 |    NA |       NA |            NA |       0.032 |    NA |     0.027 |               NA |          NA |                NA |    0.019 |     1.437 |     0.008 |
| Distance \~ p_h\_water + water_content + wr + fine_silt                                |   2.595 | -80.723 |   5 |  52 |      0.178 |    NA |         0.047 | 0.023 |       NA |            NA |          NA |    NA |     0.025 |               NA |          NA |                NA |       NA |     1.444 |     0.008 |
| Distance \~ p_h\_water + water_content + wr + finesilt_and_clay                        |   2.515 | -80.693 |   5 |  52 |      0.177 |    NA |         0.046 | 0.024 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.025 |       NA |     1.474 |     0.008 |
| Distance \~ p_h\_water + water_content + clay                                          |   1.252 | -80.688 |   4 |  52 |      0.280 |    NA |         0.056 |    NA |       NA |            NA |          NA | 0.023 |        NA |               NA |          NA |                NA |       NA |     1.479 |     0.008 |
| Distance \~ p_h\_water + ec + water_content + shannon_veg + finesilt_and_clay          |   2.567 | -80.671 |   6 |  52 |      0.128 | 0.019 |         0.056 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.033 |             0.025 |       NA |     1.496 |     0.008 |
| Distance \~ p_h\_water + water_content + nitrat_nitrit + oc_beregnet + shannon_veg     |   2.353 | -80.664 |   6 |  52 |      0.239 |    NA |         0.033 |    NA |       NA |         0.012 |       0.032 |    NA |        NA |               NA |       0.034 |                NA |       NA |     1.503 |     0.008 |
| Distance \~ p_h\_water + water_content + wr + oc_beregnet                              |   3.225 | -80.645 |   5 |  52 |      0.178 |    NA |         0.033 | 0.018 |       NA |            NA |       0.025 |    NA |        NA |               NA |          NA |                NA |       NA |     1.522 |     0.008 |
| Distance \~ p_h\_water + ec + water_content + oc_beregnet                              |   3.146 | -80.639 |   5 |  52 |      0.110 | 0.018 |         0.034 |    NA |       NA |            NA |       0.030 |    NA |        NA |               NA |          NA |                NA |       NA |     1.528 |     0.008 |
| Distance \~ p_h\_water + wr + oc_beregnet + coarse_silt_sand + shannon_veg             |   3.639 | -80.627 |   6 |  52 |      0.158 |    NA |            NA | 0.017 |       NA |            NA |       0.031 |    NA |        NA |            0.028 |       0.034 |                NA |       NA |     1.540 |     0.008 |
| Distance \~ p_h\_water + wr + oc_beregnet + shannon_veg + finesilt_and_clay            |   3.639 | -80.627 |   6 |  52 |      0.158 |    NA |            NA | 0.017 |       NA |            NA |       0.039 |    NA |        NA |               NA |       0.034 |             0.028 |       NA |     1.540 |     0.008 |
| Distance \~ p_h\_water + wr + oc_beregnet + shannon_veg                                |   3.436 | -80.622 |   5 |  52 |      0.161 |    NA |            NA | 0.018 |       NA |            NA |       0.047 |    NA |        NA |               NA |       0.033 |                NA |       NA |     1.546 |     0.008 |
| Distance \~ p_h\_water + ec + oc_beregnet + fine_silt + shannon_veg                    |   3.277 | -80.604 |   6 |  52 |      0.108 | 0.017 |            NA |    NA |       NA |            NA |       0.055 |    NA |     0.028 |               NA |       0.034 |                NA |       NA |     1.563 |     0.008 |
| Distance \~ p_h\_water + ec + oc_beregnet + shannon_veg                                |   3.187 | -80.598 |   5 |  52 |      0.110 | 0.018 |            NA |    NA |       NA |            NA |       0.055 |    NA |        NA |               NA |       0.034 |                NA |       NA |     1.569 |     0.008 |
| Distance \~ p_h\_water + ec + oc_beregnet + coarse_silt_sand + shannon_veg             |   3.215 | -80.594 |   6 |  52 |      0.108 | 0.017 |            NA |    NA |       NA |            NA |       0.048 |    NA |        NA |            0.028 |       0.034 |                NA |       NA |     1.573 |     0.008 |
| Distance \~ p_h\_water + ec + oc_beregnet + shannon_veg + finesilt_and_clay            |   3.215 | -80.594 |   6 |  52 |      0.108 | 0.017 |            NA |    NA |       NA |            NA |       0.055 |    NA |        NA |               NA |       0.034 |             0.028 |       NA |     1.573 |     0.008 |
| Distance \~ p_h\_water + water_content + oc_beregnet + coarse_silt_sand + dexter_n     |   2.651 | -80.592 |   6 |  52 |      0.249 |    NA |         0.032 |    NA |       NA |            NA |       0.033 |    NA |        NA |            0.025 |          NA |                NA |    0.019 |     1.575 |     0.008 |
| Distance \~ p_h\_water + water_content + oc_beregnet + finesilt_and_clay + dexter_n    |   2.651 | -80.592 |   6 |  52 |      0.249 |    NA |         0.032 |    NA |       NA |            NA |       0.033 |    NA |        NA |               NA |          NA |             0.025 |    0.019 |     1.575 |     0.008 |
| Distance \~ p_h\_water + water_content + oc_beregnet + clay + fine_silt                |   4.362 | -80.549 |   6 |  52 |      0.250 |    NA |         0.035 |    NA |       NA |            NA |       0.031 | 0.017 |     0.022 |               NA |          NA |                NA |       NA |     1.618 |     0.008 |
| Distance \~ p_h\_water + water_content + wr + shannon_veg + dexter_n                   |   2.721 | -80.513 |   6 |  52 |      0.164 |    NA |         0.059 | 0.028 |       NA |            NA |          NA |    NA |        NA |               NA |       0.033 |                NA |    0.015 |     1.654 |     0.008 |
| Distance \~ p_h\_water + water_content + wr + oc_beregnet + fine_silt                  |   3.349 | -80.503 |   6 |  52 |      0.176 |    NA |         0.033 | 0.017 |       NA |            NA |       0.025 |    NA |     0.026 |               NA |          NA |                NA |       NA |     1.664 |     0.008 |
| Distance \~ p_h\_water + wr + oc_beregnet + fine_silt + shannon_veg                    |   3.711 | -80.494 |   6 |  52 |      0.159 |    NA |            NA | 0.016 |       NA |            NA |       0.038 |    NA |     0.026 |               NA |       0.033 |                NA |       NA |     1.674 |     0.008 |
| Distance \~ p_h\_water + oc_beregnet + clay + fine_silt + shannon_veg                  |   3.226 | -80.483 |   6 |  52 |      0.242 |    NA |            NA |    NA |       NA |            NA |       0.060 | 0.016 |     0.018 |               NA |       0.034 |                NA |       NA |     1.684 |     0.008 |
| Distance \~ p_h\_water + wr + oc_beregnet + clay + shannon_veg                         |   3.520 | -80.478 |   6 |  52 |      0.158 |    NA |            NA | 0.018 |       NA |            NA |       0.042 | 0.026 |        NA |               NA |       0.034 |                NA |       NA |     1.689 |     0.007 |
| Distance \~ p_h\_water + water_content + shannon_veg + dexter_n                        |   1.104 | -80.477 |   5 |  52 |      0.293 |    NA |         0.066 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.029 |                NA |    0.015 |     1.691 |     0.007 |
| Distance \~ p_h\_water + oc_beregnet + shannon_veg + dexter_n                          |   1.324 | -80.466 |   5 |  52 |      0.280 |    NA |            NA |    NA |       NA |            NA |       0.066 |    NA |        NA |               NA |       0.029 |                NA |    0.016 |     1.701 |     0.007 |
| Distance \~ p_h\_water + water_content + wr + clay                                     |   2.420 | -80.462 |   5 |  52 |      0.176 |    NA |         0.048 | 0.025 |       NA |            NA |          NA | 0.023 |        NA |               NA |          NA |                NA |       NA |     1.705 |     0.007 |
| Distance \~ p_h\_water + water_content + wr + oc_beregnet + coarse_silt_sand           |   3.302 | -80.449 |   6 |  52 |      0.175 |    NA |         0.032 | 0.017 |       NA |            NA |       0.024 |    NA |        NA |            0.026 |          NA |                NA |       NA |     1.718 |     0.007 |
| Distance \~ p_h\_water + water_content + wr + oc_beregnet + finesilt_and_clay          |   3.302 | -80.449 |   6 |  52 |      0.175 |    NA |         0.032 | 0.017 |       NA |            NA |       0.025 |    NA |        NA |               NA |          NA |             0.026 |       NA |     1.718 |     0.007 |
| Distance \~ p_h\_water + water_content + ammonium + shannon_veg                        |   2.434 | -80.446 |   5 |  52 |      0.285 |    NA |         0.039 |    NA |    0.015 |            NA |          NA |    NA |        NA |               NA |       0.032 |                NA |       NA |     1.721 |     0.007 |
| Distance \~ p_h\_water + ec + water_content + coarse_silt_sand + shannon_veg           |   2.886 | -80.445 |   6 |  52 |      0.122 | 0.016 |         0.046 |    NA |       NA |            NA |          NA |    NA |        NA |            0.023 |       0.033 |                NA |       NA |     1.722 |     0.007 |
| Distance \~ p_h\_water + water_content + fine_silt + shannon_veg + dexter_n            |   1.956 | -80.443 |   6 |  52 |      0.246 |    NA |         0.059 |    NA |       NA |            NA |          NA |    NA |     0.027 |               NA |       0.029 |                NA |    0.015 |     1.724 |     0.007 |
| Distance \~ p_h\_water + ec + water_content + fine_silt                                |   2.583 | -80.442 |   5 |  52 |      0.131 | 0.019 |         0.056 |    NA |       NA |            NA |          NA |    NA |     0.028 |               NA |          NA |                NA |       NA |     1.725 |     0.007 |
| Distance \~ p_h\_water + ec + water_content                                            |   2.269 | -80.435 |   4 |  52 |      0.127 | 0.020 |         0.055 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.733 |     0.007 |
| Distance \~ p_h\_water + oc_beregnet + coarse_silt_sand + shannon_veg + dexter_n       |   2.088 | -80.412 |   6 |  52 |      0.243 |    NA |            NA |    NA |       NA |            NA |       0.045 |    NA |        NA |            0.027 |       0.030 |                NA |    0.015 |     1.756 |     0.007 |
| Distance \~ p_h\_water + oc_beregnet + shannon_veg + finesilt_and_clay + dexter_n      |   1.531 | -80.412 |   6 |  52 |      0.243 |    NA |            NA |    NA |       NA |            NA |       0.057 |    NA |        NA |               NA |       0.030 |             0.027 |    0.015 |     1.756 |     0.007 |
| Distance \~ p_h\_water + ammonium + oc_beregnet + coarse_silt_sand + shannon_veg       |   2.144 | -80.380 |   6 |  52 |      0.242 |    NA |            NA |    NA |    0.014 |            NA |       0.031 |    NA |        NA |            0.029 |       0.034 |                NA |       NA |     1.788 |     0.007 |
| Distance \~ p_h\_water + ammonium + oc_beregnet + shannon_veg + finesilt_and_clay      |   1.848 | -80.380 |   6 |  52 |      0.242 |    NA |            NA |    NA |    0.014 |            NA |       0.035 |    NA |        NA |               NA |       0.034 |             0.029 |       NA |     1.788 |     0.007 |
| Distance \~ p_h\_water + oc_beregnet + fine_silt + shannon_veg + dexter_n              |   1.563 | -80.363 |   6 |  52 |      0.241 |    NA |            NA |    NA |       NA |            NA |       0.059 |    NA |     0.027 |               NA |       0.030 |                NA |    0.015 |     1.805 |     0.007 |
| Distance \~ p_h\_water + ec + water_content + oc_beregnet + fine_silt                  |   3.436 | -80.358 |   6 |  52 |      0.105 | 0.015 |         0.031 |    NA |       NA |            NA |       0.027 |    NA |     0.025 |               NA |          NA |                NA |       NA |     1.810 |     0.007 |
| Distance \~ p_h\_water + wr + coarse_silt_sand + shannon_veg                           |   2.652 | -80.347 |   5 |  52 |      0.159 |    NA |            NA | 0.034 |       NA |            NA |          NA |    NA |        NA |            0.044 |       0.035 |                NA |       NA |     1.820 |     0.007 |
| Distance \~ p_h\_water + water_content + dexter_n                                      |   1.041 | -80.345 |   4 |  52 |      0.307 |    NA |         0.063 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |    0.019 |     1.822 |     0.007 |
| Distance \~ p_h\_water + water_content + wr + ammonium + shannon_veg                   |   2.706 | -80.335 |   6 |  52 |      0.164 |    NA |         0.038 | 0.027 |    0.013 |            NA |          NA |    NA |        NA |               NA |       0.036 |                NA |       NA |     1.833 |     0.007 |
| Distance \~ p_h\_water + oc_beregnet + dexter_n                                        |   1.305 | -80.332 |   4 |  52 |      0.292 |    NA |            NA |    NA |       NA |            NA |       0.063 |    NA |        NA |               NA |          NA |                NA |    0.021 |     1.836 |     0.007 |
| Distance \~ p_h\_water + ammonium + oc_beregnet + fine_silt + shannon_veg              |   1.853 | -80.318 |   6 |  52 |      0.240 |    NA |            NA |    NA |    0.014 |            NA |       0.035 |    NA |     0.029 |               NA |       0.034 |                NA |       NA |     1.849 |     0.007 |
| Distance \~ p_h\_water + water_content + fine_silt + dexter_n                          |   1.956 | -80.295 |   5 |  52 |      0.256 |    NA |         0.057 |    NA |       NA |            NA |          NA |    NA |     0.027 |               NA |          NA |                NA |    0.018 |     1.872 |     0.007 |
| Distance \~ p_h\_water + ec + water_content + clay + shannon_veg                       |   2.463 | -80.291 |   6 |  52 |      0.126 | 0.018 |         0.056 |    NA |       NA |            NA |          NA | 0.021 |        NA |               NA |       0.033 |                NA |       NA |     1.877 |     0.007 |
| Distance \~ p_h\_water + water_content + shannon_veg + finesilt_and_clay + dexter_n    |   1.924 | -80.284 |   6 |  52 |      0.252 |    NA |         0.056 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |       0.029 |             0.026 |    0.014 |     1.884 |     0.007 |
| Distance \~ p_h\_water + water_content + clay + fine_silt                              |   4.079 | -80.281 |   5 |  52 |      0.257 |    NA |         0.060 |    NA |       NA |            NA |          NA | 0.018 |     0.023 |               NA |          NA |                NA |       NA |     1.886 |     0.007 |
| Distance \~ p_h\_water + ec + oc_beregnet + clay + shannon_veg                         |   3.205 | -80.278 |   6 |  52 |      0.109 | 0.016 |            NA |    NA |       NA |            NA |       0.055 | 0.024 |        NA |               NA |       0.034 |                NA |       NA |     1.889 |     0.007 |
| Distance \~ p_h\_water + water_content + coarse_silt_sand + shannon_veg + dexter_n     |   2.191 | -80.260 |   6 |  52 |      0.271 |    NA |         0.043 |    NA |       NA |            NA |          NA |    NA |        NA |            0.025 |       0.030 |                NA |    0.014 |     1.907 |     0.007 |
| Distance \~ p_h\_water + ammonium + oc_beregnet + shannon_veg                          |   1.833 | -80.232 |   5 |  52 |      0.282 |    NA |            NA |    NA |    0.014 |            NA |       0.036 |    NA |        NA |               NA |       0.034 |                NA |       NA |     1.935 |     0.007 |
| Distance \~ p_h\_water + ec + water_content + finesilt_and_clay                        |   2.564 | -80.231 |   5 |  52 |      0.129 | 0.019 |         0.054 |    NA |       NA |            NA |          NA |    NA |        NA |               NA |          NA |             0.026 |       NA |     1.936 |     0.007 |
| Distance \~ p_h\_water + water_content + nitrat_nitrit + shannon_veg                   |   1.433 | -80.230 |   5 |  52 |      0.251 |    NA |         0.063 |    NA |       NA |         0.013 |          NA |    NA |        NA |               NA |       0.033 |                NA |       NA |     1.937 |     0.007 |
| Distance \~ p_h\_water + water_content + ammonium + fine_silt + shannon_veg            |   3.194 | -80.226 |   6 |  52 |      0.244 |    NA |         0.034 |    NA |    0.012 |            NA |          NA |    NA |     0.025 |               NA |       0.032 |                NA |       NA |     1.942 |     0.007 |
| Distance \~ p_h\_water + water_content + ammonium + coarse_silt_sand + shannon_veg     |   3.363 | -80.222 |   6 |  52 |      0.263 |    NA |         0.030 |    NA |    0.014 |            NA |          NA |    NA |        NA |            0.025 |       0.033 |                NA |       NA |     1.945 |     0.007 |
| Distance \~ p_h\_water + ec + water_content + wr                                       |   4.161 | -80.217 |   5 |  52 |      0.091 | 0.020 |         0.052 | 0.025 |       NA |            NA |          NA |    NA |        NA |               NA |          NA |                NA |       NA |     1.950 |     0.007 |
| Distance \~ p_h\_water + water_content + ammonium + oc_beregnet                        |   2.874 | -80.217 |   5 |  52 |      0.275 |    NA |         0.034 |    NA |    0.014 |            NA |       0.030 |    NA |        NA |               NA |          NA |                NA |       NA |     1.950 |     0.007 |
| Distance \~ p_h\_water + wr + oc_beregnet                                              |   3.109 | -80.216 |   4 |  52 |      0.178 |    NA |            NA | 0.019 |       NA |            NA |       0.048 |    NA |        NA |               NA |          NA |                NA |       NA |     1.952 |     0.007 |
| Distance \~ p_h\_water + water_content + oc_beregnet + clay + dexter_n                 |   2.476 | -80.213 |   6 |  52 |      0.258 |    NA |         0.031 |    NA |       NA |            NA |       0.033 | 0.021 |        NA |               NA |          NA |                NA |    0.019 |     1.954 |     0.007 |
| Distance \~ p_h\_water + oc_beregnet + coarse_silt_sand + dexter_n                     |   2.084 | -80.206 |   5 |  52 |      0.249 |    NA |            NA |    NA |       NA |            NA |       0.044 |    NA |        NA |            0.027 |          NA |                NA |    0.019 |     1.961 |     0.007 |
| Distance \~ p_h\_water + oc_beregnet + finesilt_and_clay + dexter_n                    |   1.486 | -80.206 |   5 |  52 |      0.249 |    NA |            NA |    NA |       NA |            NA |       0.055 |    NA |        NA |               NA |          NA |             0.027 |    0.019 |     1.961 |     0.007 |
| Distance \~ p_h\_water + nitrat_nitrit + oc_beregnet + shannon_veg                     |   1.562 | -80.197 |   5 |  52 |      0.239 |    NA |            NA |    NA |       NA |         0.013 |       0.063 |    NA |        NA |               NA |       0.034 |                NA |       NA |     1.970 |     0.007 |
| Distance \~ p_h\_water + ec + water_content + oc_beregnet + coarse_silt_sand           |   3.292 | -80.190 |   6 |  52 |      0.106 | 0.015 |         0.030 |    NA |       NA |            NA |       0.030 |    NA |        NA |            0.023 |          NA |                NA |       NA |     1.977 |     0.006 |
| Distance \~ p_h\_water + ec + water_content + oc_beregnet + finesilt_and_clay          |   3.292 | -80.190 |   6 |  52 |      0.106 | 0.015 |         0.030 |    NA |       NA |            NA |       0.027 |    NA |        NA |               NA |          NA |             0.023 |       NA |     1.977 |     0.006 |
| Distance \~ p_h\_water + oc_beregnet + fine_silt + dexter_n                            |   1.510 | -80.185 |   5 |  52 |      0.248 |    NA |            NA |    NA |       NA |            NA |       0.056 |    NA |     0.026 |               NA |          NA |                NA |    0.019 |     1.983 |     0.006 |
| Distance \~ p_h\_water + water_content + wr + oc_beregnet + clay                       |   3.247 | -80.179 |   6 |  52 |      0.174 |    NA |         0.030 | 0.018 |       NA |            NA |       0.025 | 0.023 |        NA |               NA |          NA |                NA |       NA |     1.988 |     0.006 |

Table 3.2: Best models for bacterial abundance

</details>
