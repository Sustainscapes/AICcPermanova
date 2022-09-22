######
## OBS: Lønstrup er stadig med skal ændres!!!
######


# Loading packages and data -----------------------------------------------

library(vegan)
library(dplyr)
library(readxl)
library(parallel)
library(foreach)
library(doParallel)

meta.data = read_excel("PERMANOVA_VEGETATION_ clay_silt_sand_OC_and_others_AC_Danielsen.xlsx") %>% 
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

Forms <- list()

Models <- for(i in 1:floor(nrow(env.data)/10)){
  Test <- combn(Vars, i, simplify = F)
  cl <- makeCluster(21)
  registerDoParallel(cl)
  Formulas <- foreach(j = 1:length(Test), .combine = "rbind")%dopar%{
    Dataset <- "vegetation_data_no_ID"
    DF <- data.frame(Form = NA, AICc = NA)
    Temp <- paste(Dataset,"~", paste(Test[[j]], collapse = " + ")) 
    DF$Form <- Temp
    gc()
    DF
  }
  stopCluster(cl)
  message(paste(i, "of", floor(nrow(env.data)/10), "ready", Sys.time()))
  Forms[[i]] <- Formulas
}

Forms <- Forms %>% 
  purrr::reduce(bind_rows) 

Dataset <- "vegetation_data_no_ID"

NullMod <- data.frame(Form = paste(Dataset, "~ 1", collapse = ""), AICc = NA)


Forms <- Forms %>% bind_rows(NullMod)

cl <- makeCluster(21)
registerDoParallel(cl)

Fs <- foreach(x = 1:nrow(Forms), .packages = c("vegan", "dplyr", "tidyr", "readxl"), .combine = bind_rows) %dopar% {
  meta.data = read_excel("PERMANOVA_VEGETATION_ clay_silt_sand_OC_and_others_AC_Danielsen.xlsx") %>% 
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
  Temp <- Formulas[x,]
  Temp$AICc <-  try(AICc.PERMANOVA2(adonis2(as.formula(Forms$Form[x]), data = Response, by = "terms"))$AICc, silent = T)
  Rs <- broom::tidy(adonis2(as.formula(Forms$Form[x]), data = Response, by = "terms")) %>% dplyr::filter(!(term %in% c("Residual", "Total"))) %>% dplyr::select(term, R2) %>%  pivot_wider(names_from = term, values_from = R2)
  if((x %% 100) == 0){
    sink("log.txt", append = T)
    cat(paste("finished for", x, "number of variables", Sys.time(), nrow(Formulas)))
    cat("\n")
    sink()
  }
  
  bind_cols(Temp, Rs)
}

stopCluster(cl)

saveRDS(Fs, "FS.rds")

Fs <- Fs %>% arrange(AICc)

saveRDS(Fs, "FS.rds")

# Testing the final model -------------------------------------------------

## After selecting your model you can test it 

FinalMod <- adonis2(vegetation_data_no_ID ~ habitat_type + area, data = env.data, by = "terms")

FinalMod

