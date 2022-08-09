
#Mod15 <- adonis2(vegetation_data_no_ID ~ organic_matter_content_LOI, data = env.data, by = "terms")

#Mod16 <- adonis2(vegetation_data_no_ID ~ Clay_NIR, data = env.data, by = "terms")

#Mod17 <- adonis2(vegetation_data_no_ID ~ Fine_silt_NIR, data = env.data, by = "terms")

#Mod18 <- adonis2(vegetation_data_no_ID ~ Coarse_silt_sand_NIR, data = env.data, by = "terms")




######
## OBS: Lønstrup er stadig med skal ændres!!!
######


# Loading packages and data -----------------------------------------------

library(vegan)
library(dplyr)
library(readxl)

"C:/Users/au568758/OneDrive - Aarhus Universitet/Documents/1 Paper #1/Paper-1-Bacteria-and-vegetation-in-grassland-and-heathland/Excel ark/"

vegetation_data = read_excel("C:/Users/au568758/OneDrive - Aarhus Universitet/Documents/1 Paper #1/Paper-1-Bacteria-and-vegetation-in-grassland-and-heathland/Excel ark/Presence_absence_vegetation_AC_Danielsen.xlsx")


meta.data = read_excel("C:/Users/au568758/OneDrive - Aarhus Universitet/Documents/1 Paper #1/Paper-1-Bacteria-and-vegetation-in-grassland-and-heathland/Excel ark//Metadata_mix-samples_AC_Danielsen_final.xlsx")


# Remocing the first column (ID) in vegetation-dataset and cheking there are no NAs

vegetation_data_no_ID = subset(vegetation_data, select = -Plot)

vegetation_data_no_ID

is.na(vegetation_data_no_ID)

# Removing columns from env.data that is not used in the analysis 

env.data = subset(meta.data, select = -c(ID, Barcode, date, number_sub_samples, habitat_code, 
                                         habitat, Latitude, Longitude, type, circle, ppr, areanumber,
                                         `Nitrit_mg/kg`, Species_index, Number_of_species, Ellenberg_L, 
                                         Ellenberg_T, EllenbergK, EllenbergF, EllenbergR, EllenbergN, 
                                         SeqID, Order, habitat_number_for_plotting,
                                         pH_original,
                                         `om/2`,
                                         habitat_code_ordination,
                                         areacolor,
                                         Investigator,
                                         ID_Ordination,
                                         ID_Rarefaction,
                                         SampleName,
                                         ExtractionID,
                                         LibraryID,
                                         EXT_plate,
                                         EXT_wellID,
                                         `Dilution plate ID`,
                                         WellID,
                                         SampleContent,
                                         SampleSite,
                                         SampleDate,
                                         ExtractionConc,
                                         LibraryConc,
                                         ExtractionMethod,
                                         LibraryMethod,
                                         Primer,
                                         Barcode_1,
                                         SampleStorageDate,
                                         ExtractionStorageDate,
                                         LibraryStorageDate))

str(env.data)

env.data = env.data[1:55,]

env.data

is.na(env.data)

# Making models -----------------------------------------------------------


Mod0 <- adonis2(vegetation_data_no_ID ~ 1, data = env.data, by = "terms")

Mod1 <- adonis2(vegetation_data_no_ID ~ habitat_type, data = env.data, by = "terms")

Mod2 <- adonis2(vegetation_data_no_ID ~ area, data = env.data, by = "terms")

Mod3 <- adonis2(vegetation_data_no_ID ~ pH_water, data = env.data, by = "terms")

Mod4 <- adonis2(vegetation_data_no_ID ~ water_content, data = env.data, by = "terms")

Mod5 <- adonis2(vegetation_data_no_ID ~ Ammonium, data = env.data, by = "terms")

Mod6 <- adonis2(vegetation_data_no_ID ~ NitratNitrit, data = env.data, by = "terms")

Mod7 <- adonis2(vegetation_data_no_ID ~ organic_matter_content, data = env.data, by = "terms")

Mod8 <- adonis2(vegetation_data_no_ID ~ Clay, data = env.data, by = "terms")

Mod9 <- adonis2(vegetation_data_no_ID ~ Fine_silt, data = env.data, by = "terms")

Mod10 <- adonis2(vegetation_data_no_ID ~ Coarse_silt_sand, data = env.data, by = "terms")

Mod11 <- adonis2(vegetation_data_no_ID ~ EC, data = env.data, by = "terms")

Mod12 <- adonis2(vegetation_data_no_ID ~ WR, data = env.data, by = "terms")

Mod13 <- adonis2(vegetation_data_no_ID ~ TN, data = env.data, by = "terms")

Mod14 <- adonis2(vegetation_data_no_ID ~ TC, data = env.data, by = "terms")

Mod15 <- adonis2(vegetation_data_no_ID ~ Dexter_n, data = env.data, by = "terms")

Mod16<- adonis2(vegetation_data_no_ID ~ finesilt_and_clay, data = env.data, by = "terms")





AICc.PERMANOVA <- function(adonis.model) {
  
  # check to see if object is an adonis model...
  
  if (!(adonis.model$aov.tab[1,1] >= 1))
    stop("object not output of adonis {vegan} ")
  
  # Ok, now extract appropriate terms from the adonis model
  # Calculating AICc using residual sum of squares (RSS) since I don't think that adonis returns something I can use as a liklihood function...
  
  RSS <- adonis.model$aov.tab[rownames(adonis.model$aov.tab) == "Residuals", "SumsOfSqs"]
  MSE <- adonis.model$aov.tab[rownames(adonis.model$aov.tab) == "Residuals", "MeanSqs"]
  
  k <- ncol(adonis.model$model.matrix)# + 1 # add one for error variance
  
  nn <- nrow(adonis.model$model.matrix)
  
  # AIC : 2*k + n*ln(RSS)
  # AICc: AIC + [2k(k+1)]/(n-k-1)
  
  # based on https://en.wikipedia.org/wiki/Akaike_information_criterion;
  # https://www.researchgate.net/post/What_is_the_AIC_formula;
  # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf
  
  # AIC.g is generalized version of AIC = 2k + n [Ln( 2(pi) RSS/n ) + 1]
  # AIC.pi = k + n [Ln( 2(pi) RSS/(n-k) ) +1],
  
  AIC <- 2*k + nn*log(RSS)
  AICc <- AIC + (2*k*(k + 1))/(nn - k - 1)
  
  output <- data.frame(AICc = AICc,k = k, n = nn)
  
  return(output)   
}

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

DF<- data.frame(Modelo = paste0("Mod", 0:16), AICc = NA)


DF$AICc[1] <- AICc.PERMANOVA2(Mod0)$AICc
DF$AICc[2] <- AICc.PERMANOVA2(Mod1)$AICc
DF$AICc[3] <- AICc.PERMANOVA2(Mod2)$AICc
DF$AICc[4] <-  AICc.PERMANOVA2(Mod3)$AICc
DF$AICc[5] <-  AICc.PERMANOVA2(Mod4)$AICc
DF$AICc[6] <-  AICc.PERMANOVA2(Mod5)$AICc
DF$AICc[7] <-  AICc.PERMANOVA2(Mod6)$AICc
DF$AICc[8] <-  AICc.PERMANOVA2(Mod7)$AICc
DF$AICc[9] <-  AICc.PERMANOVA2(Mod8)$AICc
DF$AICc[10] <-  AICc.PERMANOVA2(Mod9)$AICc
DF$AICc[11] <-  AICc.PERMANOVA2(Mod10)$AICc
DF$AICc[12] <-  AICc.PERMANOVA2(Mod11)$AICc
DF$AICc[13] <-  AICc.PERMANOVA2(Mod12)$AICc
DF$AICc[14] <-  AICc.PERMANOVA2(Mod13)$AICc
DF$AICc[15] <-  AICc.PERMANOVA2(Mod14)$AICc
DF$AICc[16] <-  AICc.PERMANOVA2(Mod15)$AICc
DF$AICc[17] <-  AICc.PERMANOVA2(Mod16)$AICc



Vars <- colnames(env.data)
Dataset <- "vegetation_data_no_ID"

Response = env.data

Fs <- list()


for(x in 1:(length(Vars) + 1)){
  if(x == (length(Vars) + 1)){
    Formulas <- data.frame(Form = paste(Dataset, "~ 1", collapse = ""), AICc = NA)
    Formulas$AICc[j] <- try(AICc.PERMANOVA2(adonis2(as.formula(Formulas$Form), data = Response, by = "terms"))$AICc, silent = T)
  }else{
    Test <- combn(Vars, x, simplify = F)
    Formulas <- data.frame(Form = rep(NA, length(Test)), AICc = rep(NA, length(Test)))
    for(j in 1:length(Test)){
      Temp <- paste(Dataset,"~", paste(Test[[j]], collapse = " + ")) 
      Formulas$Form[j] <- Temp
      Temp <- as.formula(Temp)
      Formulas$AICc[j] <- try(AICc.PERMANOVA2(adonis2(as.formula(Formulas$Form[j]), data = Response, by = "terms"))$AICc, silent = T) 
      gc()
    }
  }
  
  Fs[[x]] <- suppressWarnings(Formulas %>% mutate(AICc = as.numeric(AICc)) %>% dplyr::filter(!is.na(AICc)) %>% arrange(AICc))
  message(paste("finished for", x, "number of variables"))
}

Fs <- suppressWarnings(purrr::reduce(Fs, bind_rows) %>% arrange(AICc))


# Testing the final model -------------------------------------------------

## After selecting your model you can test it 

FinalMod <- adonis2(vegetation_data_no_ID ~ habitat_type + area, data = env.data, by = "terms")

FinalMod

