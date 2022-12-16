# Clear memory
rm(list = ls())

#Packages
library(car)
library(MCMCglmm)
library(dplyr)
library(lme4)

#Model validation----
load('modelinput.RData')
#Cross validation: randomly select 70% train and 30% validation

preva_trait = preva_trait %>% 
  mutate(Migration = dplyr::recode(Migration,
                                   `1`="Sedentary",
                                   `2`="PartialMigratory",
                                   `3`='FullMigratory'))%>% 
  mutate(Habitat = dplyr::recode(Habitat,
                                 `Marine`="Water",
                                 `Coastal`="Water",
                                 `Wetland`='Water',
                                 `Riverine`="Water",
                                 `Grassland`="Open",
                                 `Shrubland`='Forest',
                                 `Forest`="Forest",
                                 `Woodland`='Forest', 
                                 `Rock`="Open",
                                 `Desert`='Open',
                                 `Human Modified`='Human Modified'))


validate_result_list = list()
validate_df_list = list()

for (i in 1:length(pathogen)){
  df = preva_trait[preva_trait$Pathogen == pathogen[i],]
  validate.result = data.frame(id = numeric(), r2 = numeric(), pathogen = character())
  validate.df = df[0,]
  validate.df['pred'] = numeric()
  
  for (j in 1:99){tryCatch({
    validate.result[j,1] = j
    validate.result[j,3] = pathogen[i]
    dft = df[sample(nrow(df), round(0.7*nrow(df))),]
    dfv = subset(df, !(df$X.x %in% dft$X.x))
    
    model = MCMCglmm(cbind(Npositive,Ntested-Npositive) ~ scale(Temp)+scale(Prec)+
                       scale(Clutch_MEAN)+scale(Maximum.longevity)+
                       scale(log(Mass))+Migration+Habitat, 
                     random=~animal, 
                     pedigree = phy,
                     family ="multinomial2",
                     data=dft,
                     nitt=133000, 
                     burnin=3000, 
                     thin=100)
    
    validate = cbind(dfv, predict.MCMCglmm(model, newdata=dfv))
    names(validate)[27] = 'pred'
    validate.df = rbind(validate.df, validate)
    validate.result[j,2] = summary(lm(pred~Npositive, data = validate))$r.squared
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  validate_result_list[[pathogen[i]]] = validate.result
  validate_df_list[[pathogen[i]]] = validate.df
}


save(model_list, validate_result_list, validate_df_list, file = "models2.Rdata")