# Clear memory
rm(list = ls())

#Packages
library(car)
library(MCMCglmm)
library(dplyr)
library(lme4)

#Model input
load('modelinput.RData')
#first check the colinearity
vif(lmer(Npositive ~ scale(Temp)+scale(Prec)+
           scale(Clutch_MEAN)+scale(Maximum.longevity)+
           scale(Mass)+Habitat+as.factor(Migration)+Ntested+
           (1|animal), data=preva_trait, REML=FALSE))
#1.01 - 2.24, OK

#Recode predictors - habitat and migration
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


model_list = list()

for (i in pathogen){
  df = preva_trait[preva_trait$Pathogen == i,]
  model_list[[i]] = MCMCglmm(cbind(Npositive,Ntested-Npositive) ~ scale(Temp)+scale(Prec)+
                               scale(Clutch_MEAN)+scale(Maximum.longevity)+
                               scale(Mass)+Migration+Habitat, 
                             random=~animal, 
                             pedigree = phy,
                             family ="multinomial2",
                             data=df,
                             nitt=133000, 
                             burnin=3000, 
                             thin=100)
}
save(model_list, file = "models.Rdata")

i=1
summary(model_list[[i]])
par(mar=c(2,2,2,2))
plot(model_list[[i]]$Sol)
plot(model_list[[i]]$VCV)
R2(model_list[[i]])






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
                     scale(Mass)+Migration+Habitat, 
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


load('models.Rdata')
save(model_list, validate_result_list, validate_df_list, file = "models.Rdata")


#Prediction----
load('modelinput.Rdata')
load('models.Rdata')
load("data4predict.RData")

#Habitat classification
predict2 = predict2 %>% 
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


#Add clutch size for one species that was not in the available databases
predict2[predict2$species =="Cuculus saturatus",]$Clutch_MEAN = 1
#exclude five observations for two species without clutch size (all non-native)
predict2 = na.omit(predict2)
save.image("all4predict.Rdata")
#Now all ready for prediction
i = pathogen[2]
#Set predicted prevalence to range from 0 to 100
Npositive = 0
Ntested = 100
predictresult = cbind(predict2, predict.MCMCglmm(model_list[[i]], newdata=predict2))
