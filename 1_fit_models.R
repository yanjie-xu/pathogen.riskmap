# Clear memory
rm(list = ls())

#Packages
library(car)
library(MCMCglmm)

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

#Fit models
for (i in pathogen){
  df = preva_trait[preva_trait$Pathogen == i,]
  model_list[[i]] = MCMCglmm(cbind(Npositive,Ntested-Npositive) ~ scale(Temp)+scale(Prec)+
                               scale(Clutch_MEAN)+scale(Maximum.longevity)+
                               scale(log(Mass))+Migration+Habitat, 
                             random=~animal, 
                             pedigree = phy,
                             family ="multinomial2",
                             data=df,
                             nitt=133000, 
                             burnin=3000, 
                             thin=100)
}


save(model_list, file = "models2.Rdata")

#Check model performance
summary(model_list[[i]])
par(mar=c(2,2,2,2))
plot(model_list[[i]]$Sol)
plot(model_list[[i]]$VCV)
