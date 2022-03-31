# Clear memory
rm(list = ls())

#Packages
library(car)
library(MCMCglmm)
library(dplyr)
library(lme4)
#Function for calculate R squared for MCMCglmm
R2 = function(mod){
  fixed_eff = colMeans(mod$Sol)
  fixed_var_comp = var(as.vector(fixed_eff %*% t(mod$X)))
  all_randoms = colMeans(mod$VCV)
  residual = all_randoms[["units"]]
  random_var_comp = sum(all_randoms) - residual
  R2 = (fixed_var_comp + random_var_comp)/(sum(all_randoms) + fixed_var_comp)
  round(R2,3)
}



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





par(mar=c(2,2,2,2))
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
plot(model_list[[i]]$Sol)
plot(model_list[[i]]$VCV)
R2(model_list[[i]])


#Model selection - full model used
#library(MuMIn)
#MCMCglmm.updateable = updateable(MCMCglmm)
#global.model = MCMCglmm.updateable(cbind(Npositive,Ntested-Npositive) ~ scale(Temp)+scale(Prec)+
#                                               scale(Clutch_MEAN)+scale(Maximum.longevity)+
#                                               scale(Mass)+Migration+Habitat, 
#                                             random=~animal, 
#                                             pedigree = phy,
#                                             family ="multinomial2",
#                                             data=df,
#                                             nitt=133000, 
#                                             burnin=3000, 
#                                             thin=100)
#d = dredge(global.model, rank="DIC")
#best = subset(d, DIC==0)

load('models.Rdata')


