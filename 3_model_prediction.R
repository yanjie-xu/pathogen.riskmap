# Clear memory
rm(list = ls())

#Packages
library(car)
library(MCMCglmm)
library(dplyr)
library(lme4)

load('modelinput.Rdata')
#load('models.Rdata')
load("data4predict.RData")

#Habitat classification
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
                                 `Human Modified`='HumanModified'))


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
                                 `Human Modified`='HumanModified'))


#Add clutch size for one species that was not in the available databases
predict2[predict2$species =="Cuculus saturatus",]$Clutch_MEAN = 1
#exclude five observations for two species without clutch size (all non-native)
predict2 = na.omit(predict2)
predict2$Ntested = 100
predict2$Npositive = 0
#predict2 = predict2[,c("animal","Mass","Habitat","Migration","Clutch_MEAN",
#                             "Maximum.longevity","Temp","Prec","Ntested",
#                             "Npositive")]
#save(preva_trait, phy, predict2, pathogen, file = "all4predict.Rdata")
#Set predicted prevalence to range from 0 to 100

predictw = predict2[predict2$season == "non-breeding",]
predictb = predict2[predict2$season == "breeding",]
predictw.td = predictw %>%                            
  group_by(cell50x50) %>%
  summarise(count = n_distinct(animal))``
gridw = predictw.td[predictw.td$count > 5,]$cell50x50
predictw = predictw[predictw$cell50x50 %in% gridw,]


predictb.td = predictb %>%                            
  group_by(cell50x50) %>%
  summarise(count = n_distinct(animal))
gridb = predictb.td[predictb.td$count > 5,]$cell50x50
predictb = predictb[predictb$cell50x50 %in% gridb,]


#save(preva_trait, phy, predictb, predictw, pathogen, file = "all4predict.Rdata")
load("all4predict.Rdata")


i = pathogen[1]
df = preva_trait[preva_trait$Pathogen == i,]
#df = df[,c("animal","Mass","Habitat","Migration","Clutch_MEAN",
#                             "Maximum.longevity","Temp","Prec","Ntested",
#                             "Npositive")]
model = MCMCglmm(cbind(Npositive, Ntested-Npositive) ~ Temp + Prec +
                   Clutch_MEAN + Maximum.longevity +
                   Mass + Migration+Habitat, 
                 random=~animal, 
                 pedigree = phy,
                 family ="multinomial2",
                 data=df,
                 nitt=133000, 
                 burnin=3000, 
                 thin=100)

#predictw = predict2[predict2$season == "non-breeding",]
#predictb = predict2[predict2$season == "breeding",]
#Now all ready for prediction

predictwresult = predictw[0,]
predictwresult['predict.preva'] = double()

for (i in 1:38){
  if (i < 38) {
    predictw1 = predictw[(1+(i-1)*10000):(i*10000),]
  } else{
    predictw1 = predictw[(1+(i-1)*10000):(nrow(predictw)),]
  }
  
  predictwresult1 = cbind(predictw1, predict.MCMCglmm(model, newdata=predictw1))
  predictwresult = rbind(predictwresult,predictwresult1)
  gc()
}


predictbresult = predictb[0,]
predictbresult['predict.preva'] = double()

for (i in 1:58){
  if (i < 58) {
    predictb1 = predictb[(1+(i-1)*10000):(i*10000),]
  } else{
    predictb1 = predictb[(1+(i-1)*10000):(nrow(predictb)),]
  }
  
  predictbresult1 = cbind(predictb1, predict.MCMCglmm(model, newdata=predictb1))
  predictbresult = rbind(predictbresult,predictbresult1)
  gc()
}

predictwresult = cbind(predictw, predict.MCMCglmm(model, newdata=predictw))
predictbresult = cbind(predictb, predict.MCMCglmm(model, newdata=predictb))
save(predictwresult,predictbresult, file = paste(pathogen[1],".Rdata",sep=""))