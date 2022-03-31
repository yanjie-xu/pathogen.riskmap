# Clear memory
rm(list = ls())

#Packages
library(car)
library(MCMCglmm)
library(dplyr)
library(lme4)
#Calculate R squared for MCMCglmm
R2 <- function(mod){
  fixed_eff <- colMeans(mod$Sol)
  fixed_var_comp <- var(as.vector(fixed_eff %*% t(mod$X)))
  all_randoms <- colMeans(mod$VCV)
  residual <- all_randoms[["units"]]
  random_var_comp <- sum(all_randoms) - residual
  R2 <- (fixed_var_comp + random_var_comp)/(sum(all_randoms) + fixed_var_comp)
  round(R2,3)
}

#Prepare model input (skip to the next) ----

#Data
load('alldata.RData')
# Input data: select one tree from the 1000 treeset
phy=tree$tree_8551
plot(phy,type='fan')
#To be solved: How to get the averaged tree from the 1000 trees?

#Input data: subset prevalence data by pathogen
#Only pathogenic taxa with over 100 observations are included in this analysis
#And now bird only
preva = preva[preva$Class == "Bird",]

pathogen = preva %>% 
  group_by(Pathogen) %>%
  summarize(Nrecords = length(Season2))

pathogen = pathogen[pathogen$Nrecords >=100,]$Pathogen
#17 pathogens:
pathogen
preva = preva[preva$Pathogen %in% pathogen,]


#Merge prevalence dataset with trait data
preva$Scientific.name = gsub('_', ' ', preva$Scientific.name)

#head(trait)
preva_trait = merge(preva, trait, by.x = 'Scientific.name', by.y = 'prevalence',)
#There are duplicates because of different version names or previously merged species, remove!
preva_trait = preva_trait %>% distinct(ID, .keep_all = TRUE)
summary(preva_trait)
names(preva_trait)[19] = "animal" #Define the column name for data with multi species entries
#Now ready for modelling
#write.csv(preva_trait,'data/modelinput.csv')
save(preva_trait, pathogen, tree, phy, file = "modelinput.RData")



#Prepare predict data set (skip to the next) ----
head(ebba2)
head(birdlife)
#Get the traits
ebba2$ID = 1:nrow(ebba2)
ebba2_trait = merge(ebba2,trait,
                    by.x = 'birdlife_scientific_name', by.y = 'ebba')

ebba2_trait = ebba2_trait %>% distinct(ID, .keep_all = TRUE)

birdlife$ID = 1:nrow(birdlife)
birdlife_trait = merge(birdlife, trait,
                       by.x = 'SCINAME', by.y = 'birdlife')

birdlife_trait = birdlife_trait %>% distinct(ID, .keep_all = TRUE)



#Select useful columns of the variables for modelling
head(ebba2_trait)
ebba2_trait = ebba2_trait[c('birdlife_scientific_name','cell50x50',
                            'birdtree','TipLabel',
                            'Mass','Habitat','Migration',
                            'Clutch_MEAN','Maximum.longevity')]
names(ebba2_trait)[1] = "species"
ebba2_trait$season = "breeding"


head(birdlife_trait)
birdlife_trait = birdlife_trait[c('SCINAME','cell50x50',
                            'birdtree','TipLabel',
                            'Mass','Habitat','Migration',
                            'Clutch_MEAN','Maximum.longevity')]

names(birdlife_trait)[1] = "species"
birdlife_trait$season = "non-breeding"

#combine the two seasonal datasets
predict = rbind(ebba2_trait, birdlife_trait)
#get the climate data of each row
predict2 = merge(predict, grid, by = "cell50x50")
predict2 = predict2 %>% distinct(Lon_2, Lat_2, species, season, .keep_all = TRUE)
names(predict2)[4] = 'animal'

#predict input data ready
save(predict2, file = "data4predict.RData")


#Model2: MCMCglmm----
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




#Salmonella model
#https://stackoverflow.com/questions/51132259/phylogenetic-model-using-multiple-entries-for-each-species
df = preva_trait[preva_trait$Pathogen == pathogen[12],]
pathogen[12]


#Note: the animal column should be the Tiplabel in birdtree (XX_XX) instead of scientific name
model = MCMCglmm(cbind(Npositive,Ntested-Npositive) ~ scale(Temp)+scale(Prec)+
                   scale(Clutch_MEAN)+scale(Maximum.longevity)+
                   scale(Mass)+Migration+Habitat, 
                 random=~animal, 
                 pedigree = phy,
                 family =  "multinomial2",
                 data=df,
                 nitt=133000, 
                 burnin=3000, 
                 thin=100)
summary(model)
#Effective sample size = 1300
#Your effective sample size should be quite high (I usually aim for 1000-2000). 
#More complicated models often require more iterations to achieve a comparable effective sample size.
#The R structure is the residual structure.
#he G structure is the random effects structure

# Random effect
# Plot the posterior distribution as a histogram to check for significance and whether it's been well estimated or not
# Variance cannot be zero, and therefore if the mean value is pushed up against zero your effect is not significant
# The larger the spread of the histogram, the less well estimated the distribution is.
hist(mcmc(model$VCV)[,"animal"])
hist(mcmc(model$Sol))[,"scale(Temp)"]

#Now lets check for model convergence. 
#Fixed effect
par(mar=c(5,1,1,1))
plot(model$Sol)
#To make sure your model has converged, the trace plot should look like a fuzzy caterpillar.
#Looks OK
#Do the same for random effect
plot(model$VCV)
#again mixed well - OK

#Habitat classification probably too detailed for our case? Think about combine the habitat types? 

R2(model)
#0.588

#Check the effect of phylogeny
model_np = MCMCglmm(cbind(Npositive,Ntested-Npositive) ~ scale(Temp)+scale(Prec)+
                   scale(Clutch_MEAN)+scale(Maximum.longevity)+
                   scale(Mass)+Migration+Habitat, 
                 random=~animal, 
    #             pedigree = phy,
                 family =  "multinomial2",
                 data=df,
                 nitt=133000, 
                 burnin=3000, 
                 thin=100)
summary(model_np) #DIC = 7149 (no phylogeny) vs. 7146 (full model)
R2(model_np) #0.452 vs. 0.588 
#Not large but there is an effect of host phylogeny



#Model validation----
#Cross validation: randomly select 70% train and 30% validation
i = 12

df = preva_trait[preva_trait$Pathogen == pathogen[i],]
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
plot(log(validate$pred)~log(validate$Npositive))
summary(lm(pred~Npositive, data = validate))

bayestestR::area_under_curve(validate$pred, validate$Npositive, method = "spline")

#'use the model to predict prevalence for each species * grid combination
load("data4predict.RData")
predict2 = na.omit(predict2)
predict2 = predict2 %>% mutate(Migration = 
                                 dplyr::recode(Migration, 
                                               `1`="Sedentary",
                                               `2`="PartialMigratory",
                                               `3`='FullMigratory'))


predict2$Ntested = 100
predict2$Npositive = 0
a = predict.MCMCglmm(model, newdata=predict21)

predict21 = predict2[1:1000,]
#' Do the prediction
#' We would have to separate the predict data because of the limit of memory
#' 3/12 17:02 -3/14 18:16
#' 
predictresult = predict2[0,]
predictresult['predict.preva'] = double()

for (i in 1:95){
  if (i < 95) {
    predict21 = predict2[(1+(i-1)*10000):(i*10000),]
    summary(predict2)} else{
      predict21 = predict2[(1+(i-1)*10000):(nrow(predict2)),]
    }
  
  predictresult1 = cbind(predict21, predict.MCMCglmm(model, newdata=predict21))
  predictresult = rbind(predictresult,predictresult1)
  gc()
}

#'
#'
#'
#'
#'
# Mapping the risk----
#'First sum up the prevalence of different species occurring in each grid
predict.grid = predict2 %>% 
  group_by(season,cell50x50) %>%
  summarize(pred_preva_grid = sum(pred_preva))
#risk value range: 2 - 3525; mean = 617
ggplot(predict.grid, aes(x = pred_preva_grid)) +
  geom_histogram(aes(color = season, fill = season), 
                 position = "identity", bins = 30, alpha = 0.4) +
  theme_classic()+
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))


#50KM grid
EBBA_grid = sf::st_read("Study_area_grid_WGS84.shp") 
# Merge Name and eoagrid
EBBA_grid = EBBA_grid %>% tidyr::unite("New", Name, eoagrid, na.rm=TRUE)

#separate seasons
pred.grid.breed = predict.grid[predict.grid$season == "breeding",]
pred.grid.winter = predict.grid[predict.grid$season == "non-breeding",]
pred.grid.breed = merge(EBBA_grid,pred.grid.breed,by.x = "New", by.y = "cell50x50")
pred.grid.winter = merge(EBBA_grid,pred.grid.winter,by.x = "New", by.y = "cell50x50")
#'
#'Mapping
glgmap = get_stamenmap(bbox = c(-20, 30, 70, 75), maptype = "terrain-background", 
                       col = "bw", crop = TRUE, zoom = 4)



ggmap(glgmap) + geom_sf(aes(Lon_2, Lat_2,  fill = pred_preva_grid/100), 
                        colour = "transparent", data=pred.grid.breed)+
  scale_fill_gradient2(low = 'blue3', mid = "grey80", high = "red3", midpoint = 5)+
  xlab("Longitude")+
  ylab("Latitude")+
  labs(fill = "Predicted risk")+
  annotate(geom="text", x=-5, y=72, label = '(a) Breeding', 
           size = 4)

ggmap(glgmap) + geom_sf(aes(Lon_2, Lat_2,  fill = pred_preva_grid/100), 
                        colour = "transparent", data=pred.grid.winter)+
  scale_fill_gradient2(low = 'blue3', mid = "grey80", high = "red3", midpoint = 5)+
  xlab("Longitude")+
  ylab("Latitude")+
  labs(fill = "Predicted risk")+
  annotate(geom="text", x=-5, y=72, label = '(b) Non-breeding', 
           size = 4)
plot(pred_preva_grid~Lat_2, data = pred.grid.breed)
#The spatial pattern is ...


