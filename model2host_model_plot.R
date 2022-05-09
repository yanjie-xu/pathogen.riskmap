# Clear memory
rm(list = ls())
library(dplyr)
library(ggplot2)
library(forcats)
library(scales)
library(MCMCglmm)
library(broom.mixed)
library(cowplot)

load('models.Rdata')
load('modelinput.Rdata')



#Figure 1 - Accuracy
count= preva_trait %>%
  count(Pathogen)

accuracy = 
  bind_rows(validate_result_list) %>%
  group_by(pathogen) %>%
  summarise(sd = sd(r2, na.rm = TRUE),
            r2 = mean(r2, na.rm = TRUE))

accuracy_df = 
  bind_rows(validate_result_list)
accuracy_df = na.omit(accuracy_df)

accuracy = merge(accuracy, count, by.x = 'pathogen', by.y = 'Pathogen')
accuracy$group = c("Bacteria","Virus", "Bacteria","Bacteria","Bacteria","Bacteria",
                   "Bacteria","Protozoa","Protozoa","Protozoa","Bacteria",
                   "Bacteria","Virus","Virus","Protozoa","Virus","Virus")
accuracy_df = merge(accuracy_df, count, by.x = 'pathogen', by.y = 'Pathogen')
accuracy_df = merge(accuracy_df, accuracy, by = 'pathogen')
names(accuracy_df)[3:6]= c("r2","n","sd","r2mean")

library(tidyverse)
p = ggplot(accuracy_df, aes(pathogen, r2, color = group)) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.1) + 
  geom_errorbar(aes(ymin = r2-sd, ymax = r2+sd, label=n),data = accuracy, color = "black") +
  geom_point(aes(pathogen, r2), data = accuracy, size = 2, color = "black")+
  theme_classic()+
  ylab(expression(paste("R"^2)))+
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        axis.title.x = element_blank(),
        legend.title=element_blank())+
  geom_text(aes(y = -0.1, label = n), data = accuracy, size = 2.5)+
  scale_color_manual(values=c(muted('red'),"black", muted('blue')))

p + aes(x = fct_reorder2(pathogen,r2,group))

ggsave("fig/Fig1_accuracy.png", dpi = 600)








#pathogen models
#rank by accuracy
Pathogen = accuracy[order(-accuracy$r2),]$pathogen

s = as.data.frame(summary(model_list[[pathogen[1]]])$solutions)
s$pathogen = Pathogen[1]
s$rank = 1
row.names(s)
predictor = c("(Intercept)", "Temperature", "Precipitation","Clutch Size", "Longevity",
              "Body Mass", "Partial Migratory", "Sedentary", "Human Modified", "Open" , "Water")
s$Predictor = predictor
s$Predictor = factor(s$Predictor, levels=predictor)

s$sig = "no"
s[s$pMCMC<=0.05,]$sig = "yes"
names(s)[2:3] = c("lower","upper")

s0 = data.frame(post.mean = 0, lower = 0, upper = 0,
                eff.samp = 0, pMCMC = 0, pathogen = Pathogen, 
                rank  = 1:17, Predictor = "Forest", sig = "ref")
s00 = data.frame(post.mean = 0, lower = 0, upper = 0,
                eff.samp = 0, pMCMC = 0, pathogen = Pathogen, 
                rank  = 1:17, Predictor = "Migratory", sig = "ref")


s = rbind(s,s0,s00)


for (i in 2:17){
  s1 = as.data.frame(summary(model_list[[Pathogen[i]]])$solutions)
  s1$pathogen = Pathogen[i]
  s1$rank = i
  s1$Predictor = predictor
  s1$sig = "no"
  s1[s1$pMCMC<=0.05,]$sig = "yes"
  names(s1)[2:3] = c("lower","upper")
  s = rbind(s,s1)
}

predictor = c("(Intercept)", "Temperature", "Precipitation","Clutch Size", "Longevity",
              "Body Mass", "Migratory","Partial Migratory", "Sedentary", 
              "Forest","Human Modified", "Open" , "Water")
s$Predictor = factor(s$Predictor, levels=predictor)
s = s[s$Predictor!="(Intercept)",]
#s1 = s[s$rank<=12,]

s$pathogen[s$pathogen == "Plasmodium"] = "P[Plasmodium]"
s$pathogen[s$pathogen == "Avian influenza virus"] = "V[Influenza A]"
s$pathogen[s$pathogen == "Trichomonas"] = "P[Trichomonas]"
s$pathogen[s$pathogen == "Campylobacter"] = "B[Campylobacter]"
s$pathogen[s$pathogen == "Leucocytozoon"] = "P[Leucocytozoon]"
s$pathogen[s$pathogen == "Chlamydia"] = "B[Chlamydia]"
s$pathogen[s$pathogen == "Salmonella"] = "B[Salmonella]"
s$pathogen[s$pathogen == "Escherichia"] = "B[Escherichia]"

s$pathogen[s$pathogen == "Haemoproteus"] = "P[Haemoproteus]"
s$pathogen[s$pathogen == "West Nile virus"] = "V[West Nile]"
s$pathogen[s$pathogen == "Coxiella"] = "B[Coxiella]"
s$pathogen[s$pathogen == "Anaplasma"] = "B[Anaplasma]"
s$pathogen[s$pathogen == "Rickettsia"] = "B[Rickettsia]"
s$pathogen[s$pathogen == "Sindbis virus"] = "V[Sindbis]"
s$pathogen[s$pathogen == "Usutu virus"] = "V[Usutu]"
s$pathogen[s$pathogen == "Borrelia"] = "B[Borrelia]"
s$pathogen[s$pathogen == "Tick-borne encephalitis virus"] = "V[Tick-borne encephalitis]"

Pathogen = c("P[Plasmodium]","V[Influenza A]","P[Trichomonas]","B[Campylobacter]",
             "P[Leucocytozoon]","B[Chlamydia]","B[Salmonella]","B[Escherichia]",
             "P[Haemoproteus]","V[West Nile]","B[Coxiella]","B[Anaplasma]",
             "B[Rickettsia]","V[Sindbis]","V[Usutu]","B[Borrelia]","V[Tick-borne encephalitis]")
#Figure 2
clim = ggplot(s[s$Predictor == "Temperature"| s$Predictor == "Precipitation",],
       aes(factor(pathogen, levels=rev(Pathogen)), post.mean, color =sig)) +
  geom_point(size = 2)+
  geom_errorbar(aes(ymin=lower, ymax=upper)) +
  geom_hline(yintercept=0, lty=2, lwd=0.5, colour="grey20")+
  scale_color_manual(values=c("grey70","black"))+
  coord_flip()+
  theme_classic()+
  ylab("Posterior mean")+
  facet_wrap(~ Predictor, ncol = 2)+
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.title.y = element_blank())

trait = ggplot(s[s$Predictor == "Clutch Size"|
           s$Predictor =="Body Mass"| s$Predictor == "Longevity",],
       aes(factor(pathogen, levels=rev(Pathogen)), post.mean, color =sig)) +
  geom_point(size = 2)+
  geom_errorbar(aes(ymin=lower, ymax=upper)) +
  geom_hline(yintercept=0, lty=2, lwd=0.5, colour="grey20")+
  scale_color_manual(values=c("grey70","black"))+
  coord_flip()+
  theme_classic()+
  ylab("Posterior mean")+
  facet_wrap(~ Predictor, ncol = 3)+
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.title.y = element_blank())

habitat = ggplot(s[s$Predictor == "Forest" | s$Predictor =="Human Modified"| 
                   s$Predictor == "Open" | s$Predictor == "Water", ],
               aes(Predictor, post.mean, color =sig, shape = sig)) +
  geom_point(size = 2.5)+
  geom_errorbar(aes(ymin=lower, ymax=upper)) +
  scale_color_manual(values=c(alpha("black",.06),"black","black"))+
  scale_shape_manual(values=c(16, 1, 16))+
  theme_classic()+
  ylab("Posterior mean")+
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.title.x = element_blank())+
  ylim(-5,5)

migra = ggplot(s[s$Predictor == "Migratory" | s$Predictor =="Partial Migratory"| 
                     s$Predictor == "Sedentary", ],
                 aes(Predictor, post.mean, color =sig, shape = sig)) +
  geom_point(size = 2.5)+
  geom_errorbar(aes(ymin=lower, ymax=upper)) +
  scale_color_manual(values=c(alpha("black",.06),"black","black"))+
  scale_shape_manual(values=c(16, 1, 16))+
  theme_classic()+
  ylab("Posterior mean")+
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.title.x = element_blank())+
  ylim(-5,5)

col1 = plot_grid(habitat, migra, align = 'v',nrow = 2)
row1 = plot_grid(col1, clim, align = 'h',rel_widths=c(1.2,2), ncol = 2)

plot_grid(trait,row1, align = 'v',nrow =2)

ggsave("fig/Fig2_estimates.png", width =9, height = 7, units = "in", dpi = 600)  


#Figure S1
ggplot(s, aes(Predictor, post.mean,color =sig, shape = sig)) +
  geom_point(size = 1.5)+
  geom_errorbar(aes(ymin=lower, ymax=upper), 
                lwd=0.7, width=0) +
  geom_hline(yintercept=0, lty=2, lwd=0.5, colour="grey20")+
  scale_color_manual(values=c("grey70","black","black"))+
  scale_shape_manual(values=c(16, 1, 16))+
  coord_flip()+
  theme_classic()+
  ylab("Posterior mean")+
  facet_wrap(~ factor(pathogen, levels=Pathogen), ncol = 4)+
  theme(legend.position = "none",
        strip.background = element_blank())+
  ylim(-7,7)
  
ggsave("fig/FigS1_models.png", dpi = 600)  
  




# Mapping the risk----

load('aivmapb.Rdata')
load('aivmapw.Rdata')
names(predictwresult)[17] = 'pred_preva'
names(predictbresult)[17] = 'pred_preva'
hii = foreign::read.dbf('hii_mean.dbf')
#50KM grid
EBBA_grid = sf::st_read("Study_area_grid_WGS84.shp") 
# Merge Name and eoagrid
EBBA_grid = EBBA_grid %>% tidyr::unite("New", Name, eoagrid, na.rm=TRUE)
EBBA_grid = merge(EBBA_grid,hii, by = "ID")

#'First sum up the prevalence of different species occurring in each grid
predict.grid.w = predictwresult %>% 
  group_by(cell50x50) %>%
  summarize(overallrisk_grid = sum(pred_preva),
            communityrisk = mean(pred_preva),
            sr = n_distinct(animal))

predict.grid.b = predictbresult %>% 
  group_by(cell50x50) %>%
  summarize(overallrisk_grid = sum(pred_preva),
            communityrisk = mean(pred_preva),
            sr = n_distinct(animal))



#separate seasons
pred.grid.breed = merge(EBBA_grid,predict.grid.b,by.x = "New", by.y = "cell50x50")
pred.grid.breed$season = "Breeding"
pred.grid.winter = merge(EBBA_grid,predict.grid.w,by.x = "New", by.y = "cell50x50")
pred.grid.winter$season = "Nonbreeding"
pred.grid = rbind(pred.grid.breed, pred.grid.winter)

worldmap = map_data ("world")
range01 <- function(x){(x-min(x))/(max(x)-min(x))}


p1 = ggplot(worldmap) +
  geom_map(data = worldmap, map = worldmap, aes(x = long, y = lat, map_id = region), fill = "gray80")+
  geom_sf(aes(Lon_2, Lat_2,  fill = overallrisk_grid/100), 
                        colour = "transparent", data=pred.grid)+
  scale_fill_gradient2(low = 'blue3', mid = "grey90", high = "red3", midpoint = 10)+
  labs(fill = "Overall")+
  xlim(-22,75)+
  ylim(30,75)+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title = element_blank())+
  facet_wrap(~season, ncol = 2)
my_tag <- c("Breeding", "Non-breeding")
p1 = egg::tag_facet(p1, 
               x = -Inf, y = Inf, 
               open = "", close = "",
               fontface = 1,
               size = 4,
               tag_pool = my_tag)
  
p2 = ggplot(worldmap) +
  geom_map(data = worldmap, map = worldmap, aes(x = long, y = lat, map_id = region), fill = "gray80")+
  geom_sf(aes(Lon_2, Lat_2,  fill = communityrisk/100), 
          colour = "transparent", data=pred.grid)+
  scale_fill_gradient2(low = 'blue3', mid = "grey90", high = "red3", midpoint = 0.1)+
  labs(fill = "Community")+
  xlim(-22,75)+
  ylim(30,75)+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title = element_blank())+
  facet_wrap(~season, ncol = 2)

p3 = ggplot(worldmap) +
  geom_map(data = worldmap, map = worldmap, aes(x = long, y = lat, map_id = region), fill = "gray80")+
  geom_sf(aes(Lon_2, Lat_2,  fill = MEAN*overallrisk_grid/100), 
          colour = "transparent", data=pred.grid)+
  scale_fill_gradient2(low = 'blue3', mid = "grey90", high = "red3", midpoint = 250)+
  labs(fill = "Spillover")+
  xlim(-22,75)+
  ylim(30,75)+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title = element_blank())+
  facet_wrap(~season, ncol = 2)

#    annotate(geom="text", x=-5, y=72, label = '(b) Non-breeding', 
#             size = 4)+

plot_grid(p1,p2,p3, align = 'v',nrow =3)
#9.2 x 6.97
ggsave("fig/Fig3_mapaiv.png", width =9.2, height = 6.97, units = "in", dpi = 600)  



pd = ggplot(worldmap) +
  geom_map(data = worldmap, map = worldmap, aes(x = long, y = lat, map_id = region), fill = "gray80")+
  geom_sf(aes(Lon_2, Lat_2,  fill = sr), 
          colour = "transparent", data=pred.grid)+
  scale_fill_gradient2(low = muted('blue'), mid = "grey90", high = muted("red"), midpoint = 90)+
  #  scale_fill_distiller(palette = "RdBu")+
  labs(fill = "Diversity")+
  xlim(-22,75)+
  ylim(30,75)+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title = element_blank())+
  facet_wrap(~season, ncol = 2)
my_tag <- c("Breeding", "Non-breeding")
pd = egg::tag_facet(pd, 
                    x = -Inf, y = Inf, 
                    open = "", close = "",
                    fontface = 1,
                    size = 4,
                    tag_pool = my_tag)
#Saving 9.2 x 3.47 in image
ggsave("fig/FigS1_diversity.png", width =9.2, height = 3.47, units = "in", dpi = 600) 



phii = ggplot(worldmap) +
  geom_map(data = worldmap, map = worldmap, aes(x = long, y = lat, map_id = region), fill = "gray80")+
  geom_sf(aes(Lon_2, Lat_2,  fill = MEAN), 
          colour = "transparent", data=pred.grid)+
  scale_fill_gradient2(low = muted('blue'), mid = "grey90", high = muted("red"), midpoint = 20)+
  labs(fill = "HII")+
  xlim(-22,75)+
  ylim(30,75)+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title = element_blank())
#Saving 4.94 x 3.31 in image
ggsave("fig/FigS2_hii.png", width =4.94, height = 3.31, units = "in",  dpi = 600) 


























































#Loop MAPS----

for (i in 3:17){
load(paste(Pathogen[i],".Rdata",sep = ""))
names(predictwresult)[17] = 'pred_preva'
names(predictbresult)[17] = 'pred_preva'
hii = foreign::read.dbf('hii_mean.dbf')
#50KM grid
EBBA_grid = sf::st_read("Study_area_grid_WGS84.shp") 
# Merge Name and eoagrid
EBBA_grid = EBBA_grid %>% tidyr::unite("New", Name, eoagrid, na.rm=TRUE)
EBBA_grid = merge(EBBA_grid,hii, by = "ID")

#'First sum up the prevalence of different species occurring in each grid
predict.grid.w = predictwresult %>% 
  group_by(cell50x50) %>%
  summarize(overallrisk_grid = sum(pred_preva),
            communityrisk = mean(pred_preva),
            sr = n_distinct(animal))

predict.grid.b = predictbresult %>% 
  group_by(cell50x50) %>%
  summarize(overallrisk_grid = sum(pred_preva),
            communityrisk = mean(pred_preva),
            sr = n_distinct(animal))



#separate seasons
pred.grid.breed = merge(EBBA_grid,predict.grid.b,by.x = "New", by.y = "cell50x50")
pred.grid.breed$season = "Breeding"
pred.grid.winter = merge(EBBA_grid,predict.grid.w,by.x = "New", by.y = "cell50x50")
pred.grid.winter$season = "Nonbreeding"
pred.grid = rbind(pred.grid.breed, pred.grid.winter)

worldmap = map_data ("world")


p1 = ggplot(worldmap) +
  geom_map(data = worldmap, map = worldmap, aes(x = long, y = lat, map_id = region), fill = "gray80")+
  geom_sf(aes(Lon_2, Lat_2,  fill = overallrisk_grid/100), 
          colour = "transparent", data=pred.grid)+
  scale_fill_gradient2(low = 'blue3', mid = "grey90", high = "red3", midpoint = 10)+
  labs(fill = "Overall")+
  xlim(-22,75)+
  ylim(30,75)+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title = element_blank())+
  facet_wrap(~season, ncol = 2)
my_tag <- c("Breeding", "Non-breeding")
p1 = egg::tag_facet(p1, 
                    x = -Inf, y = Inf, 
                    open = "", close = "",
                    fontface = 1,
                    size = 4,
                    tag_pool = my_tag)

p2 = ggplot(worldmap) +
  geom_map(data = worldmap, map = worldmap, aes(x = long, y = lat, map_id = region), fill = "gray80")+
  geom_sf(aes(Lon_2, Lat_2,  fill = communityrisk/100), 
          colour = "transparent", data=pred.grid)+
  scale_fill_gradient2(low = 'blue3', mid = "grey90", high = "red3", midpoint = 0.1)+
  labs(fill = "Community")+
  xlim(-22,75)+
  ylim(30,75)+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title = element_blank())+
  facet_wrap(~season, ncol = 2)

p3 = ggplot(worldmap) +
  geom_map(data = worldmap, map = worldmap, aes(x = long, y = lat, map_id = region), fill = "gray80")+
  geom_sf(aes(Lon_2, Lat_2,  fill = MEAN*overallrisk_grid/100), 
          colour = "transparent", data=pred.grid)+
  scale_fill_gradient2(low = 'blue3', mid = "grey90", high = "red3", midpoint = 250)+
  labs(fill = "Spillover")+
  xlim(-22,75)+
  ylim(30,75)+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title = element_blank())+
  facet_wrap(~season, ncol = 2)

#    annotate(geom="text", x=-5, y=72, label = '(b) Non-breeding', 
#             size = 4)+

plot_grid(p1,p2,p3, align = 'v',nrow =3)
#9.2 x 6.97
ggsave(paste("fig/FigS_map",Pathogen[i],".png",sep = ""),
       width =9.2, height = 6.97, units = "in", dpi = 600)  

}











#The spatial pattern is ...
plot(pred_preva_grid~Lat_2, data = pred.grid.breed)
