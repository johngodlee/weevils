#Weevil damage 17/11/15
library(spatstat)
library(ggplot2)

setwd("~/Dropbox/Weevils/Analysis")
read.csv(file="damage.csv")
damage<-read.csv(file="damage.csv")
head(damage)
tail(damage)
damage$block<-as.factor(damage$block)

#Plots
plot(total_damage~population, data=damage)
plot(prev_damage~population, data=damage)
plot(curr_damage~population, data=damage)

#Basic ANOVA
pop_total_lm<-lm(total_damage~population, data=damage)
summary.aov(pop_total_lm)

#ANOVA with block as random
pop_total_wblockr_lmer<-lmer(total_damage~population + (1|block), data=damage)
summary
#Spatial analysis
damage_bl1<-damage[which(damage$block==1),]
damage_bl2<-damage[which(damage$block==2),]
damage_bl3<-damage[which(damage$block==3),]
damage_bl4<-damage[which(damage$block==4),]

rbPal <- colorRampPalette(c('red','blue'))
damage_bl1$col <- rbPal(10)[as.numeric(cut(damage_bl1$total_damage,breaks = 10))]
plot(damage_bl1$x_coord,damage_bl1$y_coord,pch=15,col = damage_bl1$col, cex=2)

qplot(x_coord,y_coord, data=damage_bl1, colour=total_damage) + scale_colour_gradient(low="red", high="blue")

symbols(damage_bl1$x_coord, damage_bl1$y_coord, circles=damage_bl1$total_damage)


#spatstat
x<-c(0,22)
y<-c(0,9)
spat<-ppp(damage_bl1$x_coord, damage_bl1$y_coord, x, y)
plot(spat)
dist_mat<-pairdist(spat)
