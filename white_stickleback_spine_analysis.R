###White Stickleback Spine Length Analysis###

##Mixed model analysis for significant difference in spine length, control for body size



#test 2 

library("nlme")
library("car")

spines.df<-read.csv("EastCoastMorphometrics2015.csv")
spines.df<-as.data.frame(spines.df)
View(spines.df)
spine.length.1<-as.numeric(spines.df$spine.1)
spine.length.2<-as.numeric(spines.df$spine.2)
spine.length.pelvic<-as.numeric(spines.df$spine.pelvic)
standard.length<-as.numeric(spines.df$std.length)
species<-as.character(spines.df$species)
population<-as.character(spines.df$population)

# this is effectively a "paired" t-test (where populations are the "subjects"),  controlling for standard length. Population is a random effect. 
first.spine<-lme(spine.length.1 ~ standard.length * species, random= ~1|population, na.action=na.omit) 
summary(first.spine) 

# type I anova
anova(first.spine) 



#mixed model spine 2
second.spine<-lme(spine.length.2 ~ standard.length * species, random= ~1|population, na.action=na.omit) 
summary(second.spine) 

# type I anova
anova(second.spine) 



#mixed model pelvic spine 
pelvic.spine<-lme(spine.length.pelvic ~ standard.length * species, random= ~1|population, na.action=na.omit) 
summary(pelvic.spine) 

# type I anova
anova(pelvic.spine) 




##Visuals
library(ggplot2)


#First Spine
Spine<-ggplot(spines.df, aes(standard.length, spines.df$spine.1, colour=factor(species)), na.rm=T) + xlab("Standard Body Length") + ylab("Spine Length (cm)")
Spine + geom_point(na.rm=T) + stat_smooth(method=lm, na.rm=T)            

#second Spine
Spine + geom_point(aes(standard.length, spines.df$spine.2), na.rm=T) + stat_smooth(method=lm, na.rm=T)               

#pelvic Spine 
Spine+ geom_point(aes(standard.length, spines.df$spine.pelvic), na.rm=T) + stat_smooth(method=lm, na.rm=T)



#dplyr
library("dplyr")

head(spines.df)

# mean spine.1 and std.length for each population
df.spine1<-spines.df[complete.cases(spines.df[,"spine.1"]),] #first, create data frame with NAs removed for spine 

mean.spine.1<- df.spine1 %>% 
  group_by(population) %>%
  summarise(mean.sp.1=mean(spine.1), std.length=mean(std.length), species=mean(species))
qplot(std.length, mean.sp.1, data=mean.spine.1, colour=factor(species)) + xlab("Average Standard Length (cm)") + ylab("Average Spine 1 Length (cm)")
#this shows species as numbers (not sure how to change or if possible) 1=Both 2=Common 3=White. Colours are same as other figures


#control for standard length with residuals 
spine1.mod <- lm(spine.1~std.length, data=df.spine1)
df.spine1$spine1.resid <- residuals(spine1.mod)


with(df.spine1, plot(std.length, spine.1)) #example residual plots
with(df.spine1, plot(std.length, spine1.resid))

#Boxplot spine 1 length residuals 
BSR<-ggplot(df.spine1, aes(population, spine1.resid, colour=factor(species))) + xlab("Population") + ylab("Residuals for Spine 1 Length")
BSR + geom_boxplot() #By Population
BSR + geom_boxplot(aes(species, spine1.resid)) #By Species



#AGAIN WITH SPINE 2...
# mean spine.2 and std.length for each population

df.spine2<-spines.df[complete.cases(spines.df[,"spine.2"]),] #first, create data frame with NAs removed for spine 

mean.spine.2<- df.spine2 %>% 
  group_by(population) %>%
  summarise(mean.sp.2=mean(spine.2), std.length=mean(std.length), species=mean(species))

qplot(std.length, mean.sp.2, data=mean.spine.2, colour=factor(species)) + xlab("Average Standard Length (cm)") + ylab("Average Spine 2 Length (cm)")

# control for std.length
spine2.mod <- lm(spine.2~std.length, data=df.spine2)
df.spine2$spine2.resid <- residuals(spine2.mod)


#Boxplot spine 2 length residuals 
BSR2<-ggplot(df.spine2, aes(population, spine2.resid, colour=factor(species))) + xlab("Population") + ylab("Residuals for Spine 2 Length")
BSR2 + geom_boxplot() #By Population
BSR2 + geom_boxplot(aes(species, spine2.resid)) #By Species



#AGAIN WITH PELVIC SPINE...
#mean spine.pelvic and std.length for each population

df.spine.pelvic<-spines.df[complete.cases(spines.df[,"spine.pelvic"]),] #first, create data frame with NAs removed for spine 

mean.spine.pelvic<- df.spine.pelvic %>% 
  group_by(population) %>%
  summarise(mean.sp.plv=mean(spine.pelvic), std.length=mean(std.length), species=mean(species))

qplot(std.length, mean.sp.plv, data=mean.spine.pelvic, colour=factor(species)) + xlab("Average Standard Length (cm)") + ylab("Average Pelvic Spine Length (cm)")

# control for std.length
spine.plv.mod <- lm(spine.pelvic~std.length, data=df.spine.pelvic)
df.spine.pelvic$spine.plv.resid <- residuals(spine.plv.mod)


#Boxplot pelvic spine length residuals 
BSRP<-ggplot(df.spine.pelvic, aes(population, spine.plv.resid, colour=factor(species))) + xlab("Population") + ylab("Residuals for Pelvic Spine Length")
BSRP + geom_boxplot() #By Population
BSRP + geom_boxplot(aes(species, spine.plv.resid)) #By Species
