###Eggs/testes size analysis###


##load raw data##

######Google sheet method not working#######
# google sheets API
library("devtools")
library("geomorph")
library("ggplot2")
library("candisc")
library("MASS")
library("dplyr")

devtools::install_github("jennybc/googlesheets")

library("googlesheets")

# list sheets and authenticate with server (will open browswer)
gs_ls()

# load in sheet from google sheets
gonad.dat.gs <- gs_title("East Coast Morphometrics - 2015- eggs-testes data")
gonad.dat <- data.frame(get_via_csv(gonad.dat.gs))

#########

gonad.dat<-read.csv("eggs_testes_data.csv")

##separate eggs and testes data##

gonad.dat<-subset(gonad.dat, egg.diameter.1!='NA' | teste.length.1 !='NA')

eggdat<-subset(gonad.dat, egg.diameter.1!='NA' & membership !='NA' & std.length !='NA') 
eggdat<-subset(eggdat, select= -c(teste.length.1, teste.length.2,
                                 teste.width.1, teste.width.2))

testedat<-subset(gonad.dat, teste.length.1!='NA' & membership !='NA' & std.length !='NA')
testedat<-subset(testedat, select= -c(egg.diameter.1, egg.diameter.2, egg.diameter.3,
                egg.diameter.4, egg.diameter.5, egg.diameter.6, egg.diameter.7,
                egg.diameter.8, egg.diameter.9, egg.diameter.10))


##get average weight and size for eggs##

eggdat$ID <- paste(eggdat$population, eggdat$individual, sep="_")

eggdat <- transform(eggdat, avg.diameter = rowMeans(eggdat[,6:15], na.rm=TRUE))

eggdat$avg.weight<-eggdat$weight/10


#graph avg egg size and weight by body size (standard length) 

library("ggplot2")
library("ggthemes")

ggplot(eggdat, aes(x=std.length, y=avg.diameter, color=factor(membership))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T) +
  labs(
    x = "Standard Length",
    y = "Avg. egg diameter",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))

ggplot(eggdat, aes(x=std.length, y=avg.weight, color=factor(membership))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T) + 
  labs(
    x = "Standard Length",
    y = "Avg. egg weight",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


#boxplots; standard length controlled through residuals of linear regression

  #diameter
egg.diam.lm <- lm(avg.diameter~std.length, data=eggdat)
eggdat$egg.diam.resid <- residuals(egg.diam.lm)

  #Diameter by Membership
ggplot(eggdat, aes(x=membership, egg.diam.resid, colour=factor(membership))) +
  geom_boxplot() +
  labs(
    x = "Membership",
    y = "Residuals for avg. egg diameter",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


  #weight
egg.weight.lm <- lm(avg.weight~std.length, data=eggdat)
eggdat$egg.weight.resid <- residuals(egg.weight.lm)

  #weight by membership
ggplot(eggdat, aes(x=membership, egg.weight.resid, colour=factor(membership))) +
  geom_boxplot() + 
  labs(
    x = "Membership",
    y = "Residuals for avg. egg weight",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


##hypothesis testing: difference in egg size & weight by membership

library("nlme")

#size
egg.test<-lme(avg.diameter ~ std.length * membership, data=eggdat, random= ~1|population, na.action=na.omit) 
summary(egg.test) 

anova(egg.test) # type I anova                         

  #linear model fit anova 
egg.diam.anova<-aov(eggdat$avg.diameter ~ eggdat$membership)
summary(egg.diam.anova)


#weight
egg.weight.test<-lme(avg.weight ~ std.length * membership, data=eggdat, random= ~1|population, na.action=na.omit) 
summary(egg.weight.test) 

anova(egg.weight.test) # type I anova                         

  #linear model fit anova 
egg.weight.anova<-aov(eggdat$avg.weight ~ eggdat$membership)
summary(egg.weight.anova)


##get average weight, size, and variance for testes##

testedat$ID <- paste(testedat$population, testedat$individual, sep="_")

testedat <- transform(testedat, avg.size = (rowMeans(testedat[,6:7], na.rm=TRUE)) *
                             rowMeans(testedat[,8:9], na.rm=TRUE))

testedat$avg.weight<-testedat$weight/2

testedat$variance <- abs((testedat[,6]*testedat[,8]) - (testedat[,7]*testedat[,9]))


#graph avg teste size, weight, and variance by body size (standard length) 

library("ggplot2")

ggplot(testedat, aes(x=std.length, y=avg.size, color=factor(membership))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T) +
  labs(
    x = "Standard Length",
    y = "Avg. teste size (length*width)",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))

ggplot(testedat, aes(x=std.length, y=avg.weight, color=factor(membership))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T) +
  labs(
    x = "Standard Length",
    y = "Avg. teste weight",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))

ggplot(testedat, aes(x=std.length, y=variance, color=factor(membership))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T)


#boxplots: size and weight differences

  #size lm to get residuals
teste.size.lm <- lm(avg.size~std.length, data=testedat, na.action=na.exclude)
testedat$teste.size.resid <- residuals(teste.size.lm)

  #size by membership
ggplot(testedat, aes(x=membership, teste.size.resid, colour=factor(membership))) +
  geom_boxplot(na.rm=T) + 
  labs(
    x = "Membership",
    y = "Residuals for avg. teste size",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))

  #weight lm for residuals
teste.weight.lm <- lm(avg.weight~std.length, data=testedat, na.action=na.exclude)
testedat$teste.weight.resid <- residuals(teste.weight.lm)

  #weight by membership
ggplot(testedat, aes(x=membership, teste.weight.resid, colour=factor(membership))) +
  xlab("Membership") + ylab("Residuals for  avg. teste weight") +
  geom_boxplot(na.rm=T) +
  labs(
    x = "Membership",
    y = "Residuals for avg. teste weight",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


##hypothesis testing: difference in teste size & weight by membership

library("nlme")

#size differences 
teste.test<-lme(avg.size ~ std.length * membership, data=testedat, random= ~1|population, na.action=na.omit) 
summary(teste.test) 

anova(teste.test)  #type I anova                        

  #fit linear model anova
teste.size.anova<-aov(testedat$avg.size ~ testedat$membership)
summary(teste.size.anova)


#weight differences

teste.weight.test<-lme(avg.weight ~ std.length * membership, data=testedat, random= ~1|population, na.action=na.omit) 
summary(teste.weight.test) 

anova(teste.weight.test)  # type I anova                        

  #fit linear model anova
teste.weight.anova<-aov(testedat$avg.weight ~ testedat$membership)
summary(teste.weight.anova)

