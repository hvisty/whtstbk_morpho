#### Final Verison: Egg and Teste analysis for White versus Common stickleback ####

## Initial uploading 

# Libraries 

library("ggplot2")
library("ggthemes")
library("nlme")
library("dplyr")


# Raw data from GitHub whtstbk_morpho
gonad.dat<-read.csv("eggs_testes_data.csv")


# Separate eggs and testes data: remove NAs, get individual data frames 

gonad.dat<-subset(gonad.dat, egg.diameter.1!='NA' | teste.length.1 !='NA')

eggdat<-subset(gonad.dat, egg.diameter.1!='NA' & membership !='NA' & std.length !='NA') 
eggdat<-subset(eggdat, select= -c(teste.length.1, teste.length.2,
                                  teste.width.1, teste.width.2))

testedat<-subset(gonad.dat, teste.length.1!='NA' & membership !='NA' & std.length !='NA')
testedat<-subset(testedat, select= -c(egg.diameter.1, egg.diameter.2, egg.diameter.3,
                                      egg.diameter.4, egg.diameter.5, egg.diameter.6, 
                                      egg.diameter.7,egg.diameter.8, egg.diameter.9, 
                                      egg.diameter.10, egg.number))



## Egg Analysis ##


# Get averages for weight and size by individual 

eggdat$ID <- paste(eggdat$population, eggdat$individual, sep="_")

eggdat <- transform(eggdat, avg.diameter = rowMeans(eggdat[,6:15], na.rm=TRUE))

eggdat$avg.weight<-eggdat$weight/10


# Average egg diameter by Std. Length: no significant difference 
ggplot(eggdat, aes(x=std.length, y=avg.diameter, color=factor(membership))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T) +
  labs(
    x = "Standard Length (cm)",
    y = "Avg. egg diameter (mm)",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


# Average egg weight by Std. Length: no significant difference 
ggplot(eggdat, aes(x=std.length, y=avg.weight, color=factor(membership))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T) + 
  labs(
    x = "Standard Length (cm)",
    y = "Avg. egg weight (mg)",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))

# Egg number by Std. Length: no significant difference 
ggplot(eggdat, aes(x=std.length, y=egg.number, color=factor(membership))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T) + 
  labs(
    x = "Standard Length (cm)",
    y = "Number of eggs",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))

# Comparison of egg size versus egg number 
ggplot(eggdat, aes(x=avg.diameter, y=egg.number, color=factor(membership))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T) + 
  labs(
    x = "Average egg diameter (mm)",
    y = "Number of eggs",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


# Boxplots by membership; Std. Length controlled using residuals of linear regression

# Egg Diameter linear regression and boxplot 
egg.diam.lm <- lm(avg.diameter~std.length, data=eggdat)
eggdat$egg.diam.resid <- residuals(egg.diam.lm)

ggplot(eggdat, aes(x=membership, egg.diam.resid, colour=factor(membership))) +
  geom_boxplot() +
  labs(
    x = "Membership",
    y = "Residuals for avg. egg diameter",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


# Egg Weight linear regression and boxplot 
egg.weight.lm <- lm(avg.weight~std.length, data=eggdat)
eggdat$egg.weight.resid <- residuals(egg.weight.lm)

ggplot(eggdat, aes(x=membership, egg.weight.resid, colour=factor(membership))) +
  geom_boxplot() + 
  labs(
    x = "Membership",
    y = "Residuals for avg. egg weight",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


# Egg Number linear regression and boxplot 
egg.num.lm <- lm(egg.number~std.length, data=eggdat)
eggdat.num<-subset(eggdat, egg.number !='NA') 
eggdat.num$egg.number.residuals<- residuals(egg.num.lm)

ggplot(eggdat.num, aes(x=membership, y=egg.number.residuals, colour=factor(membership))) +
  geom_boxplot() +
  labs(
    x = "Membership",
    y = "Residuals for number of eggs",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


# Hypothesis testing for eggs: difference in weight/size by membership
  # neither show significance in any area except for the 'intercept'

# mixed model for egg diameter
egg.test<-lme(avg.diameter ~ std.length * membership, data=eggdat, 
              random= ~1|population, na.action=na.omit) 
summary(egg.test) 

anova(egg.test) # type I anova                         


# mixed model for egg weight 
egg.weight.test<-lme(avg.weight ~ std.length * membership, data=eggdat, 
                     random= ~1|population, na.action=na.omit) 
summary(egg.weight.test) 

anova(egg.weight.test) # type I anova                         


# mixed model for egg number
egg.number.test<-lme(egg.number ~ std.length * membership, data=eggdat, 
                     random= ~1|population, na.action=na.omit) 
summary(egg.number.test) 

anova(egg.number.test) # type I anova 


## Testes Analysis ##


# Get average teste size (length*width) and weight (/2) by individual. 

testedat$ID <- paste(testedat$population, testedat$individual, sep="_")

testedat <- transform(testedat, avg.size = (rowMeans(testedat[,6:7], na.rm=TRUE)) *
                        rowMeans(testedat[,8:9], na.rm=TRUE))

testedat$avg.weight<-testedat$weight/2


# Plot average size (length*width) by Std. Length: no significant difference
ggplot(testedat, aes(x=std.length, y=avg.size, color=factor(membership))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T) +
  labs(
    x = "Standard Length (cm)",
    y = "Avg. teste size (length*width) (mm)",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))

# Plot average weight by Std. Length: significant (outliers removed)
testedat %>%
  filter(avg.weight<1) %>%
  ggplot(aes(x=std.length, y=avg.weight, color=factor(membership))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T) +
  labs(
    x = "Standard Length (cm)",
    y = "Avg. teste weight (mg)",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))

# Plot teste size vs weight 
ggplot(testedat, aes(x=avg.weight, y=avg.size, color=factor(membership))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T) +
  labs(
    x = "Avg. teste weight (mg)",
    y = "Avg. teste size (length*width) (mm)",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


# Boxplots by membership; Std. Length controlled using residuals of linear regression

# Teste size linear regression and boxplot 
teste.size.lm <- lm(avg.size~std.length, data=testedat, na.action=na.exclude)
testedat$teste.size.resid <- residuals(teste.size.lm)

ggplot(testedat, aes(x=membership, teste.size.resid, colour=factor(membership))) +
  geom_boxplot(na.rm=T) + 
  labs(
    x = "Membership",
    y = "Residuals for avg. teste size",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))

# Teste weight linear regression and boxplot 
teste.weight.lm <- lm(avg.weight~std.length, data=testedat, na.action=na.exclude)
testedat$teste.weight.resid <- residuals(teste.weight.lm)

ggplot(testedat, aes(x=membership, teste.weight.resid, colour=factor(membership))) +
  geom_boxplot(na.rm=T) +
  labs(
    x = "Membership",
    y = "Residuals for avg. teste weight",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


# Hypothesis testing for testes: difference in weight/size by membership
  # weight shows significance in the interaction between weight&membership

# mixed model for teste size 
teste.test<-lme(avg.size ~ std.length * membership, data=testedat, 
                random= ~1|population, na.action=na.omit) 
summary(teste.test) 

anova(teste.test)  #type I anova                        

#mixed model for tese weight
teste.weight.test<-lme(avg.weight ~ std.length * membership, data=testedat, 
                       random= ~1|population, na.action=na.omit) 
summary(teste.weight.test) 

anova(teste.weight.test)  # type I anova                        
