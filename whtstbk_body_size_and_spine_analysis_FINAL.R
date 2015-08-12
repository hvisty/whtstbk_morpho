#### Final version: Body Size and Spine Length analysis in White vs common stickleback #### 

## Initial uploading 

# Libraries 
library("nlme")
library("car")
library("geomorph")
library("ggplot2")
library("candisc")
library("MASS")
library("dplyr")
library("googlesheets")
library("ggplot2")
library("ggthemes")


# Read in morphometrics data from googlesheets 

# google sheets API
library(devtools)
devtools::install_github("jennybc/googlesheets")


# list your sheets and authenticate with server (will open browswer)
gs_ls()

# load in sheet from google sheets
morphodat.gs <- gs_title("East Coast Morphometrics - 2015")
morphodat <- data.frame(gs_read_csv(morphodat.gs))




### Body Size Analysis ### 

## Cut data down to body size data 

size.dat <- morphodat[,c('population', 'individual', 'sex', 'species', 'std.length',
                          'depth.1', 'depth.2')]
size.dat$species<-as.factor(size.dat$species)


# Boxplot body depth by species
ggplot(size.dat, aes(x=species, y=depth.1, color=factor(species))) +
  geom_boxplot(na.rm=T) +
  labs(
    x = "Species",
    y = "Body Depth from first dorsal spine (cm)",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue"))


# Boxplot standard length by species
ggplot(size.dat, aes(x=species, y=std.length, color=factor(species))) +
  geom_boxplot(na.rm=T) +
  labs(
    x = "Species",
    y = "Standard Length (cm)",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue"))


## Hypothesis testing 

# Create averages by population for std. length and body depth 
length.pop <- aggregate(size.dat$std.length, by=list(size.dat$population, 
              size.dat$species),FUN='mean', na.rm=T)
colnames(length.pop) <- c("population", "species", "avg.SL")

depth.pop <- aggregate(size.dat$depth.1, by=list(size.dat$population, 
              size.dat$species),FUN="mean", na.rm=T)
colnames(depth.pop) <- c("population", "species", "avg.depth")
  
size.pop <- merge(length.pop, depth.pop, by=c("population", "species"))  


# One-way anova for average population std. length by species 
std.length.anova<-aov(size.pop$avg.SL ~ size.pop$species)
summary(std.length.anova)

# One-way anova for average population body depth by species 
depth.anova<-aov(size.pop$avg.depth ~ size.pop$species)
summary(depth.anova)





### Spine Length Analysis ### 

## Cut morphodat down to spine data 

spine.dat <- morphodat[,c('population', 'individual', 'sex', 'species', 'std.length',
                          'spine.1', 'spine.2', 'spine.3', 'spine.ventral')]
spine.dat$species<-as.factor(spine.dat$species)


# Plot for spine 1 length by standard length: 'both' populations removed. No significance 
spine.dat%>%
  filter(species != 'B') %>%
  ggplot(aes(x=std.length, y=spine.1, color=factor(species))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T) +
  labs(
    x = "Standard Length (cm)",
    y = "1st Dorsal spine length (cm)",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))

# Plot for spine 2 length by standard length: 'both' populations removed. No significance 
spine.dat%>%
  filter(species != 'B') %>%
  ggplot(aes(x=std.length, y=spine.2, color=factor(species))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T) +
  labs(
    x = "Standard Length (cm)",
    y = "2nd Dorsal spine length (cm)",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))

# Plot for pelvic spine length by standard length: 'both' populations removed. No significance 
spine.dat%>%
  filter(species != 'B') %>%
  ggplot(aes(x=std.length, y=spine.ventral, color=factor(species))) +
  geom_point(na.rm=T) +
  stat_smooth(method=lm, na.rm=T) +
  labs(
    x = "Standard Length (cm)",
    y = "Pelvic spine length (cm)",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


# Boxplots: control for standard length with residuals of linear regression 

# Spine 1 linear regression and boxplot
spine1.lm <- lm(spine.1 ~ std.length, data=spine.dat)
spine1.dat <- subset(spine.dat, spine.1 != 'NA')
spine1.dat$spine1.resid <- residuals(spine1.lm)

ggplot(spine1.dat, aes(species, spine1.resid, colour=factor(species))) + 
  geom_boxplot(na.rm=T) +
  labs(
    x = "Species",
    y = "Residuals for length of 1st dorsal spine",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue"))


# Spine 2 linear regression and boxplot
spine2.lm <- lm(spine.2 ~ std.length, data=spine.dat)
spine2.dat <- subset(spine.dat, spine.2 != 'NA')
spine2.dat$spine2.resid <- residuals(spine2.lm)

ggplot(spine2.dat, aes(species, spine2.resid, colour=factor(species))) + 
  geom_boxplot(na.rm=T) +
  labs(
    x = "Species",
    y = "Residuals for length of 2nd dorsal spine",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue"))


# Pelvic Spine linear regression and boxplot
spine.p.lm <- lm(spine.ventral ~ std.length, data=spine.dat)
spine.p.dat <- subset(spine.dat, spine.ventral != 'NA')
spine.p.dat$spine.p.resid <- residuals(spine.p.lm)

ggplot(spine.p.dat, aes(species, spine.p.resid, colour=factor(species))) + 
  geom_boxplot(na.rm=T) +
  labs(
    x = "Species",
    y = "Residuals for length of pelvic spine",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue"))



## Hypothesis testing: spine length by species 


# Mixed model for length of first spine, controlling for standard length 
first.spine<-lme(spine.1 ~ std.length * species, data=spine.dat, 
                 random= ~1|population, na.action=na.omit) 
summary(first.spine) 

anova(first.spine) # type I anova


# Mixed model for length of second spine, controlling for standard length 
second.spine<-lme(spine.2 ~ std.length * species, data=spine.dat, 
                  random= ~1|population, na.action=na.omit) 
summary(second.spine) 

anova(second.spine) # type I anova


# Mixed model for length of pelvic spine, controlling for standard length 
pelvic.spine<-lme(spine.ventral ~ std.length * species, data=spine.dat, 
                  random= ~1|population, na.action=na.omit) 
summary(pelvic.spine) 

anova(pelvic.spine) # type I anova



