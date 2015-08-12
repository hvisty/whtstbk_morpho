#### Final Version: Color analysis for White versus Common Stickleback ####

## Initial uploading 

# Libraries 

library("dplyr") 
library("ggplot2")
library("ggthemes")
library("grDevices")


# RGB data frame from Github whtstbk_morpho (Raw histograms transformed into data 
# frame in initial version of analysis)

colordat <- read.csv("whtstbk_color_data.csv")



## Transformations: get separate R, G, B data frames & change histogram frequency 
##  counts to reflect 0-255 RGB scale 

#get separate R G B matrices

indx <- gsub("_.*", "", names(colordat))

list2env(
  setNames(
    lapply(split(colnames(colordat), indx), function(x) colordat[x]),
    paste('colordat', sort(unique(indx)), sep="_")), 
  envir=.GlobalEnv)

#transform R G B raw pixel values to reflect 0-255 

R.dat<-colordat_R

for(i in 1:ncol(R.dat)){
  R.dat[,i] <- R.dat[,i]*(i-1)
}

sum<-rowSums(colordat_R)
sum2<-rowSums(R.dat)
R.dat$mean<-sum2/sum


G.dat<-colordat_G

for(i in 1:ncol(G.dat)){
  G.dat[,i] <- G.dat[,i]*(i-1)
}

sum<-rowSums(colordat_G)
sum2<-rowSums(G.dat)
G.dat$mean<-sum2/sum


B.dat<-colordat_B

for(i in 1:ncol(B.dat)){
  B.dat[,i] <- B.dat[,i]*(i-1)
}

sum<-rowSums(colordat_B)
sum2<-rowSums(B.dat)
B.dat$mean<-sum2/sum



## Analysis with conversion to HSV color scale. RGB analysis on initial version. 

#convert RGB to HSV data frame

HSV.dat<-t(rgb2hsv(r=R.dat$mean, g=G.dat$mean, b=B.dat$mean))
HSV.dat<-as.data.frame(HSV.dat)

HSV.dat$species <- colordat$species
HSV.dat$population <- colordat$population
HSV.dat$ID<-colordat$ID


#get population means for Hue, Saturation, and Value 

H.pop.means<-aggregate(HSV.dat$h, by=c(colordat_population, colordat_species), FUN="mean")

S.pop.means<-aggregate(HSV.dat$s, by=c(colordat_population, colordat_species), FUN="mean")

V.pop.means<-aggregate(HSV.dat$v, by=c(colordat_population, colordat_species), FUN="mean")


# Plot population means for Hue: no pattern
ggplot(H.pop.means, aes(x= factor(population), y= H.pop.means$x, color=factor(species))) +
  geom_point(aes(size=3)) +
  labs(
    x = "Population",
    y = "Avg. Hue ",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue")) +
  theme(axis.text.x = element_text(angle=90, vjust=1))


# Plot population means for Saturation: no pattern 
ggplot(S.pop.means, aes(x= factor(population), y= S.pop.means$x, color=factor(species))) +
  geom_point(aes(size=3)) +
  labs(
    x = "Population",
    y = "Avg. Saturation ",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue")) +
  theme(axis.text.x = element_text(angle=90, vjust=1))


# Plot population means for Value (Brightness)
ggplot(V.pop.means, aes(x= factor(population), y= V.pop.means$x, color=factor(species))) +
  geom_point(aes(size=3)) +
  labs(
    x = "Population",
    y = "Avg. Value (Brightness) ",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue")) +
  theme(axis.text.x = element_text(angle=90, vjust=1))


# Plot Value (Brightness) by Species Boxplot
ggplot(HSV.dat, aes(x= species, y= v, color=factor(species))) +
  geom_boxplot(na.rm=T) +
  labs(
    x = "Species",
    y = "Value (Brightness)",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue"))


# Compare brightness value to body size (as proxy for species ID)

HSV.dat$std.length<-colordat$std.length

HSV.dat%>%
  filter(species!='B')%>%
  ggplot(aes(x=std.length, y=v, color=factor(species))) + 
  geom_point(size=3) +
  labs(
    x = "Standard Length (cm)",
    y = "Value (Brightness)",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


# Compare brightness to species in males with genetic membership data

# Read in membership IDs and join with HSV data 
gene.dat<-read.csv("eggs_testes_data.csv")

gene.dat<-subset(gene.dat, sequenced.!='NA')

gene.dat$ID<-as.character(paste(gene.dat$population, gene.dat$individual, sep="_"))
gene.dat<-gene.dat[,c('ID', 'membership', 'sex')]

genetic.HSV.dat<-left_join(gene.dat, HSV.dat, by='ID')


# Plot brightness by membership for males 
male.HSV.dat<-filter(genetic.HSV.dat, sex!='M')

ggplot(male.HSV.dat,aes(x=std.length, y=v, color=factor(membership))) + 
  geom_point(na.rm=T) +
  #stat_smooth(method=lm, na.rm=T) +
  labs(
    x = "Standard Length (cm)",
    y = "Value (Brightness) for males",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


# Boxplot brightness by membership for males
ggplot(male.HSV.dat, aes(x=membership, y=v, color=factor(membership))) + 
  geom_boxplot() +
  labs(
    x = "Membership",
    y = "Value (Brightness) for males",
    color = "Membership") +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


## Hypothesis testing: Value (Brightness) is significantly different by species

V.anova<-aov(HSV.dat$v ~ HSV.dat$species)
summary(V.anova)
TukeyHSD(V.anova) #all 3 categories are significantly different 


# anova using only genetic data (common and white membership) also shows significance
V.genetic.anova<-aov(genetic.HSV.dat$v ~ genetic.HSV.dat$membership)
summary(V.genetic.anova)
