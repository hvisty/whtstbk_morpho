### Color Analysis ###


## Get RGB data into useable form ##

library(dplyr)

########Preliminary, not needed to repeat after initial csv write########
# list files
color.data.files <- list.files("color data", full.names=T)

color.data.frame <- data.frame(row.names=color.data.files)


# build a data frame with all individuals and all RGB data 

for (i in 1:length(color.data.files)){
  
  file.R <- read.csv(file= color.data.files[i], header=T,
                     colClasses = c("NULL", NA, "NULL", "NULL"))
  file.G <- read.csv(file= color.data.files[i], header=T,
                     colClasses = c("NULL", "NULL", NA, "NULL"))
  file.B <- read.csv(file= color.data.files[i], header=T,
                     colClasses = c("NULL", "NULL", "NULL", NA))
  color.columns <- c(t(file.R), t(file.G), t(file.B))
  color.data.frame <- rbind(color.data.frame, color.columns)
  
}

#name the rows and columns, formatting, blah blah blah 

color.data.frame$ID <- list.files("color data")
color.data.frame$ID <- gsub("Histogram of", "", color.data.frame$ID)
color.data.frame$ID <- gsub(".csv", "", color.data.frame$ID)

color.data.frame$population <- 
  as.character(lapply(strsplit(color.data.frame$ID, split="_"), "[", 1))
color.data.frame$individual <- 
  as.character(lapply(strsplit(color.data.frame$ID, split="_"), "[", 2))


nameR <- paste("R", 0:255, sep="_")
colnames(color.data.frame)[1:256] <- nameR

nameG <- paste("G", 0:255, sep="_")
colnames(color.data.frame)[257:512] <- nameG

nameB <- paste("B", 0:255, sep="_")
colnames(color.data.frame)[513:768] <- nameB


#master.data.frame<-read.csv("east_coast_morphometrics_2015.csv")
#master.data.frame$ID<- paste(master.data.frame$population, master.data.frame$individual, sep="_")
#master.data.frame$species<- as.character(master.data.frame$species)
#master.data.frame$population<- as.character(master.data.frame$population)
#master.data.frame$individual<- as.character(master.data.frame$individual)

color.data.frame$species <- c(rep('C',17), rep('B', 44), rep('C', 142), rep('B', 62), rep('C', 36), rep('B', 188), rep('C',111),
             rep('W', 46), rep('C',30), rep('B',91), rep('W',179), rep('C', 67), rep('B', 196), rep('W', 22))


write.csv(color.data.frame, file="whtstbk_color_data.csv", row.names=FALSE)

######## end of preliminary ########



colordat <- read.csv("whtstbk_color_data.csv")


## BEGIN ANALYSIS ##

library(ggplot2)
library(ggthemes)

#get separate R G B matrices to work with 

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


#looking at mean R G and B by population

# R average and boxplot distribution

R.dat$species <- colordat$species
R.dat$population <- colordat$population
Rpop.means<-aggregate(R.dat$mean, by=c(colordat_population, colordat_species), FUN='mean')

ggplot(Rpop.means, aes(x= population, y= x, color=factor(species))) +
  geom_point(aes(size=3)) +
  labs(
    x = "Population",
    y = "Avg. 'R' color value",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue")) +
  theme(axis.text.x = element_text(angle=90, vjust=1))

ggplot(R.dat, aes(x= factor(R.dat$population), y= R.dat$mean, color=factor(species))) + 
  geom_boxplot() +
  labs(
    x = "Population",
    y = "'R' color value",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue")) +
  theme(axis.text.x = element_text(angle=90, vjust=1))

# G average and boxplot distribution by population


G.dat$species <- colordat$species
G.dat$population <- colordat$population
Gpop.means<-aggregate(G.dat$mean, by=c(colordat_population, colordat_species), FUN="mean")

ggplot(Gpop.means, aes(x= population, y= x, color=factor(species))) +
  geom_point(aes(size=3)) +
  labs(
    x = "Population",
    y = "Avg. 'G' color value",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue")) +
  theme(axis.text.x = element_text(angle=90, vjust=1))

ggplot(G.dat, aes(x= factor(G.dat$population), y= G.dat$mean, color=factor(species))) + 
  geom_boxplot() +
  labs(
    x = "Population",
    y = "'G' color value",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue")) +
  theme(axis.text.x = element_text(angle=90, vjust=1))

# B average and boxplot distribution by population

B.dat$species <- colordat$species
B.dat$population <- colordat$population
Bpop.means<-aggregate(B.dat$mean, by=c(colordat_population, colordat_species), FUN="mean")

ggplot(Bpop.means, aes(x= population, y= x, color=factor(species))) +
  geom_point(aes(size=3)) +
  labs(
    x = "Population",
    y = "Avg. 'B' color value",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue")) +
  theme(axis.text.x = element_text(angle=90, vjust=1))

ggplot(B.dat, aes(x= factor(B.dat$population), y= B.dat$mean, color=factor(species))) + 
  geom_boxplot() +
  labs(
    x = "Population",
    y = "'B' color value",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue")) +
  theme(axis.text.x = element_text(angle=90, vjust=1))

#averages for all 3 overlaid
ggplot(Rpop.means, aes(x= population, y= x, color=factor(species))) +
  geom_point(aes(size=3)) +
  geom_point(data= Gpop.means, aes(x= population, y= x, color=factor(species), size=3)) +
  geom_point(data= Bpop.means, aes(x= population, y= x, color=factor(species), size=3))

#boxplots by species 

ggplot(R.dat, aes(x= species, y= mean, color=factor(species))) +
  geom_boxplot() +
  labs(
    x = "Species",
    y = "'R' color value",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue")) 


ggplot(G.dat, aes(x= species, y= mean, color=factor(species))) +
  geom_boxplot() +
  labs(
    x = "Species",
    y = "'G' color value",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue"))


ggplot(B.dat, aes(x= species, y= mean, color=factor(species))) +
  geom_boxplot() +
  labs(
    x = "Species",
    y = "'B' color value",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue"))



#overall metric for brightness? sum bins/3

means<-data.frame(R.dat$mean, G.dat$mean, B.dat$mean)
brightness<-as.data.frame((rowSums(means))/3)
colnames(brightness) <- "brightness"
brightness$species<-colordat$species
brightness$population<-colordat$population
brightness$ID<-colordat$ID

ggplot(brightness, aes(x=species, y= brightness, color=factor(species))) +
  geom_boxplot() +
  labs(
    x = "Species",
    y = "'Brightness' color value",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue"))


brightness.avg<-aggregate(brightness$brightness, 
                          by=c(colordat_population, colordat_species), FUN='mean')

ggplot(brightness.avg, aes(x= population, y= x, color=factor(species))) +
  geom_point(aes(size=3)) +
  labs(
    x = "Population",
    y = "Avg. 'Brightness' value",
    color = "Species") +
  theme_classic() +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue")) +
  theme(axis.text.x = element_text(angle=90, vjust=1))

##Hypothesis Testing##

#R, G, and B; mean population averages by species

R.anova<-aov(Rpop.means$x ~ Rpop.means$species)
summary(R.anova) #F=5.904, p>.00966
TukeyHSD(R.anova) #W-B has p>.008 but neither are different from C

G.anova<-aov(Gpop.means$x ~ Gpop.means$species)
summary(G.anova) #F=6.353, p>.00731
TukeyHSD(G.anova) #W-B has p>.006 but neither are different from C

B.anova<-aov(Bpop.means$x ~ Bpop.means$species)
summary(B.anova) #F=2.29, p>.127
#no significance 

Brightness.anova<-aov(brightness.avg$x ~ brightness.avg$species)
summary(Brightness.anova) #F=6.202, p>.00802
TukeyHSD(Brightness.anova) #W-B has p>.006 but neither are different from C 


