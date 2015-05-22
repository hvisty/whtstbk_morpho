##landmark analysis in Geomorph##

install.packages("geomorph")
install.packages("candisc")
install.packages("MASS")
install.packages("candisc")

library("geomorph")
library("ggplot2")
library("candisc") # will use this later
library("MASS")
library("dplyr")

landmarks <- read.csv("Corrected Landmark Data.csv")

#Used to find outliers, plot any x/y points with indiv. ID#
#####
landmarks$ID <- paste(landmarks$population, landmarks$individual, sep="_")
ggplot(landmarks, aes(x=x.4, y=y.4, label=ID, color=population)) +geom_point() +geom_text()

#####

#Make a 3D array; only x/y cordinates, outliers removed...

landmarks <- landmarks[-c(10,125,271,253,327,360,448,544,630,641,665,877,882,887,898:900,947,959:961,997, 1061, 1172),]
landmarks <-landmarks[complete.cases(landmarks[,"x.1"]),]
landmarks.data <- landmarks[, c(5:42)]

landmark.array <- arrayspecs(landmarks.data, 19, 2) #3D array

GPA.landmarks<-gpagen(landmark.array) #Procrustes analysis

GPA.2D<-two.d.array(GPA.landmarks$coords) #2D Data frame of procrustes coordinates

#principle component stuff. 
PCA<- plotTangentSpace(GPA.landmarks$coords, groups= landmarks$population, verbose=T)
PCA$pc.scores #gives the PCA scores

#unbending; just show me PC2 vs PC3
plotTangentSpace(GPA.landmarks$coords, groups=landmarks$species, axis1=2, axis2=3)

#make a prettier plot! 
PC.scores<-as.data.frame(PCA$pc.scores)
ggplot(PC.scores, aes(x=PC2, y=PC3, color=landmarks$species)) + geom_point(size=3) 

ggplot(PC.scores, aes(x=PC5, y=PC6, color=landmarks$species)) + geom_point(size=3) 
#...nothin to show for shape differences. Dang. 

#MANOVA 
PC1 <- as.vector(PCA$pc.scores[,1]) #get a vector of PC1 to use as a covariate

procD.lm(GPA.2D~landmarks$population*landmarks$species, iter=99) #without PC1 bending accounted for
procD.lm(GPA.2D~landmarks$population+landmarks$species*PC1, iter=99) #if PC1 is a covariate...

#Linear discriminant function analysis 

lda.species <- lda(GPA.2D, landmarks$species) #what does this mean? 

# the default plot, looks meh
plot(lda.species)

# project the original values in lda space
plda <- predict(object = lda.species,
                newdata = GPA.2D)

# the percent of variance explained by the LD funcitons
prop.lda <- lda.species$svd^2/sum(lda.species$svd^2) 
prop.lda <- round(prop.lda*100)

# data frame of projected data with species names
lda.project <- data.frame(species = landmarks$species, 
                          population = landmarks$population, 
                          ld1 = plda$x[,1], 
                          ld2 = plda$x[,2])

# plot with ggplot
lda.project %>%
  ggplot(aes(color = species, x = ld1,y = ld2))+
    geom_point(size = 3) +
    labs(x = paste0("LD1 (", prop.lda[1], "%)"),
         y = paste0("LD2 (", prop.lda[2], "%)"))

# no "Both" locations
# ooo pretty
lda.project %>%
  filter(species != "B") %>%
  ggplot(aes(color = species, x = ld1,y = ld2))+
    geom_point(size = 3) +
    labs(x = paste0("LD1 (", prop.lda[1], "%)"),
         y = paste0("LD2 (", prop.lda[2], "%)"))
