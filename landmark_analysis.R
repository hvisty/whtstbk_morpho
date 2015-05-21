##landmark analysis in Geomorph##

install.packages("geomorph")
library("geomorph")
library("ggplot2")

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
ggplot(PC.scores, aes(x=PC2, y=PC3, color=landmarks$species)) + geom_point() 
#...nothin to show for shape differences. Dang. 

#MANOVA 
PC1 <- as.vector(PCA$pc.scores[,1]) #get a vector of PC1 to use as a covariate


procD.lm(GPA.2D~landmarks$population*landmarks$species, iter=99) #without PC1 bending accounted for
procD.lm(GPA.2D~landmarks$population*landmarks$species*PC1, iter=99) #if PC1 is a covariate...

#Linear discriminant function analysis 

install.packages("MASS")
library("MASS")

lda(GPA.2D, landmarks$species) #what does this mean? 
