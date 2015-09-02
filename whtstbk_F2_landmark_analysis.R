#### Final Version: Landmark morphometrics analysis for White and Common Stickleback ####


## Initial uploading 

# Libraries 
library("nlme")
library("car")
library("geomorph")
library("ggplot2")
library("candisc")
library("MASS")
library("dplyr")
library("ggplot2")
library("ggthemes")


# Read in corrected landmark file from Github/ whtstbk_morpho 
#   corrections applied to remove severe outliers (see postscript for removal code)

F2.data <- read.csv("Corrected F2 landmarks.csv")
# for analysis with outliers in place, read in "F2 white stickleback morphometrics.csv"


select<- dplyr::select #allows dplyr to use select when MASS package is loaded



## GEOMORPH Analysis ## 

# Create a 3D array with just the x/y coordinates 
landmarks.data <- F2.data %>%
  select(matches("^[x,y]{1}."))

landmark.array <- arrayspecs(landmarks.data, 19, 2) #3D array

# Procrustes analysis
GPA.landmarks <- gpagen(landmark.array) 

GPA.2D <- two.d.array(GPA.landmarks$coords) #2D Data frame of procrustes coordinates



## Principle Component Analysis ##

PCA<- plotTangentSpace(GPA.landmarks$coords, groups= factor(landmarks.data$family), verbose=T)
PCA$pc.scores #gives the PCA scores

# 'Unbending' ; just show me PC2 vs PC3
plotTangentSpace(GPA.landmarks$coords, groups=factor(landmarks.data$family), axis1=3, axis2=4)

# Make a prettier plot 
PC.scores<-as.data.frame(PCA$pc.scores)
ggplot(PC.scores, aes(x=PC2, y=PC3, color=factor(F2.data$family))) +
  geom_point(size=5) +
  labs( x= "PC2",
        y= "PC3", 
        color= "Family") +
  theme_classic()




## MANOVA 

PC1 <- as.vector(PCA$pc.scores[,1]) #get a vector of PC1 to use as a covariate

#without PC1 bending accounted for
procD.lm(GPA.2D~F2.data$family, iter=99)

#if PC1 is a covariate...
procD.lm(GPA.2D~F2.data$family*PC1, iter=99)




## Linear discriminant function analysis ##

lda.family <- lda(GPA.2D, F2.data$family) 

# the percent of variance explained by the LD funcitons
prop.lda <- lda.family$svd^2/sum(lda.family$svd^2) 
prop.lda <- round(prop.lda*100)

# project the original values in lda space
plda <- predict(object = lda.family,
                newdata = GPA.2D)

# data frame of projected data with species names
lda.project <- data.frame(family = F2.data$family,  
                          individual = F2.data$individual,
                          ID = F2.data$ID,
                          ld1 = plda$x[,1], 
                          ld2 = plda$x[,2])

# Plot LDA 
ggplot(lda.project, aes(color = factor(family), x = ld1,y = ld2))+
  geom_point(size = 5) +
  labs(x = paste0("LD1 (", prop.lda[1], "%)"),
       y = paste0("LD2 (", prop.lda[2], "%)")) +
  theme_classic()

# Do LDA categories match color or size differentiation? 

lda.project$melanocytes <- F2.data$melanocytes

#   ld1 versus number of melanocytes
ggplot(lda.project, aes(x=ld1, y=melanocytes, color=factor(family))) +
  geom_point(size=5) +
  theme_classic()

#   ld2 versus number of melanocytes 
ggplot(lda.project, aes(x=ld2, y=melanocytes, color=factor(family))) +
  geom_point(size=5) +
  theme_classic()



# Centroid size by family 

# extract centroids
landmark.centroids <- data.frame(family = F2.data$family, 
                                 id = F2.data$ID,  
                                 centroid = GPA.landmarks$Csize)

# Plot centroid sizes by population 

ggplot(landmark.centroids, aes(x = family, y = centroid, color=factor(family))) +
  geom_point(size=5) + 
  labs( x= "Family",
        y="Centroid size",
        color= "Family") +
  theme_classic()



#### Post-Script ####

# Finding and removing OUTLIERS from PCA/GPA data

#PCA<- plotTangentSpace(GPA.landmarks$coords, groups= factor(landmarks.data$family), verbose=T)
#PC.scores <- as.data.frame(PCA$pc.scores) #gives the PCA scores
#PC.scores$ID <- paste(F2.data$family, F2.data$individual, sep="-")


#determined outliers from PC plot....
#PC.scores <- filter(PC.scores, PC1> -0.10 & PC1< 0.05 & PC2> -0.1 & PC2< 0.1)
#ggplot(PC.scores, aes(x=PC1, y=PC2, label=ID)) +geom_point() +geom_text() #filter works

#create updated landmarks data frame, outliers removed 
#F2.data$ID<- gsub("F2-", "", F2.data$ID)

#landmarks <- left_join(PC.scores, F2.data, by=c("ID"))
#landmarks <- subset(landmarks, select=-c(PC1:PC38, X, X.1))


#write.csv(landmarks, "Corrected F2 landmarks.csv")
####