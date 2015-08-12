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
#   corrections applied to remove severe outliers (see initial analysis)

landmarks <- read.csv("corrected_landmark_data.csv")

# Reclassify SF and SR as 'W' populations
for (i in 1:nrow(landmarks)){
  if(landmarks$population[i] =="SF2014" | landmarks$population[i] == "SR2014") {
    landmarks$species[i] <- "W"
  }
}

select<- dplyr::select #allows dplyr to use select when MASS package is loaded



## GEOMORPH Analysis ## 

# Create a 3D array with just the x/y coordinates 
landmarks.data <- landmarks %>%
  select(matches("^[x,y]{1}."))

landmark.array <- arrayspecs(landmarks.data, 19, 2) #3D array

# Procrustes analysis
GPA.landmarks <- gpagen(landmark.array) 

GPA.2D <- two.d.array(GPA.landmarks$coords) #2D Data frame of procrustes coordinates


## Principle Component Analysis ##

PCA<- plotTangentSpace(GPA.landmarks$coords, groups= factor(landmarks$population), verbose=T)
PCA$pc.scores #gives the PCA scores

# 'Unbending' ; just show me PC2 vs PC3
plotTangentSpace(GPA.landmarks$coords, groups=factor(landmarks$species), axis1=2, axis2=3)

# Make a prettier plot 
PC.scores<-as.data.frame(PCA$pc.scores)
ggplot(PC.scores, aes(x=PC2, y=PC3, color=landmarks$species)) +
  geom_point(size=3) 

ggplot(PC.scores, aes(x=PC3, y=PC4, color=landmarks$species)) +
  geom_point(size=3) 
#...nothin to show for shape differences. Dang. 


## MANOVA 

PC1 <- as.vector(PCA$pc.scores[,1]) #get a vector of PC1 to use as a covariate

#without PC1 bending accounted for
procD.lm(GPA.2D~landmarks$population*landmarks$species, iter=99)

#if PC1 is a covariate...
procD.lm(GPA.2D~landmarks$population+landmarks$species*PC1, iter=99)


## Linear discriminant function analysis ##

lda.species <- lda(GPA.2D, landmarks$species) 

# the percent of variance explained by the LD funcitons
prop.lda <- lda.species$svd^2/sum(lda.species$svd^2) 
prop.lda <- round(prop.lda*100)

# project the original values in lda space
plda <- predict(object = lda.species,
                newdata = GPA.2D)

# data frame of projected data with species names
lda.project <- data.frame(species = landmarks$species, 
                          population = landmarks$population, 
                          individual = landmarks$individual,
                          ID = landmarks$ID,
                          ld1 = plda$x[,1], 
                          ld2 = plda$x[,2])

# Plot LDA without "Both" locations
lda.project %>%
  filter(species != "B") %>%
  ggplot(aes(color = species, x = ld1,y = ld2))+
  geom_point(size = 3) +
  labs(x = paste0("LD1 (", prop.lda[1], "%)"),
       y = paste0("LD2 (", prop.lda[2], "%)")) +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


## Cross Validation of LDA ##

# Method 1: take half of the values and create the LDA, 
#   then plot the rest to check model 

# LDA with CV function (automatic cross validation)
GPA.CW <- GPA
GPA.CW$species <- landmarks$species
GPA.CW$centroid = GPA.landmarks$Csize

GPA.CW<- GPA.CW %>%
  filter(species != "B") #remove "Both" category

landmarks.CW<- landmarks %>%
  filter(species != "B")

GPA.CW <- GPA.CW[-c(39,40)]

lda.species.CV <- lda(GPA.CW, landmarks.CW$species, CV=T)

# percent correct (accuracy of species prediction)
ct <- table(landmarks.CW$species, lda.species.CV$class)
diag(prop.table(ct, 1))

# total percent correct
sum(diag(prop.table(ct)))


# LDA validation Method 2: sample a random half of individuals
#   with equal representation of W and C 

# build a data frame with ID info and procrustes x/y
GPA.50<- data.frame(species = landmarks$species, 
                    population = landmarks$population, 
                    individual = landmarks$individual,
                    ID = landmarks$ID,
                    centroid = GPA.landmarks$Csize)

GPA.50<-merge(GPA.50, GPA, by.x="ID")

GPA.parse<- GPA.50 %>% 
  filter(species != "B") %>%
  group_by(species) 

# sample a random half of the data with equal W and C numbers
GPA.parse <- sample_n(GPA.parse, 403)
GPA.parse1<- GPA.parse[c(1:201, 403:604),]
GPA.parse2 <- GPA.parse[c(202:402,605:806),]

GPA.1 <- GPA.parse1[-c(1:4)]
GPA.2 <- GPA.parse2[-c(1:4)]


# run LDA for each 1/2 sample
lda.CV1 <- lda(GPA.1, as.character(GPA.parse1$species))
lda.CV2 <- lda(GPA.2, as.character(GPA.parse2$species))


#subtract scaling between the two halves and look at differences 
lda.CV <- as.vector(lda.CV1$scaling - lda.CV2$scaling)
plot(lda.CV) #weird LDA units 

# look at the plda and overplotting for each half

# project in lda space
plda.CV1 <- predict(object = lda.CV1,
                    newdata = GPA.1)

plda.CV2 <- predict(object = lda.CV2,
                    newdata = GPA.2)

# data frames of projected data with species names
lda.project.CV1 <- data.frame(species = GPA.parse1$species, 
                              population = GPA.parse1$population, 
                              individual = GPA.parse1$individual,
                              ID = GPA.parse1$ID,
                              ld1 = plda.CV1$x)

lda.project.CV2 <- data.frame(species = GPA.parse2$species, 
                              population = GPA.parse2$population, 
                              individual = GPA.parse2$individual,
                              ID = GPA.parse2$ID,
                              ld1 = plda.CV2$x)

# plot both halves on one plot... looks like they match 
ggplot(lda.project.CV1,aes(fill = species, x = LD1, alpha=0.2)) +
  geom_histogram() +
  geom_histogram(aes(x = lda.project.CV2$LD1, fill = species)) +
  theme_classic() +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))


## END cross validation ##


### Genetics VS Morphology ###

# extract centroids
landmark.centroids <- data.frame(population = landmarks$population, 
                                 id = landmarks$ID, 
                                 species = landmarks$species, 
                                 centroid = GPA.landmarks$Csize,
                                 sex= landmarks$sex)

# Plot centroid sizes by population 
landmark.centroids %>%
  ggplot(aes(x = population, y = centroid, color=species)) +
  geom_point()
  facet_grid(~ species, scale = "free", space = "free")


#read in genetic cluster data

cluster.dat <- read.table(file=list.files(pattern="structure",full.names = TRUE), stringsAsFactors = FALSE)
cluster.dat.reduced <- cluster.dat %>%
  filter(k.value.run == 2) %>%
  select(id,X1,X2,membership)

landmark.centroids.rename <- landmark.centroids
landmark.centroids.rename$id <- landmark.centroids$id %>% 
  gsub("2014","",.) %>% gsub("_","",.)

cluster.morph <- left_join(cluster.dat.reduced, landmark.centroids.rename)

# Plot centroid size by population again
cluster.morph  %>%
  filter(!is.na(population)) %>%
  ggplot(aes(x = population, y = centroid, color=factor(membership))) +
  geom_point(size=3) +
  facet_grid(~ species, scale = "free", space = "free")

