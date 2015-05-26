##landmark analysis in Geomorph##

#### LIBRARIES ####

# required libraries
install.packages("geomorph")
install.packages("candisc")
install.packages("MASS")
install.packages("candisc")
install.packages("devtools")

# google sheets API
library(devtools)
devtools::install_github("jennybc/googlesheets")

# load libraries
library("geomorph")
library("ggplot2")
library("candisc")
library("MASS")
library("dplyr")
library("googlesheets")

#### LIBRARIES ####



#### LOAD RAW DATA ####

# google sheets method

## list your sheets and authenticate with server (will open browswer)
gs_ls()

## load in sheet from google sheets
morphodat.gs <- gs_title("East Coast Morphometrics - 2015")
morphodat <- data.frame(get_via_csv(landmarks.gs))


## write csv to file (for working offline)
# write.csv(morphodat, file="east_coast_morphometrics_2015.csv")
# morphodat <- read.csv("east_coast_morphometrics_2015.csv")

# build 'corrected' landmark file
# pop, id, sex, species, x.1,y.1,x.2,y.2, etc.

## select landmark columns, remove nas
landmarks.full <- morphodat %>%
  select(population, individual, sex, species, correction, contains("lndmrk")) %>%
  filter(!is.na(lndmrk.x1))

## extract ids
landmark.ids <- landmarks.full %>%
  select(population, individual, sex, species)

## fix names
names(landmarks.full) <- gsub("lndmrk.","",names(landmarks.full)) %>% gsub ("^x", "x.", .) %>% gsub ("^y", "y.", .)

## strip whitespace (...) and apply correction
landmarks.matrix <- landmarks.full %>%
  select(matches("^[x,y]{1}.")) %>%
  as.list %>%
  lapply(function(x)gsub(" ","",x)) %>%
  lapply(as.numeric) %>%
  as.data.frame %>%
  as.matrix

## apply correction
landmarks.matrix <- data.frame(sweep(landmarks.matrix, MARGIN=1, landmarks.full$correction, `/`))

## rearrange so x.1,y.1, not x.1,x.2
landmarks.order <- landmarks.matrix  %>%
  names %>%
  gsub ("^x.", "", .) %>%
  gsub ("^y.", "", .) %>%
  as.numeric %>%
  order

landmarks.raw <- cbind(landmarks.ids, landmarks.matrix[,landmarks.order])

## manual method
## landmarks <- read.csv("Corrected Landmark Data.csv")

#### LOAD RAW DATA ####


#### GEOMORPH ANALYSIS ####

## find outliers, plot any x/y points with indiv. ID#

# plot outliers (top left cluster)
landmarks.raw %>%
  mutate(ID=paste(population,individual, sep="_")) %>%
  ggplot(aes(x=x.1, y=y.1, label=ID, color=population))+
    geom_point()+ 
    geom_text()

# remove outliers
landmarks <- landmarks.raw %>%
  filter(y.4<4) %>%
  filter(!is.na(x.1))

landmarks %>% 
  mutate(ID=paste(population,individual, sep="_"))%>%
  ggplot(aes(x=x.5, y=y.5, label=ID, color=population))+
    geom_point()+
    geom_text()

# how did we determine these outliers before?
# landmarks <- landmarks[-c(10,125,271,253,327,360,448,544,630,641,665,877,882,887,898:900,947,959:961,997, 1061, 1172),]
#landmarks <-landmarks[complete.cases(landmarks[,"x.1"]),]

#Make a 3D array; only x/y cordinates, outliers removed...

landmarks.data <- landmarks %>%
  select(matches("^[x,y]{1}."))

landmark.array <- arrayspecs(landmarks.data, 19, 2) #3D array

GPA.landmarks <- gpagen(landmark.array) #Procrustes analysis

GPA.2D <- two.d.array(GPA.landmarks$coords) #2D Data frame of procrustes coordinates
  
#principle component stuff. 
PCA<- plotTangentSpace(GPA.landmarks$coords, groups= factor(landmarks$population), verbose=T)
PCA$pc.scores #gives the PCA scores

#unbending; just show me PC2 vs PC3
plotTangentSpace(GPA.landmarks$coords, groups=factor(landmarks$species), axis1=2, axis2=3)

#make a prettier plot! 
PC.scores<-as.data.frame(PCA$pc.scores)
ggplot(PC.scores, aes(x=PC3, y=PC4, color=landmarks$species)) + geom_point(size=3) 

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

# so what do the lda axes mean?

# the loadings of the axes
lda.species$scaling

# make a data frame with their original names (x.1, y.1 ,etc.)
landmark.names <- names(landmarks)[-c(1:4)][-39]

lm.loadings <- data.frame(lm.name = landmark.names, 
           ld1.loading = lda.species$scaling[,1],
           ld2.loading = lda.species$scaling[,2])

head(lm.loadings)

######## GENETICS VS. MORPHOLOGY ########

# extract centroids
landmark.centroids <- data.frame(population = landmarks$population, 
                                 id = landmarks$ID, 
                                 species = landmarks$species, 
                                 centroid = GPA.landmarks$Csize)
# plot centroid sizes
landmark.centroids %>%
  ggplot(aes(x = population, y = centroid, color=species)) +
  geom_point()+
  facet_grid(~ species, scale = "free", space = "free")

#read in genetic cluster data

cluster.dat <- read.table(file=list.files(pattern="structure",full.names = TRUE), stringsAsFactors = FALSE)
cluster.dat.reduced <- cluster.dat %>%
  filter(k.value.run == 2) %>%
  select(id,X1,X2,membership)

landmark.centroids.rename <- landmark.centroids
landmark.centroids.rename$id <- landmark.centroids$id %>% gsub("2014","",.) %>% gsub("_","",.)

cluster.morph <- left_join(cluster.dat.reduced, landmark.centroids.rename)

cluster.morph  %>%
  filter(!is.na(population)) %>%
  ggplot(aes(x = population, y = centroid, color=factor(membership))) +
  geom_point(size=3)+
  facet_grid(~ species, scale = "free", space = "free")


