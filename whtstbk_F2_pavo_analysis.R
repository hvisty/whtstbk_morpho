####spectral analysis of f2 data
####using 'pavo' library by r. maia
####KMS 2014 --- edited HV 2015 

# required libraries
library(pavo)

# set this to where ever you choose to put these files
setwd("/Users/Hannah/Documents/spec_analysis/")

# read in spectral measurements
# separated uv (200-300 nm) from rest of data
specs <- getspec(where = "spec_data", ext = "txt", lim = c(300, 700))
specs.uv <- getspec(where ="spec_data", ext = "txt", lim = c(200, 300))



### NON UV DATA

# plot first 200 groups of 3 (measurements taken 3 times per landmark)
explorespec(specs[, 1:200], by = 3, lwd = 1)

#average every 3 measures
mspecs <- aggspec(specs, by = 3, FUN = mean)

#take a gander at the averaged data
explorespec(mspecs[, 1:200], by = 1, lwd = 1)

#plot spec data + 4 levels of smoothing to assess need for smoothing
plotsmooth(mspecs[, 1:200], minsmooth = 0.05, maxsmooth = 0.5, curves = 4, ask = F)

#apply smoothing
spec.sm <- procspec(mspecs, opt='smooth', span=0.4)

#take a gander
explorespec(spec.sm, by=1,lwd=1)

#aggplot(spec.sm,by=5)

#find peak reflectance measures
# B3 = height of reflectance peak
# H1 = location (in nm) of peak
# FWHM = width of peak (in nm)
peaks <- peakshape(spec.sm)



### UV DATA
# plot first 200 groups of 3 (measurements taken 3 times per landmark)
explorespec(specs.uv[, 1:200], by = 3, lwd = 1)

#average every 3 measures
mspecs.uv <- aggspec(specs.uv, by = 3, FUN = mean)

#take a gander at the averaged data
explorespec(mspecs.uv[, 1:200], by = 1, lwd = 1)

#plot spec data + 4 levels of smoothing to assess need for smoothing
plotsmooth(mspecs.uv[, 1:200], minsmooth = 0.05, maxsmooth = 0.5, curves = 4, ask = F)

#apply smoothing
spec.sm.uv <- procspec(mspecs.uv, opt='smooth', span=0.2)

#take a gander
explorespec(spec.sm.uv, by=1,lwd=1)

#aggplot(spec.sm.uv,by=5)

#find peak reflectance measures
# B3 = height of reflectance peak
# H1 = location (in nm) of peak
# FWHM = width of peak (in nm)
peaks.uv <- peakshape(spec.sm.uv, lim=c(200, 300))


### Process to compare spec data with melanophore data 

## Change spec data into RGB, then HSV color scale for easier understanding 
# spec data goes to hexidecimal with pavo 
spec.hex <- spec2rgb(spec.sm)

# change to RGB then HSV 
spec.rgb <- col2rgb(spec.hex)
spec.hsv <- rgb2hsv(spec.rgb)


## Split color data IDs into matching format for merging data 

# transform data and split IDs into family, individual, spec number
spec.hsv <- as.data.frame(t(spec.hsv))

spec.hsv$ID <- rownames(spec.hsv)
rownames(spec.hsv) <- NULL

spec.hsv$ID <- sub("[.]", "-", spec.hsv$ID)
spec.hsv$ID <- sub("[.]", "-", spec.hsv$ID)

spec.hsv$family <- as.character(lapply(strsplit(spec.hsv$ID, split="-"), "[", 2))
spec.hsv$individual <- as.character(lapply(strsplit(spec.hsv$ID, split="-"), "[", 3))
spec.hsv$individual <- gsub("[.]\\d", "", spec.hsv$individual)

spec.hsv$spec.number <- as.numeric(lapply(strsplit(spec.hsv$ID, split="[.]"), "[", 2))
spec.hsv$spec.number[is.na(spec.hsv$spec.number)] <- 0

spec.hsv$ID <- gsub("[.]\\d", "", spec.hsv$ID)

spec.hsv<-spec.hsv[,c(4:7, 1:3)]

# transform each spec measurement point (0:4) into its own column for H, S, V
spec.hsv0 <- subset(spec.hsv, spec.number==0, select=c(ID, family, individual, h, s, v))
spec.hsv1 <- subset(spec.hsv, spec.number==1, select=c(ID, h, s, v))
spec.hsv2 <- subset(spec.hsv, spec.number==2, select=c(ID, h, s, v))
spec.hsv3 <- subset(spec.hsv, spec.number==3, select=c(ID, h, s, v))
spec.hsv4 <- subset(spec.hsv, spec.number==4, select=c(ID, h, s, v))

spec.hsv <- merge(spec.hsv0, spec.hsv1, by="ID")
spec.hsv <- merge(spec.hsv, spec.hsv2, by="ID")
spec.hsv <- merge(spec.hsv, spec.hsv3, by="ID")
spec.hsv <- merge(spec.hsv, spec.hsv4, by="ID")

colnames(spec.hsv) <- c("ID", "family", "individual",
                        "H.0", "S.0", "V.0",
                        "H.1", "S.1", "V.1",
                        "H.2", "S.2", "V.2",
                        "H.3", "S.3", "V.3",
                        "H.4", "S.4", "V.4")


# read in morphometric and melanocyte data frame 
F2.morphodat <- read.csv("/Users/Hannah/Documents/GitHub/whtstbk_morpho/F2 white stickleback morphometrics.csv")
F2.morphodat <- F2.morphodat[,c(2:11)]


# merge spec (hsv dat) with morphodat 

spec.hsv$family <- as.numeric(spec.hsv$family)
F2.dat<-left_join(F2.morphodat, spec.hsv, by=c("ID", "family", "individual"))

F2.dat$melanophores<-as.numeric(F2.dat$melanophores)



## Plots: HSV data compared to body size and melanophore counts 

# Value (brightness) by family
# V.0 -V.4 for different measurement locations on fish
ggplot(spec.hsv, aes(x=family, y=V.1, color=factor(family))) +
  geom_point(na.rm=T, size=5) + 
  labs( x= "Family",
        y= "Value (brightness)",
        color= "Family") +
  theme_classic()

# Hue by family
# H.0 -H.4 for different measurement locations on fish
ggplot(spec.hsv, aes(x=family, y=H.0, color=factor(family))) +
  geom_point(na.rm=T, size=5) + 
  labs( x= "Family",
        y= "Hue",
        color= "Family") +
  theme_classic()

# Standard Length by Value (brightness) 
# V.0 -V.4 for different measurement locations on fish
ggplot(F2.dat, aes(x=standard.length, y=V.2, color=factor(family))) + 
  geom_point(na.rm=T, size=5) + 
  labs( x= "Standard Length", 
        y= "Value (brightness)", 
        color= "Family") + 
  theme_classic()

# Body Depth by Value (brightness) 
# V.0 -V.4 for different measurement locations on fish
ggplot(F2.dat, aes(x=body.depth, y=V.4, color=factor(family))) + 
  geom_point(na.rm=T, size=5) + 
  labs( x= "Body Depth", 
        y= "Value (brightness)", 
        color= "Family") + 
  theme_classic()

# Value (brightness) vs melanophores
# V.0 -V.4 for different measurement locations on fish
ggplot(F2.dat, aes(x=(V.4*S.4), y=melanophores, color=factor(family))) +
  geom_point(na.rm=T, size=5) + 
  labs( x= "Value (brightness)",
        y= "# melanohphores",
        color= "Family") +
  theme_classic()

# Hue measures (color/wavelength) vs melanophores
# H.0-H.4 for different measurement points on fish
ggplot(F2.dat, aes(x=H.4, y=melanophores, color=factor(family))) +
  geom_point(na.rm=T, size=5) + 
  labs( x= "Hue",
        y= "# melanohphores",
        color= "Family") +
  theme_classic()




#### UV DATA 
### Process to compare spec data with melanophore data 

# split uv data IDs into matching format for merging data 
peaks.uv$ID <- as.character(peaks.uv$id)
peaks.uv$ID <- sub("[.]", "-", peaks.uv$id)
peaks.uv$ID <- sub("[.]", "-", peaks.uv$ID)

peaks.uv$family <- as.character(lapply(strsplit(peaks.uv$ID, split="-"), "[", 2))
peaks.uv$individual <- as.character(lapply(strsplit(peaks.uv$ID, split="-"), "[", 3))
peaks.uv$individual <- gsub("[.]\\d", "", peaks.uv$individual)

peaks.uv$spec.number <- as.numeric(lapply(strsplit(peaks.uv$ID, split="[.]"), "[", 2))
peaks.uv$spec.number[is.na(peaks.uv$spec.number)] <- 0

peaks.uv$ID <- gsub("[.]\\d", "", peaks.uv$ID)

peaks.uv<-peaks.uv[,c(8:11, 2:7)]


# transform each measurement point into its own column for B3, H1, FWHM
uv.dat0 <- subset(peaks.uv, spec.number== 0, select=c(ID, family, individual, B3, H1, FWHM))
uv.dat1 <- subset(peaks.uv, spec.number==1, select=c(ID, B3, H1, FWHM))
uv.dat2 <- subset(peaks.uv, spec.number==2, select=c(ID, B3, H1, FWHM))
uv.dat3 <- subset(peaks.uv, spec.number==3, select=c(ID, B3, H1, FWHM))
uv.dat4 <- subset(peaks.uv, spec.number==4, select=c(ID, B3, H1, FWHM))

uv.dat <- merge(uv.dat0, uv.dat1, by="ID")
uv.dat <- merge(uv.dat, uv.dat2, by="ID")
uv.dat <- merge(uv.dat, uv.dat3, by="ID")
uv.dat <- merge(uv.dat, uv.dat4, by="ID")

colnames(uv.dat) <- c("ID", "family", "individual", 
                        "B3.0", "H1.0", "FWHM.0",
                        "B3.1", "H1.1", "FWHM.1",
                        "B3.2", "H1.2", "FWHM.2",
                        "B3.3", "H1.3", "FWHM.3",
                        "B3.4", "H1.4", "FWHM.4")

# merge with morphodat 

uv.dat$family <- as.numeric(uv.dat$family)

F2.uvdat<-left_join(F2.morphodat, uv.dat, by=c("ID", "family", "individual"))
F2.uvdat$melanophores<-as.numeric(F2.uvdat$melanophores)


## Plots: HSV data compared to body size and melanophore counts 

# B3 (brightness measure) by family 
# B3.0 -B3.4 for different measurement locations on fish
ggplot(uv.dat, aes(x=family, y=B3.2, color=factor(family))) +
  geom_point(na.rm=T, size=5) + 
  labs( x= "Family",
        y= "Brightness measure",
        color= "Family") +
  theme_classic()


# Standard Length by Brightness 
# B3.0 -B3.4 for different measurement locations on fish
ggplot(F2.uvdat, aes(x=standard.length, y=B3.4, color=factor(family))) + 
  geom_point(na.rm=T, size=5) + 
  labs( x= "Standard Length", 
        y= "Brightness", 
        color= "Family") + 
  theme_classic()


# Brightness vs melanophores
# B3.0 -B3.4 for different measurement locations on fish
ggplot(F2.uvdat, aes(x=(B3.4*FWHM.4), y=melanophores, color=factor(family))) +
  geom_point(na.rm=T, size=5) + 
  labs( x= "brightness (peak height) * peak width",
        y= "# melanophores",
        color= "Family") +
  theme_classic()






##### Postscript: Analysis with spec data untransformed ##### 
### Compare spec data with melanocyte data 


# split spec data IDs into matching format for merging data 
peaks$ID <- as.character(peaks$id)
peaks$ID <- sub("[.]", "-", peaks$id)
peaks$ID <- sub("[.]", "-", peaks$ID)

peaks$family <- as.character(lapply(strsplit(peaks$ID, split="-"), "[", 2))
peaks$individual <- as.character(lapply(strsplit(peaks$ID, split="-"), "[", 3))
peaks$individual <- gsub("[.]\\d", "", peaks$individual)

peaks$spec.number <- as.numeric(lapply(strsplit(peaks$ID, split="[.]"), "[", 2))
peaks$spec.number[is.na(peaks$spec.number)] <- 0

peaks$ID <- gsub("[.]\\d", "", peaks$ID)

peaks<-peaks[,c(8:11, 2:7)]


# transform each measurement point into its own column for B3, H1, FWHM
peak.dat <- subset(peaks, spec.number== 0, select=c(ID, B3, H1, FWHM))
peak.dat1 <- subset(peaks, spec.number==1, select=c(ID, B3, H1, FWHM))
peak.dat2 <- subset(peaks, spec.number==2, select=c(ID, B3, H1, FWHM))
peak.dat3 <- subset(peaks, spec.number==3, select=c(ID, B3, H1, FWHM))
peak.dat4 <- subset(peaks, spec.number==4, select=c(ID, B3, H1, FWHM))

peak.dat <- merge(peak.dat, peak.dat1, by="ID")
peak.dat <- merge(peak.dat, peak.dat2, by="ID")
peak.dat <- merge(peak.dat, peak.dat3, by="ID")
peak.dat <- merge(peak.dat, peak.dat4, by="ID")

colnames(peak.dat) <- c("ID", "B3.0", "H1.0", "FWHM.0",
                        "B3.1", "H1.1", "FWHM.1",
                        "B3.2", "H1.2", "FWHM.2",
                        "B3.3", "H1.3", "FWHM.3",
                        "B3.4", "H1.4", "FWHM.4")

# merge with morphodat 
F2.dat<-left_join(F2.morphodat, peak.dat, by="ID")

F2.dat$melanophores<-as.numeric(F2.dat$melanophores)