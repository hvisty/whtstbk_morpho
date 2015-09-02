#### F2 white stickleback body size and melanocyte analysis #### 

### Initial loading 

data <- read.csv("F2 white stickleback morphometrics.csv")

data$family <- as.factor(data$family)

data$melanocytes <- as.integer(data$melanocytes)



### Body Size: Length and Depth

morphodat <- data[,c(2:4, 8:11)]

## Graphs 

# Boxplot body depth by family
ggplot(morphodat, aes(x=family, y=body.depth)) +
  geom_boxplot(na.rm=T) +
  labs(
    x = "F2 Family",
    y = "Body Depth from first dorsal spine (cm)") +
  theme_classic() 

# Boxplot standard length by family
ggplot(morphodat, aes(x=family, y=standard.length)) +
  geom_boxplot(na.rm=T) +
  labs(
    x = "F2 Family",
    y = "Standard Length (cm)") +
  theme_classic() 



## Hypothesis testing  

# One-way anova for average population std. length by species 
std.length.lm<-lm(morphodat$standard.length ~ morphodat$family)
anova(std.length.lm)

# One-way anova for average population body depth by species 
depth.lm <- lm(morphodat$body.depth ~ morphodat$family)
anova(depth.lm)





### melanocyte counts analysis ### 

## Graphs 

# melanocyte counts by family
ggplot(data, aes(x=family, y=melanopcytes)) + 
  geom_boxplot(na.rm=T) + 
  labs(
    x = "F2 Family",
    y = "Number of melanocytes") +
  theme_classic() 



## Hypothesis testing: melanocyte counts by family 
melan.lm<-lm(data$melanocytes ~ data$family)
anova(melan.lm)


