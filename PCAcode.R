dat<- read.table("pastmatrix4.11.txt", header=TRUE, sep="\t", dec=".", strip.white=TRUE) # loads the data
info<- read.table("info.txt", header=TRUE, sep="\t", dec=".", strip.white=TRUE) # loads the data

library("stats", lib.loc="C:/Program Files/R/R(-3.2.0/library")
library("dplyr", lib.loc="~/R/win-library/3.2")
library("devtools", lib.loc="~/R/win-library/3.2")
library("sinkr")

#merge data
datinfo<-full_join(dat,info, by="sample")

# log transform 
log.dat <- log(datinfo[1:89, 2:13])
dat.species <- datinfo[1:89, 14]
dat.lake <- datinfo[1:89, 15]
dat.date <- datinfo[1:89, 16]
dat.type<-datinfo[1:89, 17]
dat.type<-as.factor(dat.type)

### Interpolation, EOF, and Reconstruction ###
#Interpolation with DINEOF
log.dat<-as.matrix(log.dat)

din <- dineof(log.dat)
Xa <- din$Xa 
image(Xa)

# EOF
Ed <- eof(Xa) # obs, dineof + lseof

###Reconstruction
Rd <- eofRecon(Ed)

# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
dat.pca <- prcomp(Rd,data=dat.species,
                 center = TRUE,
                 scale. = TRUE, na.rm=TRUE) 
plot(dat.pca)

install_github("vqv/ggbiplot")

library(ggbiplot)
g <- ggbiplot(dat.pca, obs.scale = 1,var.scale=1, groups=dat.species, labels=dat.lake,labels.size=4, ellipse = TRUE, circle=TRUE, var.axes=FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',legend.position = 'top')
print(g)
