## broennimann et al 2012 nice overlap functions, wrapped into more
# straight forward functionality.

# My interpretation of the script that they work through in the example...
# I will make this into functions as and when they are needed.

# Here are the libraries that are needed according to script.R

library(BIOMOD)
library(ade4)
library(adehabitat)
library(sp)
library(gam)
library(MASS)
library(mvtnorm)
library(gbm)
library(dismo)

# Function build here.
# Some of this stuff is implemented in ecospat, but I don't think it's
# the best... e.g. some functions should be within others, and also all
# the functions start with ecospat... weird.
# I am also not sure this deals with multidimensional spatial data?

fitBoennimann <- function(clim1, clim2, extent) {
  # climatic data provided as raster, and extent given.
  # Crop rasters to extent, and then turn into a format as in the 
  # broennimann script.
  # Perhaps test what the format of the climate data is to allow a data.frame
  # or matrix, skip the extract step.
  clim1 <- crop(clim1, extent)
  clim2 <- crop(clim2, extent)
  env1 <- raster::extract(clim1, coordinates(clim1))
  env2 <- raster::extract(clim2, coordinates(clim2))

}

# Start with climate data for each othe study sites. What are these data?
# One trait, or many?
# climate data is x, y (coords), then the values are columns. Data for all
# sites of the study area (i.e. all pixels).
# first read in some climate data, and also just get the bee data (while I am
# at it...)
climate_root <- "/home/hfg/Documents/projects/advent/data/bioclim"
bioclim_now <- stack(paste0(climate_root, "/bioclim_now_cropped.grd"))
layers <- c("bio_1", "bio_4", "bio_5", "bio_13", "bio_15")
bioc <- bioclim_now[[which(names(bioclim_now) %in% layers)]]

# since these bees are at the same place clim a and clim b are the same.
clim1 <- clim2 <- extract(bioc, coordinates(bioc))


bees <- read.csv("/home/hfg/Documents/projects/advent/data/d1c/step_bumblebees/CANPOLIN_2014_05_13_ungrided.csv")


# Start with the functions that are provided in the ESM

occ.desaggragation <-function(df,colxy,colvar=NULL,min.dist,plot=T){
  initial<-df
  train<-initial
  xx<-colxy[1]
  yy<-colxy[2]
  kept<-0 ;out<-0; keep<-c()
  x11(2,2,pointsize = 12); par(mar=c(0,0,0,0)); plot.new()

  while(nrow(train)>0){
      
    i<-sample(1:nrow(train),1)
      
    if(sum(sqrt((train[,xx]-train[i,xx])^2 + (train[,yy]-train[i,yy])^2)<=min.dist)>1) {
      out<-out+1
      plot.new(); text(0.5,0.8,paste("# initial:",nrow(initial))); text(0.5,0.5,paste("# kept: ",kept)); text(0.5,0.2,paste("# out: ",out))
    }
    else {
      keep<-c(keep,row.names(train[i,]))
      kept<-kept+1
      plot.new(); text(0.5,0.8,paste("# initial:",nrow(initial))); text(0.5,0.5,paste("# kept: ",kept)); text(0.5,0.2,paste("# out: ",out))
    }

    train<-train[-i,]
  }
  keep.row<-rep(F,nrow(initial))

  for(k in 1:nrow(initial)){
    if( sum(row.names(initial)[k]==keep)==1) keep.row[k]<-T
  }
  dev.off()

  if(is.null(colvar))final<-initial[keep.row,colxy]
  if(ncol(df)==2)final<-initial[keep.row,colxy]
  if(!is.null(colvar)&ncol(df)>2)final<-initial[keep.row,c(colxy,colvar)]

  if(plot==T){
    x11()
    plot(initial[,colxy],main="distribution of occurences",sub=paste("# initial (black):",nrow(initial)," | # kept (red): ",kept),pch=19,col="black",cex=0.2)
    points(final[,colxy],pch=19,col="red",cex=0.2)
  }
  return(final)
}

##################################################################################################
##written by Olivier Broennimann. Departement of Ecology and Evolution (DEE). 
##October 09. University of Lausanne. Switzerland
##
##DESCRIPTION
##
## add environmental values to a species dataframe.
## the xy (lat/long) coordinates of the species occurrences are compared to those of the environment dataframe
## and the value of the closest pixel is added to the species dataframe. 
## when the closest environment pixel is more distant than resolution, NA is added instead of the value.
## (similar to sample() in ArcGIS)

##ARGUMENTS
##dfsp: species dataframe with x, y and optional other variables
##colspxy: the range of columns for x and y in dfsp
##colspkept: the columns of dfsp that should be kept in the final dataframe (by default: xy )
##dfvar: environmental dataframe with x, y and environmental variables
##colvarxy: the range of columns for x and y in dfvar
##colvar: the range of enviromental variables columns in dfvar. (by default: all exept xy )
##resolution: distance between x,y of species and environmental datafreme after which values shouldn't be added 
##(typically, the resolution of the data in dfvar)

sample.sp.globvar <-function(dfsp,colspxy,colspkept="xy",dfvar,colvarxy,colvar="all",resolution){
  if(sum(colspkept=="xy")==1)colspkept<-colspxy
  if(sum(colvar=="all")==1) {
    if(!is.null(colspkept)) colvar<-(1:ncol(dfvar))[-colvarxy]
    if(is.null(colspkept))  colvar<-(1:ncol(dfvar))
  }
  colspx<-colspxy[1];colspy<-colspxy[2];colvarx<-colvarxy[1];colvary<-colvarxy[2]

  x<-dfsp[,colspx]
  X<-dfvar[,colvarx]
  y<-dfsp[,colspy]
  Y<-dfvar[,colvary]

  train<-data.frame(matrix(nrow=nrow(dfsp),ncol=length(colvar)))
  names(train)<-names(dfvar)[colvar]

  x11(2,2,pointsize = 12); par(mar=c(0,0,0,0));
  for (i in 1:nrow(dfsp)){
    dist<-sqrt((X-x[i])^2 + (Y-y[i])^2)
    min<-min(dist)
    if(min<=resolution){
      if(length(colvar)>1)train[i,]<-dfvar[dist==min,colvar][1,]
      if(length(colvar)==1) train[i,]<-dfvar[dist==min,colvar][1]
    }
    plot.new(); text(0.5,0.5,paste(paste("sampling:","\n","runs to go: ",nrow(dfsp)-i))); 
  }
  dev.off()

  if(!is.null(colspkept))final<-cbind(dfsp[,colspkept],train)
  if(is.null(colspkept))final<-train

  return(final)
}

##################################################################################################
##written by Olivier Broennimann. Departement of Ecology and Evolution (DEE). 
##October 09. University of Lausanne. Switzerland
##
##DESCRIPTION
##Investigate spatial autocorrelation by drawing a mantel Correlogram (autocorrelation vs distance)
##
##ARGUMENTS
##df: dataframe with x, y, and variables
##colxy: the range of columns for x and y in df
##colvar: the range of columns for variables in df
##n: number of random occurences used for the test (computation time increase tremendiously when using more than 500occ.)   
##max: maximum distance to be computed in the correlogram
##nclass: number of class of distance to be computed in the correlogram
##nperm: number of permutation in the randomization process

mantel.correlogram <- function(df,colxy,n,colvar,max,nclass,nperm){
  #TODO: proper references here, instead of library
  library(ecodist)

  envnorm<-data.frame(t((t(df[,colvar])-mean(df[,colvar]))/sd(df[,colvar])))
  row.rand<-sample(1:nrow(df),n,replace=T)
  envdist<-dist(envnorm[row.rand,])
  geodist<-dist(df[row.rand,colxy])
  b<- seq(from = min(geodist), to = max, length.out = nclass)
  crlg<-mgram(envdist,geodist,breaks=b,nperm=nperm)
  plot(crlg)
  abline(h=0)
}

##################################################################################################
##written by Olivier Broennimann. Departement of Ecology and Evolution (DEE). 
##October 09. University of Lausanne. Switzerland
##
##DESCRIPTION
##randomly sample pseudoabsences from an environmental dataframe covering the study area
##A minimum distance from presences can be set.
##ARGUMENTS
##nbabsences: number of pseudoabsences desired 
##glob: environmental dataframe covering the study area to sample, with x,y 
##colxyglob: the range of columns for x and y in glob
##colvar: the range of columns for x and y in glob. colvar="all" keeps all the variables in glob in the final dataframe. colvar=NULL keeps only x and y
##presence: occurence dataframe 
##colxypresence: the range of columns for x and y in presence
##mindist: minimum distance from prensences closer to wich pseudoabsences shouldn't be drawn (buffer distance around presences)

rand.pseudoabsences<-function(nbabsences, glob, colxyglob,colvar="all", presence, colxypresence, mindist){
  colxglob<-colxyglob[1]
  colyglob<-colxyglob[2]
  colxpresence<-colxypresence[1]
  colypresence<-colxypresence[2]

  keep<-c()

  no.i<-1
  while(no.i <= nbabsences){
    ki<-sample(1:nrow(glob),1)
    if(sum(((glob[ki,colxglob]- presence[,colxpresence])^2 + (glob[ki,colyglob]- presence[,colypresence])^2) <= mindist^2)==0) {
      keep[no.i]<-ki
      no.i<-no.i+1
    }
  }
  if(sum(colvar=="all")==1) colvar<-(1:ncol(glob))[-colvarxy]
  if(!is.null(colvar))pseudoabs<-glob[keep,c(colxyglob,colvar)]
  if(is.null(colvar))pseudoabs<-glob[keep,colxyglob]

  return(pseudoabs)
}

## Written by Olivier Broennimann. Departement of Ecology and Evolution (DEE). 
## University of Lausanne. Switzerland. October 09.
##
## DESCRIPTION
##
## functions to perform measures of niche overlap and niche equivalency/similarity tests as described in Broennimann et al. (submitted)
## 
## list of functions:
##
## grid.clim(glob,glob1,sp,R) 
## use the scores of an ordination (or SDM predictions) and create a grid z of RxR pixels 
## (or a vector of R pixels when using scores of dimension 1 or SDM predictions) with occurrence densities
## Only scores of one, or two dimensions can be used 
## sp= scores for the occurrences of the species in the ordination, glob = scores for the whole studies areas, glob 1 = scores for the range of sp 
##
## niche.overlap(z1,z2,cor)
## calculate the overlap metrics D and I (see Warren et al 2008) based on two species occurrence density grids z1 and z2 created by grid.clim
## cor=T correct occurrence densities of each species by the prevalence of the environments in their range
##
## niche.equivalency.test(z1,z2,rep)
## runs niche equivalency test(see Warren et al 2008) based on two species occurrence density grids
## compares the observed niche overlap between z1 and z2 to overlaps between random niches z1.sim and z2.sim.
## z1.sim and z2.sim are built from random reallocations of occurences of z1 and z2
## rep is the number of iterations
##
## niche.similarity.test(z1,z2,rep)
## runs niche similarity test(see Warren et al 2008) based on two species occurrence density grids
## compares the observed niche overlap between z1 and z2 to overlaps between z1 and random niches (z2.sim) in available in the range of z2 (z2$Z) 
## z2.sim have the same patterns as z2 but their center are randomly translatated in the availabe z2$Z space and weighted by z2$Z densities
## rep is the number of iterations
##
## plot.niche(z,title,name.axis1,name.axis2)
## plot a niche z created by grid.clim. title,name.axis1 and name.axis2 are strings for the legend of the plot
##
## plot.contrib(contrib,eigen)
## plot the contribution of the initial variables to the analysis. Typically the eigen vectors and eigen values in ordinations
##
## plot.overlap.test(x,type,title)
## plot an histogram of observed and randomly simulated overlaps, with p-values of equivalency and similarity tests. 
## x must be an object created by niche.similarity.test or niche.equivalency.test.
## type is either "D" or "I". title is the title of the plot

##################################################################################################

grid.clim<-function(glob,glob1,sp,R){

  # glob: global background dataset for the whole study area, 
  # glob1: background for sp1
  # sp: occurrence dataset
  # R: resolution of the grid
  l<-list()

  if (ncol(glob)>2) stop("cannot calculate overlap with more than two axes")

  if(ncol(glob)==1){                      #if scores in one dimension (e.g. LDA,SDM predictions,...)
    xmax<-max(glob[,1])
    xmin<-min(glob[,1])
    sp.dens<-density(sp[,1],kernel="gaussian",from=xmin,to=xmax,n=R,cut=0)    # calculate the density of occurrences in a vector of R pixels along the score gradient
                              # using a gaussian kernel density function, with R bins.
    glob1.dens<-density(glob1[,1],kernel="gaussian",from=xmin,to=xmax,n=R,cut=0)  # calculate the density of environments in glob1
    x<-sp.dens$x                      # breaks on score gradient
    z<-sp.dens$y*nrow(sp)/sum(sp.dens$y)              # rescale density to the number of occurrences in sp
                              # number of occurrence/pixel
    Z<-glob1.dens$y*nrow(glob)/sum(glob1.dens$y)            # rescale density to the number of sites in glob1
    z[z<max(z)/1000]<-0                     # remove infinitesimally small number generated by kernel density function
    Z[Z<max(Z)/1000]<-0                     # remove infinitesimally small number generated by kernel density function

    z.uncor<-z/max(z)                     # rescale between [0:1] for comparison with other species
    z<-z/Z                        # correct for environment prevalence 
    z[is.na(z)]<-0                      # remove n/0 situations
    z[z=="Inf"]<-0                      # remove 0/0 situations
    z.cor<-z/max(z)                     # rescale between [0:1] for comparison with other species
    l$x<-x;l$z.uncor<-z.uncor;l$z.cor<-z.cor;l$Z<-Z;l$glob<-glob;l$glob1<-glob1;l$sp<-sp 
  }

  if(ncol(glob)==2){ #if scores in two dimensions (e.g. PCA)          


    library(adehabitat)
        xmin<-min(glob[,1]);xmax<-max(glob[,1]);ymin<-min(glob[,2]);ymax<-max(glob[,2])     # data preparation  
    glob1r<-data.frame(cbind((glob1[,1]-xmin)/abs(xmax-xmin),(glob1[,2]-ymin)/abs(ymax-ymin)))  # data preparation  
    spr<-data.frame(cbind((sp[,1]-xmin)/abs(xmax-xmin),(sp[,2]-ymin)/abs(ymax-ymin)))       # data preparation
    mask<-ascgen(cbind((1:R)/R,(1:R)/R),nrcol=R,count=F)                # data preparation
    sp.dens<-kernelUD(spr[,1:2],grid=mask,kern="bivnorm")         # calculate the density of occurrences in a grid of RxR pixels along the score gradients
                              # using a gaussian kernel density function, with RxR bins.
    sp.dens<-asc2spixdf(sp.dens[[1]]$UD)              # data manipulation                                   
    glob1.dens<-kernelUD(glob1r[,1:2],grid=mask,kern="bivnorm")
    glob1.dens<-asc2spixdf(glob1.dens[[1]]$UD)
    x<-seq(from=min(glob[,1]),to=max(glob[,1]),length.out=R)        # breaks on score gradient 1 
    y<-seq(from=min(glob[,2]),to=max(glob[,2]),length.out=R)        # breaks on score gradient 2 
    z<-matrix(sp.dens$var*nrow(sp)/sum(sp.dens$var),nrow=R,ncol=R,byrow=F)      #rescale density to the number of occurrences in sp
    Z<-matrix(glob1.dens$var*nrow(glob1)/sum(glob1.dens$var),nrow=R,ncol=R,byrow=F)   #rescale density to the number of sites in glob1
    z[z<max(z)/1000]<-0                       # remove infinitesimally small number generated by kernel density function
    Z[Z<max(Z)/1000]<-0                       # remove infinitesimally small number generated by kernel density function

    z.uncor<-z/max(z)                     # rescale between [0:1] for comparison with other species 
    z<-z/Z                        # correct for environment prevalence
    z[is.na(z)]<-0                      # remove n/0 situations
    z[z=="Inf"]<-0                      # remove n/0 situations
    z.cor<-z/max(z)                     # rescale between [0:1] for comparison with other species 
    l$x<-x;l$y<-y;l$z.uncor<-z.uncor;l$z.cor<-z.cor;l$Z<-Z;l$glob<-glob;l$glob1<-glob1;l$sp<-sp
  }
  return(l)
}

##################################################################################################
niche.overlap<-function(z1,z2,cor){ 

  # z1 = species 1 occurrence density grid created by grid.clim
  # z2 = species 2 occurrence density grid created by grid.clim
  # cor=T correct occurrence densities of each species by the prevalence of the environments in their range

  l<-list()

  if(cor==F){   
  p1<-z1$z.uncor/sum(z1$z.uncor) # rescale occurence densities so that the sum of densities is the same for both species
  p2<-z2$z.uncor/sum(z2$z.uncor) # rescale occurence densities so that the sum of densities is the same for both species
  }

  if(cor==T){
  p1<-z1$z.cor/sum(z1$z.cor)  # rescale occurence densities so that the sum of densities is the same for both species
  p2<-z2$z.cor/sum(z2$z.cor)  # rescale occurence densities so that the sum of densities is the same for both species
  }

  D <-1-(0.5*(sum(abs(p1-p2))))       # overlap metric D
  I <-1-(0.5*(sqrt(sum((sqrt(p1)-sqrt(p2))^2))))  # overlap metric I
  l$D<-D
  l$I<-I
  return(l)
}

##################################################################################################

niche.equivalency.test<-function(z1,z2,rep){
  R<-length(z1$x)
  l<-list()
  x11(2,2,pointsize = 12); par(mar=c(0,0,0,0));

  obs.o<-niche.overlap(z1,z2,cor=T)                 #observed niche overlap
  sim.o<-data.frame(matrix(nrow=rep,ncol=2))              #empty list of random niche overlap
  names(sim.o)<-c("D","I")
  for (i in 1:rep){
    plot.new(); text(0.5,0.5,paste("runs to go:",rep-i+1))
    
    if(is.null(z1$y)){ #overlap on one axis
    
      occ1.sim<-sample(z1$x,size=nrow(z1$sp),replace=T,prob=z1$z.cor)     #random sampling of occurrences following the corrected densities distribution
      occ2.sim<-sample(z2$x,size=nrow(z2$sp),replace=T,prob=z2$z.cor)
    
      occ.pool<-c(occ1.sim,occ2.sim)   # pool of random occurrences
      rand.row<-sample(1:length(occ.pool),length(occ1.sim),replace=T)     # random reallocation of occurrences to datasets
      sp1.sim<-occ.pool[rand.row]
      sp2.sim<-occ.pool[-rand.row]
    
      z1.sim<-grid.clim(z1$glob,z1$glob1,data.frame(sp1.sim),R) # gridding
      z2.sim<-grid.clim(z2$glob,z2$glob1,data.frame(sp2.sim),R) 
    }

    if(!is.null(z1$y)){                     #overlap on two axes
      coordinates<-which(z1$z.cor>0,arr.ind=T)            # array of cell coordinates ((1,1),(1,2)...)
      weight<-z1$z.cor[z1$z.cor>0]                #densities in the same format as cells
      coordinates.sim1<-coordinates[sample(1:nrow(coordinates),size=nrow(z1$sp),replace=T,prob=weight),] #random sampling of coordinates following z1$z.cor distribution
      occ1.sim<-cbind(z1$x[coordinates.sim1[,1]],z1$y[coordinates.sim1[,2]])  # random occurrences following the corrected densities distribution
    
      coordinates<-which(z2$z.cor>0,arr.ind=T)            # array of cell coordinates ((1,1),(1,2)...)
      weight<-z2$z.cor[z2$z.cor>0]                #densities in the same format as cells
      coordinates.sim2<-coordinates[sample(1:nrow(coordinates),size=nrow(z2$sp),replace=T,prob=weight),] #random sampling of coordinates following z1$z.cor distribution
      occ2.sim<-cbind(z2$x[coordinates.sim2[,1]],z2$y[coordinates.sim2[,2]])  # random occurrences following the corrected densities distribution
    
      occ.pool<-rbind(occ1.sim,occ2.sim) # pool of random occurrences
      rand.row<-sample(1:nrow(occ.pool),nrow(occ1.sim),replace=T)       # random reallocation of occurrences to datasets
      sp1.sim<-occ.pool[rand.row,]
      sp2.sim<-occ.pool[-rand.row,]
    
      z1.sim<-grid.clim(z1$glob,z1$glob1,data.frame(sp1.sim),R)
      z2.sim<-grid.clim(z2$glob,z2$glob1,data.frame(sp2.sim),R) 
    }

    o.i<-niche.overlap(z1.sim,z2.sim,cor=F)             # overlap between random and observed niches
    sim.o$D[i]<-o.i$D                     # storage of overlaps
    sim.o$I[i]<-o.i$I
  }

  dev.off()
  l$sim<-sim.o  # storage
  l$obs<-obs.o  # storage
  l$p.D<- min((sum(sim.o$D <= obs.o$D ) + 1),(sum(sim.o$D >= obs.o$D ) + 1))*2/(length(sim.o$D) + 1)  # storage of p-values
  l$p.I<- min((sum(sim.o$I <= obs.o$I ) + 1),(sum(sim.o$I >= obs.o$I ) + 1))*2/(length(sim.o$I) + 1)  # storage of p-values

  return(l)
}

##################################################################################################

niche.similarity.test<-function(z1,z2,rep){
  R<-length(z1$x)
  x11(2,2,pointsize = 12); par(mar=c(0,0,0,0));           # countdown window
  l<-list()
  obs.o<-niche.overlap(z1,z2,cor=T)                 #observed niche overlap
  sim.o<-data.frame(matrix(nrow=rep,ncol=2))            #empty list of random niche overlap
  names(sim.o)<-c("D","I")

  for (k in 1:rep){
    plot.new(); text(0.5,0.5,paste("similarity tests:","\n","runs to go:",rep-k+1)) # countdown
    
    if(is.null(z2$y)){
      center<-which(z2$z.cor==1,arr.ind=T)          # define the centroid of the observed niche
      Z<-z2$Z/max(z2$Z)
      rand.center<-sample(1:R,size=1,replace=F,prob=Z)        # randomly (weighted by environment prevalence) define the new centroid for the niche
      
      xshift<-rand.center-center              # shift on x axis
      z2.sim<-z2
      z2.sim$z.cor<-rep(0,R)                # set intial densities to 0
      for(i in 1:R){
        i.trans<-i+xshift
        if(i.trans>R|i.trans<0)next()           # densities falling out of the env space are not considered
        z2.sim$z.cor[i.trans]<-z2$z.cor[i]          # shift of pixels
      }
      z2.sim$z.cor<-(z2$Z!=0)*1*z2.sim$z.cor          # remove densities out of existing environments
    }

    if(!is.null(z2$y)){
      centroid<-which(z2$z.cor==1,arr.ind=T)[1,]        # define the centroid of the observed niche
      Z<-z2$Z/max(z2$Z)
      rand.centroids<-which(Z>0,arr.ind=T)          # all pixels with existing environments in the study area
      weight<-Z[Z>0]
      rand.centroid<-rand.centroids[sample(1:nrow(rand.centroids)
        ,size=1,replace=F,prob=weight),]          # randomly (weighted by environment prevalence) define the new centroid for the niche
      xshift<-rand.centroid[1]-centroid[1]          # shift on x axis
      yshift<-rand.centroid[2]-centroid[2]          # shift on y axis
      z2.sim<-z2
      z2.sim$z.cor<-matrix(rep(0,R*R),ncol=R,nrow=R)        # set intial densities to 0
      for(i in 1:R){
        for(j in 1:R){
          i.trans<-i+xshift
          j.trans<-j+yshift
          if(i.trans>R|i.trans<0)next()         # densities falling out of the env space are not considered
          if(j.trans>R|j.trans<0)next()
          z2.sim$z.cor[i.trans,j.trans]<-z2$z.cor[i,j]    # shift of pixels
        }
      }
      z2.sim$z.cor<-(z2$Z!=0)*1*z2.sim$z.cor          # remove densities out of existing environments
    }
    o.i<-niche.overlap(z1,z2.sim,cor=T)             # overlap between random and observed niches
    sim.o$D[k]<-o.i$D                   # storage of overlaps
    sim.o$I[k]<-o.i$I
  }
  dev.off()
  l$sim<-sim.o                      # storage
  l$obs<-obs.o                      # storage
  l$p.D<- min((sum(sim.o$D <= obs.o$D ) + 1),(sum(sim.o$D >= obs.o$D ) + 1))*2/(length(sim.o$D) + 1)  # storage of p-values
  l$p.I<- min((sum(sim.o$I <= obs.o$I ) + 1),(sum(sim.o$I >= obs.o$I ) + 1))*2/(length(sim.o$I) + 1)  # storage of p-values

  return(l)
}

##################################################################################################

plot.niche<-function(z,title,name.axis1,name.axis2,cor=F){

  if(is.null(z$y)){
    R<-length(z$x)
    x<-z$x
    xx<-sort(rep(1:length(x),2))

    if(cor==F)y1<-z$z.uncor/max(z$z.uncor)
    if(cor==T)y1<-z$z.cor/max(z$z.cor)
    Y1<-z$Z/max(z$Z)
    yy1<-sort(rep(1:length(y1),2))[-c(1:2,length(y1)*2)]
    YY1<-sort(rep(1:length(Y1),2))[-c(1:2,length(Y1)*2)]

    plot(x,y1,type="n",xlab=name.axis1,ylab="density of occurrence")
    polygon(x[xx],c(0,y1[yy1],0,0),col="grey")
    lines(x[xx],c(0,Y1[YY1],0,0))
  }

  if(!is.null(z$y)){
    if(cor==F)image(z$x,z$y,z$z.uncor,col = gray(100:0 / 100),zlim=c(0.000001,max(z$z.uncor)),xlab=name.axis1,ylab=name.axis2)
    if(cor==T)image(z$x,z$y,z$z.cor,col = gray(100:0 / 100),zlim=c(0.000001,max(z$z.cor)),xlab=name.axis1,ylab=name.axis2)
    contour(z$x,z$y,z$Z,add=T,levels=quantile(z$Z[z$Z>0],c(0,0.5)),drawlabels=F,lty=c(1,2))
  }
  title(title)
}

plot.contrib<-function(contrib,eigen){

  if(ncol(contrib)==1){
        h<-c(unlist(contrib))
        n<-row.names(contrib)
        barplot(h,space=0,names.arg=n)
        title(main="variable contribution")
    }
  if(ncol(contrib)==2){
          s.corcircle(contrib[,1:2]/max(abs(contrib[,1:2])),grid = F)
          title(main="correlation circle", sub=paste("axis1 = ",round(eigen[1]/sum(eigen)*100,2),"%","axis2 = ",round(eigen[2]/sum(eigen)*100,2),"%"))
    }
}


plot.overlap.test<-function (x,type,title) {
  if(type=="D") {
    obs <- x$obs$D
    sim <- x$sim$D
    p<-x$p.D
  }
  if(type=="I") {
    obs <- x$obs$I
    sim <- x$sim$I
    p<-x$p.I
  }
  r0 <- c(sim, obs)
      l0 <- max(sim) - min(sim)
      w0 <- l0/(log(length(sim), base = 2) + 1)
      xlim0 <- range(r0) + c(-w0, w0)
      h0 <- hist(sim, plot = FALSE, nclass = 10)
      y0 <- max(h0$counts)
      hist(sim, plot = TRUE, nclass = 10, xlim = xlim0, col = grey(0.8), main= title, xlab=type,sub = paste("p.value = ",round(p,5)))
      lines(c(obs, obs), c(y0/2, 0),col="red")
      points(obs, y0/2, pch = 18, cex = 2,col="red")
      invisible()
}

##################################################################################################



