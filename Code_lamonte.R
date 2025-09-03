#spatial analysis of shibantan trace fossils Lamonte trevallis
library(fitdistrplus)
library(mclust)
library(methods)
library(plyr)
library(selectspm)
library(spatstat)
library(spatstat.utils)
library(VGAM)
library(XML)
library(rlist)
library(ecespa)


##input data in. svg file with lines on
##run dex filaments from Katie Delahooke: https://github.com/kmdelahooke/dex/blob/main/dex_filaments.R
#--------------------------------------------data preparation ------------------------------------------------
#all data supporting this paper are in supplementary and uploaded to the github
#xlsx file ---> csv files 
dataLine<-read.csv(file.choose())##segments position
dataAll<-read.csv(file.choose())#trace fossil data
edge<-read.csv("Edge.csv")#edge of the surface


##Divide data into groups
data.iso<-dataAll[dataAll[,2]=="isolated",]#isolated vertical burrows
data.con<-dataAll[dataAll[,2]=="connected",]#exclude isolated #connected vertical burrows
data.trail<-dataAll[dataAll[,2]!="isolated",]#all lines


# Plot points for Ediacaran window creation. 
plot(dataLine[,c(1,3)])
points(dataLine[,c(2,4)],add=TRUE)

points(data.iso[,c(4,5)],add=TRUE)
points(data.con[,c(4,5)],add=TRUE)

win.tr.final<-clickpoly(add=TRUE) #create a new window shape
win.tr.final<-owin(poly=edge) #edge file

#create psp line 
psp.line<-psp(dataLine[,1],dataLine[,3],dataLine[,2],dataLine[,4],win.tr.final)

#create ppp of lines' midpoints
ppp.line.len<-ppp(coords(midpoints.psp(psp.line))[,1],coords(midpoints.psp(psp.line))[,2],window=win.tr.final, marks=lengths_psp(psp.line))
ppp.line.ang<-ppp(coords(midpoints.psp(psp.line))[,1],coords(midpoints.psp(psp.line))[,2],window=win.tr.final, marks=angles.psp(psp.line))
ppp.line<-ppp(coords(midpoints.psp(psp.line))[,1],coords(midpoints.psp(psp.line))[,2],window=win.tr.final)

#ppp of vertical burrows
ppp.iso<-ppp(data.iso[,4],data.iso[,5],window=win.tr.final)
ppp.con<-ppp(data.con[,4],data.con[,5],window=win.tr.final)
ppp.trail.len<-ppp(data.trail$DiscX,data.trail$DiscY,window=win.tr.final,marks=data.trail$Length)
ppp.trail.wid<-ppp(data.trail$DiscX,data.trail$DiscY,window=win.tr.final,marks=data.trail$Width)



#--------------------------------------------length width and directions-------------------------------------
  ##list of the Segment lengths (each trace). 
tr.sh.line<-(lengths_psp(psp.line))
tr.line.ang<-(angles.psp(psp.line))

par(mfrow=c(2,2))
hist(marks(ppp.trail.len),main="Trace length",xlab="Trace length (cm)",col="#EDD6D7",xaxs="i",yaxs="i",font.lab=2,border="#ffffff")
hist(marks(ppp.trail.wid),main="Trace width",xlab="Trace width (cm)",col="#EDD6D7",xaxs="i",yaxs="i",font.lab=2,border="#ffffff")
hist(marks(ppp.line.len),main="Segent length",xlab="Segment length (cm)",col="#EDD6D7",xaxs="i",yaxs="i",font.lab=2,border="#ffffff")

  ##roseplot calibrated to the real orientation
angle<-(marks(ppp.line.ang)*180+310)/pi
for (i in 1:2213)
{if (angle[i]>180)
{angle[i]<-angle[i]-180}}
rose(angle,breaks=36,col="#EDD6D7",main="Trace angle")



##--------------------------------------------density plots--------------------------------------------------
  #nice plots #using automation to determine the appropriate density. 
  ##line density
tr.den<-density(psp.line,10)#plot(density(psp.line,20)
plot(tr.den,main="Midpoint density")
plot(psp.line,lwd=1.5,add=TRUE)


  ##points density
tr.den.iso<-density(ppp.iso,bw.diggle(ppp.iso))
plot(density(ppp.iso,bw.diggle(ppp.iso)),legend=FALSE,xaxs="i",yaxs="i",font.lab=2,main="density of isolated vertical burrows")
plot(ppp.iso,add=TRUE)

tr.den.con<-density(ppp.con,bw.diggle(ppp.con))
plot(density(ppp.con,bw.diggle(ppp.con)),legend=FALSE,xaxs="i",yaxs="i",font.lab=2,main="density of connected vertical burrows")
plot(ppp.con,add=TRUE)


#--------------------------------------------PCF-------------------------------------------------------------
#uni PCF####
  ##midpoint and intersection pcf analsyes
plot(envelope(ppp.line.len,pcf,nsim=999, nrank=49),legend=FALSE,xaxs="i",yaxs="i",font.lab=2,main="length pcf",ylim=c(0,3),xlab="Length (cm)",ylab="Midpoint PCF")

  ##isolated plugs pcf
plot(envelope(ppp.iso,pcf,nsim=999, nrank=49),legend=FALSE,xaxs="i",yaxs="i",font.lab=2,main="Isolated vertical burrows PCF",ylim=c(0,3),xlab="Length (cm)",ylab="PCF")

  ##connected plugs pcf
plot(envelope(ppp.con,pcf,nsim=999, nrank=49),legend=FALSE,xaxs="i",yaxs="i",font.lab=2,main="Connected vertical burrows PCF",ylim=c(0,3),xlab="Length (cm)",ylab="PCF")

#Bivariate PCF ####
ppp.iso_line<-superimpose(rest=unmark(ppp.iso), lines1=unmark(ppp.line))
ppp.iso_con<-superimpose(rest=unmark(ppp.iso), lines1=unmark(ppp.con))
ppp.con_line<-superimpose(rest=unmark(ppp.con),lines1=unmark(ppp.line))

plot(envelope(ppp.iso_line,pcfcross,nsim=999, nrank=49),legend=FALSE,xaxs="i",yaxs="i",font.lab=2,ylab="Bivariate PCF",xlab="Distance (cm)",main="Bivariate (isolated vs line) PCF")
plot(envelope(ppp.iso_con,pcfcross,nsim=999, nrank=49),legend=FALSE,xaxs="i",yaxs="i",font.lab=2,ylab="Bivariate PCF",xlab="Distance (cm)",main="Bivariate (isolated vs connected) PCF")
plot(envelope(ppp.con_iso,pcfcross,nsim=999, nrank=49),legend=FALSE,xaxs="i",yaxs="i",font.lab=2,ylab="Bivariate PCF",xlab="Distance (cm)",main="Bivariate (connected vs isolated) PCF")

par(mfrow=c(1,2))

#comparisons with CSR
fit.line.an2<-envelope(ppm(unmark(ppp.line.ang), trend=~1,statistic="pcf"),Gest,savefuns=TRUE,nsim=999,rank=50)
fit.con_line<-envelope(ppm(unmark(ppp.con_line), trend=~1,statistic="pcf"),Gest,savefuns=TRUE,nsim=999,rank=50)
fit.iso_line<-envelope(ppm(unmark(ppp.iso_line), trend=~1,statistic="pcf"),Gest,savefuns=TRUE,nsim=999,rank=50)

LF.gof(fit.line.an2)
LF.gof(fit.con_line)
LF.gof(fit.iso_line)


#--------------------------------- mark correlations functions-------------------------


env.ppp.line.ang<-envelope(ppp.line.ang,markcorr,nsim=999,nrank=49,savefuns=TRUE)
env.ppp.line.len<-envelope(ppp.line.len,markcorr,nsim=999,nrank=49,savefuns=TRUE)

plot(env.ppp.line.ang,legend=FALSE,xaxs="i",yaxs="i",font.lab=2,main="MCF segemnt angles",xlab="Length (cm)",ylab="MCF")
plot(env.ppp.line.len,legend=FALSE,xaxs="i",yaxs="i",font.lab=2,main="MCF segment lengths",xlab="Length (cm)",ylab="MCF")

LF.gof(env.ppp.line.ang)#$p#,ginterval=c(0,50)
LF.gof(env.ppp.line.len)#$p#,ginterval=c(0,50)

plot(env.ppp.line.len,legend=FALSE,xaxs="i",yaxs="i",font.lab=2,main="MCF lengths",xlab="Length (cm)",ylab="MCF",ylim=c(0.5,1.5))
plot(markcorr(ppp.line.len),lty=2,add=TRUE)

plot(env.ppp.line.ang,legend=FALSE,xaxs="i",yaxs="i",font.lab=2,main="MCF angles",xlab="Length (cm)",ylab="MCF",ylim=c(0.5,1.5))
plot(markcorr(ppp.line.ang),lty=2,add=TRUE)


#Spatial correlation----------------------

##autocorrelation corrected Pearson's correlation####
###autocorrelation test
library(spdep)
library(sf)
library(tmap)
morantest<-function(den1){
  ##dataframe without NA
  x_coords <- den1$xcol 
  y_coords <- den1$yrow
  value1 <- den1$v
  grid <- expand.grid(x = x_coords, y = y_coords)
  df1 <- data.frame(
    x = grid$x,
    y = grid$y,
    value = as.vector(value1)
  )
  
  df.clean<-na.omit(df1) ##remove NA
  ##moran test
  data_sf <- st_as_sf(df.clean, coords = c("x", "y"))
  coords <- st_coordinates(data_sf)
  knn <- knearneigh(coords, k = 5)
  nb <- knn2nb(knn)
  listw <- nb2listw(nb, style = "W")
  moran_test <- moran.test(df.clean$value, listw = listw)##Moran I statistics: >0 positive autocorrelation 
  print(moran_test)
  
}


###clifford corrected Pearson's correlation
library(SpatialPack)
clf.corr<-function(den1,den2)
{
  set1<-den1$v
  set2<-den2$v
  coords <- expand.grid(x = 1:nrow(set1), y = 1:ncol(set2))
  vec1 <- c(set1)
  vec2 <- c(set2)
  spatial_test <- modified.ttest(vec1, vec2, coords = coords)
  print(spatial_test)
}

#Overlapping extent: Schoener's D & Hellinger's I tests#### 
D.I <- function(den1, den2) {
  p1<-den1$v
  p2<-den2$v
  if (any(p1 < 0, na.rm = TRUE) || any(p2 < 0, na.rm = TRUE)) {
    p1[p1<0]<-NA
    p2[p2<0]<-NA
  }
  sum_p1 <- sum(p1, na.rm = TRUE)
  sum_p2 <- sum(p2, na.rm = TRUE)
  if (sum_p1 == 0 || sum_p2 == 0) {
    stop("values must not be 0/NA")
  }
  p1 <- p1 / sum_p1
  p2 <- p2 / sum_p2
  
  #remove NA
  mask <- !is.na(p1) & !is.na(p2)
  p1_masked <- p1[mask]
  p2_masked <- p2[mask]
  
  D <- 1 - 0.5 * sum(abs(p1_masked - p2_masked)) #Schoener's D
  
  I <- sum(sqrt(p1_masked * p2_masked)) #Hellinger's I
  
  return(list(D = D, I = I)) #[0,1] higher value=higher overlapping
}


# homogeneity tests------------------------

#HW test, LH*
homtest(unmark(ppp.line.ang),nsim=999)
sum(Gest(ppp.line.ang)$han-Gest(ppp.line.ang)$theo)

descdist(tr.sh.line)
descdist(tr.line.ang)

homtest(unmark(ppp.iso),nsim=999)
sum(Gest(ppp.iso)$han-Gest(ppp.iso)$theo)

homtest(unmark(ppp.con),nsim=999)
sum(Gest(ppp.con)$han-Gest(ppp.con)$theo)


####mclust of angle cohorts
mcangle<-Mclust(angle,G=2) #G=4
plot(mcangle, what = "BIC")
summary(mcangle,parameters=TRUE)
