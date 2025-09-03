library(spatstat)
library(ecespa)

# Load demo data
data_demo <- read.csv(file.choose())#demo_dataset
edge_demo <- read.csv(file.choose())#demo_edge

# Create window from edge
win.demo <- owin(poly=list(x=c(0,20,20,0), y=c(0,0,20,20)))

# Create point pattern
ppp.demo <- ppp(data_demo$X, data_demo$Y, window = win.demo, marks = data_demo$lenth)

# Plot points
plot(ppp.demo, main = "Demo trace fossil distribution")

# Density estimation
den.demo <- density(ppp.demo, sigma = 1)
plot(den.demo, main = "Demo density plot")
plot(ppp.demo, add = TRUE)

# Pair correlation function (PCF) with small number of simulations
plot(envelope(ppp.demo,pcf,nsim=99, nrank=4),legend=FALSE,xaxs="i",yaxs="i",font.lab=2,main="length pcf",ylim=c(0,5),xlab="Length (cm)",ylab="Midpoint PCF")
fit<-envelope(ppm(unmark(ppp.demo), trend=~1,statistic="pcf"),Gest,savefuns=TRUE,nsim=99,rank=4)
LF.gof(fit)

# Mark correlation funciton (MCF) with small number of simulations
env<-envelope(ppp.demo,markcorr,nsim=99,nrank=4,savefuns=TRUE)
plot(env,legend=FALSE,xaxs="i",yaxs="i",font.lab=2,main="MCF segemnt angles",xlab="Length (cm)",ylab="MCF")
