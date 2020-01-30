library(stats)
library(ISLR)
library(MASS)

Sys.getlocale()
# Sys.setlocale("LC_ALL", "C")
getwd()
setwd("C:/Users/sin17h/Downloads")
ls()

#ggplots drag and drop plot builder
#esquisse::esquisser()
#tidylog
#heatmaply

#w = read.csv('who.csv')

at = ISLR::Auto

x=rnorm (50)
y=x + rnorm(50, mean=50, sd=.1)
plot(rnorm(50, mean=50, sd=.1))

t = rnorm(50, mean=50, sd=.1)
t[t>50.1]
t[t<49.9]
sd(t)

3:7
x=seq(-pi ,pi ,length =50)
 y=x
 f=outer(x,y,function (x,y)cos(y)/(1+x^2))
 contour (x,y,f)
 contour (x,y,f,nlevels =45, add=T)
 fa=(f-t(f))/2
 contour (x,y,fa,nlevels =15)

 #replace contour with image func above to get heatmaps and with persp to get 3-dim plot
contour (x,y,fa)
 image(x,y,fa)
  persp(x,y,fa)
  persp(x,y,fa ,theta =30)
  persp(x,y,fa ,theta =30, phi =20)
  persp(x,y,fa ,theta =30, phi =70)
  persp(x,y,fa ,theta =30, phi =40)
 
A=matrix (1:16 ,4 ,4)
A[c(1,3) ,c(2,4) ]
A[-c(1,3) ,-c(1,3)]
dim(A[-c(1,3) ,-c(1,3)])
dim(A[-c(1,3) ,-c(1,3,4)])

Auto = ISLR::Auto
fix(Auto)
