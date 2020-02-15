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
str(Auto)
summary(Auto)
attach((Auto))
boxplot(Auto$cylinders , Auto$mpg )
cylinders =as.factor (cylinders )
plot(cylinders , mpg)
plot(cylinders , mpg , col ="red ")
plot(cylinders , mpg , col ="red", varwidth =T)
plot(cylinders , mpg , col ="red", varwidth =T,horizontal =T)
plot(cylinders , mpg , col ="red", varwidth =T, xlab=" cylinders ",
     ylab ="MPG ")

hist(mpg)
hist(mpg ,col =2, breaks =15)
pairs(Auto) #facet plot?
pairs(~ mpg + displacement + horsepower + weight + acceleration , Auto)

plot(horsepower ,mpg)
identify (horsepower ,mpg ,name)

dummyData=data.frame(rnorm (10000),rnorm (10000))

college = ISLR::College
fix(college)
head(college)
rownames(college)
#college = college[,-1]
fix(college)
summary(college)
pairs(college[,1:10])
plot(college$Private,college$Outstate)
college$Private= as.factor(college$Private)

college$Elite = as.factor(college$Top10perc>50)
summary(college$Elite)
plot(college$Elite, college$Outstate)

hist(college$Apps)
hist(college$perc.alumni, col=2)
hist(college$S.F.Ratio, col=3, breaks=10)
hist(college$Expend, breaks=100)

par(mfrow=c(1,1))
plot(college$Outstate, college$Grad.Rate)
# High tuition correlates to high graduation rate.
plot(college$Accept / college$Apps, college$S.F.Ratio)
# Colleges with low acceptance rate tend to have low S:F ratio.
plot(college$Top10perc, college$Grad.Rate)
# Colleges with the most students from top 10% perc don't necessarily have
# the highest graduation rate. Also, rate > 100 is erroneous!

auto = ISLR::Auto
auto = na.omit(auto)
summary(auto)
sapply(auto[, 1:7], range)
sapply(auto[, 1:7], mean)
sapply(auto[, 1:7], sd)

auto = auto[-(10:85),]
dim(auto)
sapply(auto[, 1:7], range)
sapply(auto[, 1:7], mean)
sapply(auto[, 1:7], sd)

# (e)
pairs(Auto)
plot(Auto$mpg, Auto$weight)
# Heavier weight correlates with lower mpg.
plot(Auto$mpg, Auto$cylinders)
# More cylinders, less mpg.
plot(Auto$mpg, Auto$year)
# Cars become more efficient over time.

Boston = MASS::Boston
?MASS::Boston
pairs(bost)

plot(Boston$age, Boston$crim)
# Older homes, more crime
plot(Boston$dis, Boston$crim)
# Closer to work-area, more crime
plot(Boston$rad, Boston$crim)
# Higher index of accessibility to radial highways, more crime
plot(Boston$tax, Boston$crim)
# Higher tax rate, more crime
plot(Boston$ptratio, Boston$crim)
# Higher pupil:teacher ratio, more crime

par(mfrow=c(1,3))
hist(Boston$crim[Boston$crim>1], breaks=25)
# most cities have low crime rates, but there is a long tail: 18 suburbs appear
# to have a crime rate > 20, reaching to above 80
hist(Boston$tax, breaks=25)
# there is a large divide between suburbs with low tax rates and a peak at 660-680
hist(Boston$ptratio, breaks=25)
# a skew towards high ratios, but no particularly high ratios

length(Boston$chas[Boston$chas==1])
dim(subset(Boston, chas == 1))
summary(Boston$ptratio)
subset(Boston, medv == min(Boston$medv))
t(subset(Boston, medv == min(Boston$medv)))

length(Boston$rm[Boston$rm>7])
dim(subset(Boston,rm>7))
summary(subset(Boston,rm>8))
summary(Boston)
pairs(subset(Boston,rm>8))

fix(Boston)
bos = Boston
write.csv(bos,file = 'c:/temp/boston.csv')
model = lm(medv~lstat,bos)
model
summary(model)

auto = ISLR::Auto
model = lm(mpg~horsepower, auto)
model
summary(model)
summary(model)$r.sq
predict(model, data.frame(horsepower=c(98)), interval="confidence")
predict(model, data.frame(horsepower=c(98)), interval="prediction")

par(mfrow=c(1,1))
plot(auto$horsepower, auto$mpg)
abline(model)

par(mfrow=c(2,2))
plot(model)

pairs(auto)

cor(auto)
auto_sb = subset(auto, select = -c(name))
names(auto_sb)
mat = cor(auto_sb)
mat>.4
mat[mat>.4]
which(mat>.4, arr.ind = T)

model = lm(mpg~.,auto_sb) #this is equivalent to following
model = lm(mpg~.-name,auto)
model
summary(model)
plot(model)
plot(predict(model), rstudent(model))

model2 = with(auto,lm(mpg~cylinders*displacement+displacement*weight))
summary(model2)

model3 = lm(mpg~log(weight)+sqrt(horsepower)+acceleration+I(acceleration^2), auto)
summary(model3)
par(mfrow=c(2,2))
plot(model3)

model4<-lm(log(mpg)~cylinders+displacement+horsepower+weight+acceleration+year+origin,data=auto)
summary(model4)

cs = ISLR::Carseats
model = lm(Sales~Price+Urban+US,cs)
summary(model)
plot(model)
str(cs$Urban)
as.double(cs$Urban)

model = lm(Sales~Price+as.double(Urban)+as.double(US),cs)
summary(model)

model = lm(Sales~Price+US,cs)
summary(model)
plot(model)

confint(model)

set.seed (1)
x=rnorm (100)
y=2*x+rnorm (100)

model = lm(y~x+0)#linear regression without inercept.
summary(model)
