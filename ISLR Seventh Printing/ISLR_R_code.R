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
y=2*x+rnorm (10)

par(mfrow=c(2,2))
model = lm(y~x+0)#linear regression without inercept.
summary(model)
confint(model)
plot(x,y) + abline(model)
plot(model)

model = lm(x~y+0)#linear regression without inercept.
summary(model)
plot(y,x) + abline(model)
plot(model)


set.seed(1)
x = rnorm(100)
y = 2*x
lm.fit = lm(y~x+0)
lm.fit2 = lm(x~y+0)
summary(lm.fit)
summary(lm.fit2)
confint(lm.fit)

set.seed(1)
x <- rnorm(100)
y <- -sample(x, 100)
length(unique(y))#it is 100 because the item once removed from the sample by sampling is not put back again. It is a parameter of this func  

# y <- -sample(x, 100, replace = TRUE)#by replacing the sample item back into the bin, the probability of same item getting sampled is not 0
# length(unique(y))#it is < 100 

sum(x^2)
sum(y^2)
lm.fit <- lm(y~x+0)
lm.fit2 <- lm(x~y+0)
summary(lm.fit)
summary(lm.fit2)

set.seed(1)
x = rnorm(100, 0, sqrt(1))
mean(x)
sd(x)

eps = rnorm(100, 0, sqrt(0.25))
y = -1 + .5*x + eps
length(y)
plot(x,y)

lm.fit = lm(y~x)
summary(lm.fit)
abline(lm.fit, lwd=3, col=2)
abline(-1, 0.5, lwd=3, col=3)
legend(-1, legend = c("model fit", "pop. regression"), col=2:3, lwd=3)


lm.fit2 = lm(y~x + I(x^2))
summary(lm.fit2)
abline(lm.fit2, lwd=3, col=2)

#less noise
set.seed(1)
x = rnorm(100, 0, sqrt(1))
mean(x)
sd(x)

eps = rnorm(100, 0, sqrt(0.025))
y = -1 + .5*x + eps
length(y)
plot(x,y)

lm.fit = lm(y~x)
summary(lm.fit)
abline(lm.fit, lwd=3, col=2)
abline(-1, 0.5, lwd=3, col=3)
legend(-1, legend = c("model fit", "pop. regression"), col=2:3, lwd=3)


lm.fit2 = lm(y~x + I(x^2))
summary(lm.fit2)
abline(lm.fit2, lwd=3, col=2)

#more noise
set.seed(1)
x = rnorm(100, 0, sqrt(1))
mean(x)
sd(x)

eps = rnorm(100, 0, sqrt(2.5))
y = -1 + .5*x + eps
length(y)
plot(x,y)

lm.fit = lm(y~x)
summary(lm.fit)
abline(lm.fit, lwd=3, col=2)
abline(-1, 0.5, lwd=3, col=3)
legend(-1, legend = c("model fit", "pop. regression"), col=2:3, lwd=3)


lm.fit2 = lm(y~x + I(x^2))
summary(lm.fit2)
abline(lm.fit2, lwd=3, col=2)

set.seed (1)
x1=runif (100)
x2 =0.5* x1+rnorm (100) /10
y=2+2* x1 +0.3* x2+rnorm (100)
plot(x1,x2)
cor(x1,x2)

lm.fit = lm(y~x1+x2)
summary(lm.fit)

lm.fit2 = lm(y~x1)
summary(lm.fit2)

lm.fit3 = lm(y~x2)
summary(lm.fit3)


x1 = c(x1, 0.1)
x2 = c(x2, 0.8)
y = c(y, 6)
lm.fit1 = lm(y~x1+x2)
summary(lm.fit1)
lm.fit2 = lm(y~x1)
summary(lm.fit2)
lm.fit3 = lm(y~x2)
summary(lm.fit3)

library(MASS)
summary(Boston)
Boston$chas <- factor(Boston$chas, labels = c("N","Y"))
summary(Boston)
attach(Boston)
lm.zn = lm(crim~zn)
summary(lm.zn) # yes
lm.indus = lm(crim~indus)
summary(lm.indus) # yes
lm.chas = lm(crim~chas) 
summary(lm.chas) # no
lm.nox = lm(crim~nox)
summary(lm.nox) # yes
lm.rm = lm(crim~rm)
summary(lm.rm) # yes
lm.age = lm(crim~age)
summary(lm.age) # yes
lm.dis = lm(crim~dis)
summary(lm.dis) # yes
lm.rad = lm(crim~rad)
summary(lm.rad) # yes
lm.tax = lm(crim~tax)
summary(lm.tax) # yes
lm.ptratio = lm(crim~ptratio)
summary(lm.ptratio) # yes
lm.black = lm(crim~black)
summary(lm.black) # yes
lm.lstat = lm(crim~lstat)
summary(lm.lstat) # yes
lm.medv = lm(crim~medv)
summary(lm.medv) # yes

lm.all = lm(crim~., data=MASS::Boston)
summary(lm.all)

x = c(coefficients(lm.zn)[2],
      coefficients(lm.indus)[2],
      coefficients(lm.chas)[2],
      coefficients(lm.nox)[2],
      coefficients(lm.rm)[2],
      coefficients(lm.age)[2],
      coefficients(lm.dis)[2],
      coefficients(lm.rad)[2],
      coefficients(lm.tax)[2],
      coefficients(lm.ptratio)[2],
      coefficients(lm.black)[2],
      coefficients(lm.lstat)[2],
      coefficients(lm.medv)[2])
y = coefficients(lm.all)[2:14]
plot(x, y)ta


#chapter 4 classification
sm = ISLR::Smarket
attach(ISLR::Smarket)
summary(sm)
head(sm,2)
str(sm)

pairs(cor(sm[sapply(sm , is.numeric)]))
a = cor(sm[sapply(sm , is.numeric)])
(a>-.3)&(a>.3) 
a[(a>-.3)&(a>.3)]=1
a[!((a>-.3)&(a>.3))] = 0
a#only year and volume shows a worthwhile correlation

plot(sm$Year, sm$Volume)
plot(as.factor(sm$Year), sm$Volume)
boxplot(as.factor(sm$Year), sm$Volume)#it doesn't work but above works. So use plot for boxplots?? 

#train the model on all the given data
plot(sm$Volume)
glm.fits=glm(Direction~Lag1+Lag2+Lag3+Lag4+Lag5+Volume, data=ISLR::Smarket, family=binomial)
summary(glm.fits)

#test the data using the same complete data set that was used for training. Not good practise
glm.probs =predict(glm.fits,type ="response")
glm.probs[1:10]
contrasts(Direction)

glm.pred=rep ("Down" ,1250)
glm.pred[glm.probs >.5]="Up"
table(glm.pred ,Direction)#confidence matrix

(507+145) /1250
mean(glm.pred == Direction)#fraction of days for which prediction was correct. Same as calculated above from confusion matrix.

train =(Year <2005)
Smarket.2005= ISLR::Smarket[!train,]#can be done using subset as well 
dim(subset(sm, !Year<2005))
dim(Smarket.2005)
Direction.2005= Direction[!train]

glm.fits=glm(Direction~Lag1+Lag2+Lag3+Lag4+Lag5+Volume ,
               data=sm ,family =binomial ,subset =train )
glm.probs =predict(glm.fits,Smarket.2005 , type="response")

glm.pred = rep("Down", 252)
glm.pred[glm.probs>.5] = "Up"
table(glm.pred,Direction.2005)#confusion matrix to calculate the error rate. Correct prediction rate in this case is (507+145)/1250.
mean(glm.pred == Direction.2005)#this can also be used to do same
mean(glm.pred != Direction.2005)#error rate is .519 which is worse than random guesses

#retrain the model on the training subset on only the predictors which have a low value(hence more significant). But why is Volume not included? 
glm.fits=glm(Direction~Lag1+Lag2, data=ISLR::Smarket, family=binomial, subset = train)
summary(glm.fits)

#now predict
glm.probs = predict(glm.fits, Smarket.2005, type="response")
glm.pred = rep("Down", 252)
glm.pred[glm.probs>.5] = "Up"
table(glm.pred,Direction.2005)
mean(glm.pred == Direction.2005)# Correct prediction rate has increased. it is now 56%

predict (glm.fits,newdata =data.frame(Lag1=c(1.2 ,1.5),Lag2=c(1.1 , -0.8) ),type ="response")


#linear discriminant analysis
library(MASS)
lda.fit=lda(Direction~Lag1+Lag2 ,data=ISLR::Smarket ,subset =train)
lda.fit
summary(lda.fit)
plot(lda.fit)

#predict
lda.pred=predict (lda.fit , Smarket.2005)
names(lda.pred)

lda.class = lda.pred$class
table(lda.class,Direction.2005)
mean(lda.class == Direction.2005)#correct prediction rate

sum(lda.pred$posterior [ ,1] >=.5)
sum(lda.pred$posterior [,1]<.5)

lda.pred$posterior[1:15 ,1]
lda.class[1:15]

#quadratic disciminant analysis
qda.fit=qda(Direction~Lag1+Lag2, data=ISLR::Smarket, subset=train)
qda.fit

qda.class =predict(qda.fit, Smarket.2005)$class
table(qda.class ,Direction.2005)

#now predict
qda.class = predict (qda.fit, Smarket.2005)$class
table(qda.class, Direction.2005)
mean(qda.class == Direction.2005)#this quadratic correct prediction rate is greater than linear LDA and Logistic regression.

#knn (k nearest neighbours) classification
train =(Year <2005)
library (class)
train.X=cbind(Lag1 ,Lag2)[train ,]
test.X=cbind (Lag1 ,Lag2)[!train ,]
train.Direction =Direction [train]

#predict
set.seed(1)
knn.pred=knn(train.X,test.X,train.Direction ,k=1)#k=1 the number of nearest neighbours to be used by the classifier
table(knn.pred ,Direction.2005)
mean(knn.pred == Direction.2005)#only 50% correct prediction rate

knn.pred=knn(train.X,test.X,train.Direction, k=3)#perhaps increasing k to 3 might improve prediction rate
table(knn.pred ,Direction.2005)
mean(knn.pred == Direction.2005)#slight improvement in correct prediction rate.

#in the end, for this classification excecise, out of lda,knn,qda, logistic regression, it is the qda which performs the best

Caravan = ISLR::Caravan
dim(Caravan)
attach(Caravan)
summary(Caravan$Purchase)

#scale or standardize the data. There is a subtle difference between them but standarization is preferred (sd of 1 and mean of 0).
standardized.X=scale(Caravan[,-86])
var(Caravan[,1])
var(standardized.X[,1])

#now predict using standardized data.
test =1:1000
train.X=standardized.X[-test ,]
test.X=standardized.X[test ,]
train.Y=Purchase [-test]
test.Y=Purchase [test]
set.seed (1)
knn.pred=knn(train.X,test.X,train.Y,k=1)
mean(test.Y!= knn.pred)
mean(test.Y!=" No")

knn.pred=knn(train.X,test.X,train.Y,k=3)
table(knn.pred ,test.Y)

knn.pred=knn(train.X,test.X,train.Y,k=5)
table(knn.pred ,test.Y)
#page 167 bookmark

#excercise 10 page 171
Weekly = ISLR::Weekly
summary(Weekly)
pairs(Weekly)

#Weekly[,lapply(Weekly, is.numeric)] # this wont work because we want the apply function to return strings(as in columns names) 
#which satisfy the condition. So we have to use the sapply equivalent.
num_weekly = Weekly[,sapply(Weekly, is.numeric)]
pairs(num_weekly)
cors = cor(num_weekly)
cors[(cors>.09 | cors < -.09 )] = 1
cors[!(cors>.09 | cors < -.09 )] = 0
cors
plot(as.factor(num_weekly$Year), num_weekly$Volume)

#the only correlation we see is between volume and year

#b.create a model using full data set
glm.model = glm(Direction~Lag1+Lag2+Lag3+Lag4+Lag5+Volume, data=ISLR::Weekly, family=binomial)
glm.model
summary(glm.model)

glm.pred = predict(glm.model,Weekly, type = "response")
glm.pred.Directions = rep("Down",nrow(Weekly))

glm.pred.Directions[glm.pred>.5] = "Up"

table(glm.pred.Directions, Weekly$Direction)#
mean(glm.pred.Directions == Weekly$Direction)# correct prediction rate is 56%

#total weeks the market goes up is 48+557 out of which we have predicted 557 times. 557/(557+48)= 92.1%.
#54/(430+54) = 11.2% of the time market goes up but model predicts wrong.

#this time divide the data set into train/test subsets
train_cond = Weekly$Year %in% c(1990:2008)
train_ds = Weekly[train_cond,]
test_ds =  Weekly[!train_cond,]

#train the logistic regression model
glm.model = glm(Direction~Lag2,data=train_ds,family = binomial)

#test the model
glm.pred = predict(glm.model, test_ds, type = "response")
glm.pred.Directions = rep("Down", nrow(test_ds))
glm.pred.Directions[glm.pred>.5] = "Up"
table(glm.pred.Directions, test_ds$Direction)
mean(glm.pred.Directions == test_ds$Direction)#correct prediction rate is 62.5%

#train the lda model
lda.model = lda(Direction~Lag2,data=Weekly ,subset = train_cond)
lda.pred = predict(lda.model,test_ds)

#verify the prediction rates
lda.pred.Directions = lda.pred$class
table(lda.pred.Directions, test_ds$Direction)
mean(lda.pred.Directions == test_ds$Direction)#correct prediction rate is 62.5%

#now train using qda. Note that train_cond is supplied
qda.model = qda(Direction~Lag2, data=Weekly, subset = train_cond)

qda.pred = predict(qda.model,newdata = test_ds)
qda.pred.class = qda.pred$class  
table(qda.pred.class,test_ds$Direction)
mean(qda.pred.class == test_ds$Direction)#58%. All the predictions were "Up"


#model using KNN=1
library (class)
train.X=as.matrix(Lag2[train_cond])#could have also written train_ds$Lag2 but the length of vectors was different! How so?
test.Y=as.matrix(Lag2[!train_cond])#could have also written test_ds$Lag2 
train.Direction = Direction[train_cond]

set.seed(1)
knn.pred=knn(train = train.X, test = test.Y, train.Direction, k = 1)
table(knn.pred, test_ds$Direction)
mean(knn.pred == test_ds$Direction)#prediction rate is 57%

#excercise 11 page 171
Auto = ISLR::Auto
Auto$mpg01[Auto$mpg > median(Auto$mpg)] = 1
Auto$mpg01[!(Auto$mpg > median(Auto$mpg))] = 0
Auto$mpg01

auto_cor = cor(Auto[sapply(Auto,is.numeric)])
pairs(Auto)
positive_cor_cond = (auto_cor>.3)
negative_cor_cond = (auto_cor< -.3)

positive_cor_cond | negative_cor_cond    

hist(Auto$year)#look at the data distribution to come up with a decent breakpoint
is.element(FALSE, is.na(Auto$year))# FALSE %in% is.na(Auto$year)

train_cond = Auto$year<76 # a better condition was (Auto$year %% 2 == 0) so that it sample the complete data set
train_ds = subset(Auto,subset = train_cond)# similar to Auto[train_cond,] but i prefer subset
test_ds = Auto[!train_cond,]

#11.d lets use .5 as the cutoff for correlation. 
positive_cor_cond = (auto_cor>.6)
negative_cor_cond = (auto_cor< -.6)
positive_cor_cond | negative_cor_cond    

#this gives us 4 variables with high cor
lda.model = lda(mpg01~cylinders+weight+displacement+horsepower,data = Auto, subset = train_cond)
lda.pred = predict(lda.model,test_ds)

table(lda.pred$class, test_ds$mpg01)
mean(lda.pred$class != test_ds$mpg01)# 10.8% error rate

#now lets change the condition used to split the dataset into train and test subsets
train_cond = (Auto$year %% 2 == 0) #better sampling condition
train_ds = subset(Auto,subset = train_cond)# similar to Auto[train_cond,] but i prefer subset
test_ds = Auto[!train_cond,]


#train and predict again
#this gives us 4 variables with high cor
lda.model = lda(mpg01~cylinders+weight+displacement+horsepower,data = Auto, subset = train_cond)
lda.pred = predict(lda.model,test_ds)

table(lda.pred$class, test_ds$mpg01)
mean(lda.pred$class != test_ds$mpg01)# 12.6% error rate. So the error rate increased. hmm!

#11.f logistic regression prediction. We are using the even year dataset splitting condition
glm.model = glm(mpg01~cylinders+weight+displacement+horsepower, data = Auto, subset = train_cond, family = binomial)

glm.probs = predict(glm.model,test_ds, type="response")
glm.pred = rep("0",length(glm.probs))
glm.pred[glm.probs>.5] = "1"

mean(glm.pred != test_ds$mpg01)#12% error rate


#11.g KNN
train_cond = (Auto$year %% 2 == 0) #better sampling condition
train_ds = subset(Auto, subset = train_cond, select = c(cylinders, weight, displacement, horsepower))# The same can be done using cbind(cylinders, weight, displacement, horsepower)[train,]
test_ds = subset(Auto, subset = !train_cond, select = c(cylinders, weight, displacement, horsepower))#cbind(cylinders, weight, displacement, horsepower)[test,]
label_ds_train = subset(Auto, subset = train_cond, select = c(mpg01))# this should have been equivalent to Auto$mpg01[train_cond] but it is not??
label_ds_test = subset(Auto, subset = !train_cond, select = c(mpg01))
#length(Auto$mpg01[train_cond])#210 length whereas above the length is 397. Why the difference
#label_ds = Auto$mpg01[train_cond]

set.seed(1)
knn.pred = knn(train = train_ds, test = test_ds, cl = label_ds_train$mpg01, k=1)
mean(knn.pred!=label_ds_test$mpg01)#error rate 15%

knn.pred = knn(train = train_ds, test = test_ds, cl = label_ds_train$mpg01, k=10)
mean(knn.pred!=label_ds_test$mpg01)#error rate 16.4%

knn.pred = knn(train = train_ds, test = test_ds, cl = label_ds_train$mpg01, k=100)
mean(knn.pred!=label_ds_test$mpg01)#error rate 14.3%


#excercise:12 page:172
Power = function(x,a){
  return(x^a)
}
Power(2,2)
x = 1:10000
plot(x,Power(x,2))
x = 1:10
plot(x, Power(x, 2),  log="x" )
plot(x, Power(x, 2),  log="y" )
plot(x, Power(x, 2),  log="xy", ylab="Log of y = x^2", xlab="Log of x",
     main="Log of x^2 versus Log of x")

PlotPower = function(x,y){
  plot(x,Power(x,y))  
}
PlotPower(1:10,4)

#page 191 chapter 5
set.seed(1)
Auto = ISLR::Auto
dim(Auto)#has 392 rows
row_indexes = sample(x=nrow(Auto),size=196)#this gives us random row indices that we can use into the df
Auto[row_indexes,]
#subset(x = Auto, subset = row_indexes) # since row_indexes is not a logical condition, we can't use it here. But we can use it as logical statement
#in the statement below??
lm.fit = lm(mpg~horsepower, data=Auto, subset = row_indexes)# same as Auto[row_indexes,]. So row_indexes is a conditions that can be used in subset 

#now predict the response/label/output MSE by  exluding the training rows
mean((Auto$mpg - predict (lm.fit ,Auto))[-row_indexes ]^2)#Since it is a regression and not a classification problem, we talk in terms of MSE and not error rate
#in this case, the MSE is 26.14

#now test using quadratic and cubic polynomial linear regressions (linear in terms of coefficients of the features)
lm.fit = lm(mpg~poly(horsepower,2), data=Auto, subset = row_indexes)
mean((Auto$mpg - predict (lm.fit ,Auto))[-row_indexes ]^2)#MSE for quadratic is 19.82

lm.fit = lm(mpg~poly(horsepower,3), data=Auto, subset = row_indexes)
mean((Auto$mpg - predict (lm.fit ,Auto))[-row_indexes ]^2)#MSE for quadratic is 19.78

#now we can resample again (by setting a different seed this time) to get a different training set and running the model again would give different MSE
set.seed(2)
row_indexes = sample(x=nrow(Auto),size=196)#this gives us random row indices that we can use into the df

lm.fit = lm(mpg~horsepower, data=Auto, subset = row_indexes)
mean((Auto$mpg - predict (lm.fit ,Auto))[-row_indexes ]^2)#MSE for quadratic is 23.29

lm.fit = lm(mpg~poly(horsepower,2), data=Auto, subset = row_indexes)
mean((Auto$mpg - predict (lm.fit ,Auto))[-row_indexes ]^2)#MSE for quadratic is 18.90

lm.fit = lm(mpg~poly(horsepower,3), data=Auto, subset = row_indexes)
mean((Auto$mpg - predict (lm.fit ,Auto))[-row_indexes ]^2)#MSE for quadratic is 19.25


#glm() function performs logistic regression when used with argument family="binomial". Without it, it performs linear regression
#just like the lm() function. GLM stands for generalized linear model whereas lm stands for linear model  
#cv.glm() function, part of boot package, gives cross validation results
library (boot)
glm.fit=glm(mpg~horsepower ,data=Auto)
cv.err =cv.glm(Auto ,glm.fit)
cv.err$delta

cv.error=rep (0,5)
for (i in 1:5){
   glm.fit=glm(mpg~poly(horsepower ,i),data=Auto)
   cv.error[i]=cv.glm (Auto, glm.fit)$delta [1]
  }
cv.error#decrease in ross validation error when we moved to quadratic but thereafter no decrease in error as we used higher polynomials


set.seed (17)
cv.error.10= rep (0 ,10)
for (i in 1:10) {
  glm.fit=glm(mpg~poly(horsepower ,i),data=Auto)
  cv.error.10[i]=cv.glm (Auto, glm.fit, K=10) $delta [1]
  }
cv.error.10

#page 194
