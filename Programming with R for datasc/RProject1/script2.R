#vector is homogenous in terms of its member types (like array in c++). Same in case of r matrix
#list is like a generic collection where members can have different types. Same in case of r df
library(ggplot2)
c(1,2,3) #vector
1:5 #sequence which is internally a vector?
? typeof
? sample

? rnorm #normal distribution using random no.s rnorm(5,mean=2, sd=1)
? rgamma #gamma distribution using random no.s rgamma(5,shape=2, rate=2)
? rbinom #binomial distribution using random no.s rbinom(5,size=100, prob=.3)
? rmultinom #multi nominal distribution using random no.s rmultinom(5,size=100, prob=c(.3,.2,.5))
? runif #uniform distribution using random no.s runif(5,min=1, max=2). Normal distribution 
#is gaussian distribution

? c #this is the concatenate function to produce a vector, not cat
? cat # catalogue function like print. Can also redirect output to a file
? paste # is for concatenating cahracter vectors. Like a '+' for string
? lapply #output is list. also tapply given below
? sapply #output is vector or matrix. sapply(df,class). also tapply given below
#data.frame(5:7, c(14, 12, 13))  note that what u specify in dataframe constructor is equal length columns. 
#Matrix can also be used to instantiate data.frame if the constructor of matrix is more useful for a use case.
? seq
? rep #replicates the values
? replicate # replicate(n=100, simmie.1()) # repeatedly call the function n times.
? rbind # this takes inputs which becomes rows
? cbind # this takes inputs which becomes columns
? matrix # this takes one long vector and then you supply shape to transform the vector into the matrix
#include the drop=FALSE argument in all your matrix manipulation code
#this dimension reduction problem is also prevalent in dataframes
#the other problem in dataframes is auto conversion of string to factors. For that use stringsAsFactors=FALSE
#or as.is=TRUE
? mode ;mode("Name");mode(2)#gives you storage mode. character,numeric, etc.
? length
? dim
? nrow
? ncol
? dimnames # retreive dimension names
? rownames # retreive dimension names
? colnames # retreive dimension names
? getwd
? setwd
? read.table #in columns of logical, integer, numeric and complex values that R interprets a blank space as a missing value
? read.csv # they are views on read.table 
? read.csv2 # they are views on read.table 
? readLines
? scan # read into a vector
? file #close
#use of [] and [[]] ???
? write.table
? write.csv # they are views on write.table
? write.csv2 # they are views on write.table 
? readLines
? writeLines
? sink #redirects output to a file instead of command prompt
? dump #serializes the R objects to a file. But this is ascii serialization.
? source # deserialize the r objects from a file. But this is ascii de-serialization.
? dput #dput and dget also do something similar as dump and source but with some differences. 
#But this is ascii serialization.
? dget #dput and dget also do something similar as dump and source but with some differences. 
#But this is ascii serialization.
? save # serialize in r-binary format
? load # de-serialize in r-binary format
? ls #list r objects from workspace
? rm # remove r objects from workspace
? library(RODBC) # sql data connections. odbc library is more efficient that rodbc
? merge #merge dataframes..sort of like inner join
? names #names(df) gets col names
? class #class of df column, i.e., numeric or factor, etc.
? str #str(df) gets object structure 
?fix
#logical indexing df[df$temp>90,c("temp","rain")]. can also use subset for this
?subset #subset(df, temp>80, select=c(temp,rain))
#subset(df, temp>80, select = -temp)
#subset(df, select=temp:rain)
#subset(df, !is.na(temp))
? summary #descriptive stats
? grep
? grepl
? sub
? table
? log  
? tapply #tapply(df$height, df$gender, mean). could also be done by 
#aggregate(height~gender, df, mean) AND by(df$height, df$gender, mean)

#tapply(df$height, list(df$gender, df$tmt), mean). could also be done by aggregate
#(height~gender+tmt, df, mean) AND by(df$height, list(df$gender, df$tmt), mean).
#but the way the output is returned is different for tapply, aggregate and by. But for large 
#datasets, tapply is the fastest and aggregate the slowest. Note the difference between
#these 3 functions with table. Also similar is xtabs 

? aggregate(height~gender, df, mean) #see tapply
? by(df$height, df$gender, mean)  #see tapply
? attach #attached datasets so that you do have to refer them by name. But it has it's issues 
#around hiding names. So avoid
? detach
? with #temporary attach
? searchpaths
? table #table(airquality$Ozone>80, airquality$Month). compare it with xtabs, tapply, aggregate and by
? xtabs # compare it with table, tapply, aggregate and by
? ftable
? dplyr  
? par # par(mfrow=c(3,3))
?example #runs the code in help section for the specified keyword/operator
?ifelse #one line vectorized conditional statement. It can be nested as well:ifelse inside of ifelse 

example(seq)  
example("&")
example("hist")
example(persp)

x <- rnorm(10000)#normal distribution(different from Uniform/Gaussian)
fix(x)
View(x)
head(x)
mean(x)
var(x)
#p<-ggplot2::ggplo
hist(x)
x + ggplot2::geom_histogram()
ggplot2::geom_histogram(x)
t(x) + ggplot2::geom_histogram()

mean(exp(rnorm(100000)))


exp(0 - 1 / 2)
exp(0 + 1 / 2)
ls()
rm()

a <- 1
b <- 3
c <- 1

my.vector <- c((-b + sqrt(b ^ 2 - 4 * a * c)) / 2 * a, (-b - sqrt(b ^ 2 - 4 * a * c)) / 2 * a)
my.vector
c(-0.4, -2.6) / my.vector * 100
?? "kmeans+"

x <- rnorm(100, mean = .5, sd = .3)
mean(x)
sd(x)
hist(x)

set.seed(1234)
x <- rnorm(100, mean = .5, sd = .3)
mean(x)
sd(x)
hist(x,axes=FALSE)
hist(x, axes = FALSE, ylab="")
axis(4)
axis(1)

set.seed(1)
x <- rnorm(100)
head(x)

?tail
help(tail)
? "tail"

x <- 5:6
x

ls()

a <- 1:10
a<- seq(1,10,2)
a<-(1:5)*2 -1
a <- c(1:20)[c(T, F)]
a

A <- seq(1, by = 2, len = 5)
B <- mean(A)
X <- seq(2, by = 2, len = 5)
Z <- A + X

as.complex(-3)
as.complex(-.3)
as.complex(-.00003)

paste("R session", 1)
data.frame(5:7, c(14, 12, 13)) + 2
data.frame(5:7,2:4, c(14,12,13)) + 2
data.frame(5:7, c("rr", "ss", "ff")) + 2

myfunc <- function(v) {
    return (v - mean(v)) / sd(v)
}

myfunc(1:5)
a <- rnorm(1000)
mean(a)
var(a)
s <- myfunc(a)

fofx <- function(x) {
    3 * sin(x / 2) + x
}

?plot.function
example("plot.function")
# there are 2 ways to plot a function. One is the evaluate the function first and then plot it's
#results. Other is to pass the function and the parameters it expects to the plot function and
#let the plot evaluate the func
g <- function(t) { return ((t^2+1)^0.5) }  # define g()
x <- seq(0,5,length=10000)  # x = [0.0004, 0.0008, 0.0012,..., 5]
y <- g(x)  # y = [g(0.0004), g(0.0008), g(0.0012), ..., g(5)]
plot(x,y,type="l")

plot.function(g,from = 0, to=5, n=50000)#n is the number of steps between 'from' and 'to'


plot(fofx, from = -7, to = 7)
plot(fofx, from = -7, to = 7, n=5)#the plot is not smooth
plot(fofx, from = -7, to = 7, n=15)#the plot is now smooth as number of steps is increased
ggplot(data = data.frame(x = 0), mapping=aes(x=x)) +
  stat_function(fun=fofx, color="blue") + xlim(-7,7)

g1 <- function(x) return(sin(x))
g2 <- function(x) return(sqrt(x^2+1))
g3 <- function(x) return(2*x-1)
plot(c(0,1),c(-1,1.5))  
# prepare the graph, specifying X and Y ranges> 
for (f in c(g1,g2,g3)) plot(f,0,1,add=T)  # add plot to existing graph

?nlm #gives you the function minimum. You can also plot the function to visually look the min
#note that the way you get the minimum of a function is to set the derivative of the function to 0(the
#point of 0 slope). If there are more than 1 vars in the func, take partial derivatives wrt them and 
#set them to 0.
func<-function(x) return(x^2-sin(x))
nlm(func,8)#minimum value was found to be approximately −0.23, occurring at x = 0.45
plot.function(func,from = -100,to = 100)
plot.function(func,from = -100,to = 100, n=100000)#n is the number of steps between from and to. 
#Default is 101 steps i think. Increasing the n would make the plot smoother if it is wavy.
dx2x <- deriv(~ x^2, "x") ; dx2x

#Method1. you can also plot a 3d surface. Uses expand.grid
library(lattice)
a <- 1:10
b <- 1:15
eg <- expand.grid(x=a,y=b)
eg$z <- eg$x^2 + eg$x * eg$y
wireframe(z ~ x+y, eg)
wireframe(z ~ x+y, eg, shade=T)
cloud(z ~ x+y, eg)

#method 2. another way 3d surface plot it can be plotted. Without expand.grid with nested loops to fill in the grid
# Grid over which we will calculate J
theta0_vals <- seq(-10, 10, length.out=100)
theta1_vals <- seq(-2, 4, length.out=100)

# initialize J_vals to a matrix of 0's
J_vals <- matrix(data=0,nrow=length(theta0_vals), ncol=length(theta1_vals))

#this is the cost func
J <- function(X, y, theta) {
  m = nrow(X)#length(y) or length(theta) would have been same
  (2*m)^-1 * sum((H(X,theta) - y)^2)
}

# Fill out J_vals
for (i in 1:length(theta0_vals)) {
  for (j in 1:length(theta1_vals)) {
    J_vals[i,j] <- J(X, y, c(theta0_vals[i], theta1_vals[j]))
  }
}
wireframe(J_vals, drape=T, col.regions=rainbow(100))

#method 3: 3d surface plot. does not use expand.grid or nested loop but uses outer to populate the grid
fdejong <- function (x, y) {
  return (x^2 + y^2)
}
x <- seq(-10, 10, length= 30)
y <- x
z <- outer(x, y, fdejong)
z[is.na(z)] <- 1
require(lattice)
wireframe(z, drape=T, col.regions=rainbow(100))

#method 4: 3d surface plot that you can move around(interactive). 
# Grid over which we will calculate J
theta0_vals <- seq(-10, 10, length.out=100)
theta1_vals <- seq(-2, 4, length.out=100)
# initialize J_vals to a matrix of 0's
J_vals <- matrix(0,length(theta0_vals), length(theta1_vals))
# Fill out J_vals
for (i in 1:length(theta0_vals)) {
  for (j in 1:length(theta1_vals)) {
    J_vals[i,j] <- J(X, y, c(theta0_vals[i], theta1_vals[j]))
  }
}
#interactive 3D surface plot
#install.packages("rgl")
library(rgl) 
open3d()

nbcol = 100
color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
J_vals_col  = cut(J_vals, nbcol)

persp3d(theta0_vals, theta1_vals, J_vals,col = color[J_vals_col],
        xlab=expression(theta_0),ylab=expression(theta_1),
        zlab="Cost",main = "Gradient Descent")
points3d(theta_history[, 1], theta_history[, 2], J_history+10, 
         col="red",size=3.5)
lines3d(theta_history[, 1], theta_history[, 2], J_history+10, col="red")

#method 5: this is for a 3d interactive scatter plot(not surface)
library(rgl)

open3d()
plot3d(X[,2],X[,3],y, 
       xlab= "sq-ft of room", ylab="#bedroom", zlab="price", col="blue",  
       type="s",size=1.5, main="Result of Gradient Descent")

xx <- seq(0,5000,length.out=25)
yy <- seq(1,5,length.out = 25)
zz <- matrix(0,length(xx),length(yy))

for (i in 1:length(xx))
  for (j in 1:length(yy))
    zz[i,j] <- cbind(1, (xx[i]-mu[1])/sigma[1],(yy[j]-mu[2])/sigma[2]) %*% theta

#MATLAB Like plane
nbcol = 100
color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
zcol  = cut(zz, nbcol)
persp3d(xx,yy,zz, add = TRUE, col=color[zcol],alpha=.6)



matrix(c(5, 4, 3, 2, 1, 0) + 2, nrow = 2) < 5
sin
? cat #concatenate and print
? density
? hist
? "if"
? plot
? summary

summary(a)
sd(a)
??skewed

my.display <- function(x, display = FALSE, type, prob) {
    cat("Summary of input: \n")
    if (display==TRUE) {
        if (type == "") {

        } else {

        }
    } else {

    }
    return(summary(x))
}
my.display(a)
set.seed(1234)
my.data <- rnorm(200)
my.display(my.data)

set.seed(1786)
data.exercise.3.1 <- exp(matrix(rnorm(2000), nrow = 100))
index1.temp <- sample(1:100, 10)
index2.temp <- sample(1:20, 10)
for (i in 1:10) {
    data.exercise.3.1[index1.temp[i], index2.temp[i]] <- -1
}

my.data <- data.exercise.3.1

for (i in row(my.data)) {
    if (any(my.data[1, ] < 0) == TRUE) a else b
    #my.data[1,my.data[1,]<0]
    }
? any

k <- 10
y <- matrix(rnorm(k ^ 2), nrow = k)
z <- 0 * y
head(y)
#loop:
time1 <- as.numeric(Sys.time())
for (i in 1:k) {
    #loop:
    for (j in 1:k) {
        z[i, j] <- y[i, j] ^ 2
    }
}

time2 <- as.numeric(Sys.time())
# using object form in R:
time3 <- as.numeric(Sys.time())
# using object form in R:
z <- y ^ 2
time4 <- as.numeric(Sys.time())
# run time increase factor:
(time2 - time1) / (time4 - time3)

? numeric
my.dimensions <- c(10, 20, 50, 100, 200, 500, 800, 1000)
my.runtime.factors <- numeric(8)

r <- 8
k <- my.dimensions[r]
y <- matrix(rnorm(k ^ 2), nrow = k)
z <- 0 * y
time1 <- as.numeric(Sys.time())
#loop:
for (i in 1:k) {
    for (j in 1:k) {
        z[i, j] <- y[i, j] ^ 2
    }
}
time2 <- as.numeric(Sys.time())
time3 <- as.numeric(Sys.time())
# using object form in R:z<-y^2
time4 <- as.numeric(Sys.time())
# run time increase factor: 
my.runtime.factors[r] <- (time2 - time1) / (time4 - time3)

plot(my.dimensions ^ 2, my.runtime.factors, log = "xy", xlab = "Number of operations")


k1 <- 10
k2 <- 100000
my.data <- as.data.frame(matrix(rnorm(k1 * k2), nrow = k1))

mean1 <- numeric(k2)
mean2 <- numeric(k2)

time1 <- as.numeric(Sys.time())
for (i in 1:k2) {
    mean1[i] <- mean(my.data[, i])
}

time2 <- as.numeric(Sys.time())
time3 <- as.numeric(Sys.time())
mean2 <- sapply(my.data, mean)
time4 <- as.numeric(Sys.time())
(time2 - time1) / (time4 - time3)

x <- matrix(1:12, 4)
x[cbind(c(1, 3, 2), c(3, 3, 2))]
x[c(1, 3, 2), c(3, 3, 2)]


row <- matrix(rep(1:100, 100), nrow = 100)
column <- matrix(rep(1:100, 100), nrow = 100, byrow = T)
A <- 3 * column ^ 3 / (1 + row * column)

dim(A)
sum(A)

sum(sum(A[, 1]))
sum(A[row <= column]) # A, row and column are 3 different matrices. c how 
#other matrices themselves are used to index into a separate matrix(that too
#in a condition!)

x <- 1:3
y <- seq(4, 8)
z <- rep(9:10, 1)

set.seed(9852)
my.data <- list()
for (i in 1:100) {
    my.data[[i]] <- matrix(rnorm(16), nrow = 4)
}

my.index <- list()
for (i in 1:100) {
    #my.index[[i]] <- my.data[[1]] < 0
    my.index[[i]] <- (my.data[[i]] < 0)
}

my.negatives <- matrix(rep(0, 16), nrow = 4)
for (i in 1:100) {
    my.negatives <- my.negatives + my.index[[i]]
}

my.negatives 
sum(my.negatives)

my.negative.values <- numeric(0)
for (i in 1:100) {
    my.negative.values <- c(my.negative.values, my.data[[i]][my.index[[i]]])
}
length(my.negative.values)
summary(my.negative.values)

my.data.frame <- read.table("data.exercise5.1.dat", header = TRUE, skip = 1)
my.data.frame <- read.table("data.exercise5.2.dat", skip = 2, header = TRUE, 
                            sep = ";", dec = ",")

my.data.frame <- read.csv("Exercise 5.3.csv", na.strings = c("","NA"), skip = 2)[, -1]
is.na(my.data.frame[4, 8])
is.na(my.data.frame[4, 7])
is.na(my.data.frame[3, 2])

my.data.frame <- read.table("Exercise 5.4a.txt", header = TRUE, skip = 2)
dimnames(my.data.frame)[[2]] <- c("", "", "")

l <- levels(unlist(read.table("Exercise 5.4a.txt", header = F, skip = 1, nrows = 1)))
dimnames(my.data.frame)[[2]] <- levels(unlist(read.table("Exercise 5.4a.txt", header = F,
                                                         skip = 1, nrows = 1)))[2:length(l)]

f1 <- file("Exercise 5.4a.txt", open = "r")
my.names <- scan(f1, what = "", nlines = 1, skip = 1)
my.data <- read.table(f1)
close(f1)

my.names <- paste(my.names[c(1, 3, 5)], my.names[c(2, 4, 6)])
names(my.data) <- my.names

my.data <- list()
my.names <- character(2)
f1 <- file("Exercise 5.4b.txt", open = "r")
my.names[1] <- scan(f1, what = "", nlines = 1, skip = 1)
my.data[[1]] <- scan(f1, nlines = 1)
my.names[2] <- scan(f1, what = "", nlines = 1)
my.data[[2]] <- matrix(scan(f1), nrow = 4, byrow = T)
close(f1)
names(my.data) <- my.names
my.data

set.seed(9007)
my.data <- data.frame(x = rnorm(10), y = rnorm(10) + 5, z = rchisq(10, 1))
additional.data <- data.frame(x = rnorm(3), y = rnorm(3) + 5, z = rchisq(3, 1))
write.table(my.data, "Exercise 6.1.txt", row.names = FALSE, col.names = FALSE)
write.table(additional.data, "Exercise 6.1.txt", row.names = FALSE, 
            col.names = FALSE, append = T)

set.seed(45)
my.data <- data.frame(x = rnorm(10), y = rnorm(10), z = rnorm(10))
write.csv(my.data, "Exercise 6.2.csv")

my.data<-data.frame(a=LETTERS[1:5],b=1:5)
write.table(my.data,file="Exercise 6.3a.txt", sep=";",row.names=FALSE)

my.text<-"TITLE extra line\n2 3 5 7\n11 13 17 \nOne more line"
writeLines(my.text,con="Exercise 6.3b.txt")

ls()

set.seed(45)
my.data<-data.frame(x=rnorm(10),y=rnorm(10),z=rnorm(10))
save(my.data,file="Exercise 6.4.Rdata")
rm(my.data)
my.data
load("Exercise 6.4.Rdata")

write.csv(mtcars,file = "mtcars.txt")
write.table(mtcars,file = "mtcars1.txt",sep = ",")

set.seed(9007)
my.data<-data.frame(x=rnorm(10),y=rnorm(10)+5,z=rchisq(10,1))
write.table(round(my.data,digits = 2),"Assignment 6a.txt",row.names=FALSE) 
write.table(format(my.data,digits = 20),"Assignment 6b.txt",row.names=FALSE)
#options(digits=20)

library(RODBC)#use the new odbc driver as it is faster
connStr <- paste(
  "Server=msedxeus.database.windows.net",
  "Database=DAT209x01",
  "uid=RLogin",
  "pwd=P@ssw0rd",
  "Driver={SQL Server}",
  sep=";"
)
conn <- odbcDriverConnect(connStr)
sqlTables(conn, tableType = "TABLE")
sqlColumns(conn,"bi.sentiment")[c("COLUMN_NAME","TYPE_NAME")]
sqlQuery(conn,"SELECT COUNT(*) FROM bi.sentiment")

my.data.frame<- sqlQuery(conn,
                         "SELECT AVG(Score),Date
                         FROM bi.sentiment
                         WHERE State='WA'
                         GROUP BY Date"
)
names(my.data.frame)<-c("Average Score","Date")
my.data.frame
my.data.frame <- sqlQuery(conn,
                         "SELECT SUM(Revenue), SUM(Units), ProductID
                         FROM bi.salesFact
                         WHERE Date > '2013-12-31' AND Date < '2015-01-01'
                         GROUP BY ProductID"
)
names(my.data.frame)<-c("SUM(Revenue)","SUM(Units)","ProductID") 
order(my.data.frame["ProductID"],decreasing = T)
my.data.frame$ProductID[order(my.data.frame$"SUM(Units)",decreasing=TRUE)][1:5]
my.data.frame$ProductID[order(my.data.frame$"SUM(Revenue)",decreasing=TRUE)][1:5]
my.data.frame[order(my.data.frame$"SUM(Revenue)",decreasing=TRUE),]

data.frame.x<-data.frame(names=c("Gretha","Robert","John","Heather"),
                         age=c(30,18,25,70),
                         nickname=c("Quicksilver","The Man","Nifty","Starlight"))

data.frame.y<-data.frame("Person_name"=c("William","Nancy","Charlotte","Henry"),
                         age=c(15,75,32,51),
                         "pet_dog"=c("King","Whity","Captain Vom","Doggie"))

data.frame.z<-merge(data.frame.x,data.frame.y)
data.frame.z<-merge(data.frame.y,data.frame.x,
                    by.x=c("Person_name","age"),
                    by.y=c("names","age"),all=TRUE)#outer join as all=TRUE

#just keep in mind complete rows are being removed here just like a SQL where clause
s<-subset(airquality,!is.na(Ozone) & !is.na(Solar.R))
colMeans(s)
mean(airquality$Solar.R, na.rm = T)


names(iris)
levels(iris$Species)
median(iris$Sepal.Length)
setosa.data<-subset(iris, 
                    Species == "setosa" & Sepal.Length<median(Sepal.Length),
                    select = -Species) 
summary(setosa.data)

my.text<-"Over the last decade, bluetongue virus have spread northwards from the mediterranean area. Initially this was ascribed to climate changes, but it has since been realized that a major contributing factor has been new transmitting vectors, culicoides obsoletus and culicoides pulicaris, which have the ability to aquire and transmit the disease. Recently, schmallenberg virus has emerged in northern europe, transmitted by biting midges as well."
my.lowercase<-c("bluetongue","culicoides","europe","mediterranean",
                "northern","schmallenberg")

my.uppercase<-c("Bluetongue","Culicoides","Europe","Mediterranean",
                "Northern","Schmallenberg")
my.new.text<-my.text
for(i in 1:length(my.lowercase)){
  my.new.text<-gsub(my.lowercase[i],my.uppercase[i],my.new.text)
}

Set.seed(885)
my.posixct<-as.POSIXct(sample((60*60*24*365*50):(60*60*24*365*55),20), 
                       origin = as.Date("1960-01-01"))
my.posixct2 <- my.posixct+9010
head(data.frame(my.posixct,my.posixct2))

set.seed(449)
your.dates<-as.Date(sample(18000:20000,20), origin = "1960-01-01")

your.days<-c(julian(your.dates,origin=as.Date("1960-01-01"))) 
             
set.seed(119)
my.days<-sample(18000:20000,20)
library(chron)
my.days.structure<-month.day.year(my.days,origin=c(1,1,1960)) 
my.dates<-as.Date(my.days, origin = "1960-01-01") 
my.date.info<-data.frame(Weekday=weekdays(my.dates),my.days.structure) 

log.airquality <- log(airquality)

#3 different ways to do same
tapply(df$height, df$gender, mean)
aggregate(height~gender, df, mean)
by(df$height, df$gender, mean)

#presidents is timeseries, not df
presidents
cycle(presidents) 
tapply(presidents,cycle(presidents),mean,na.rm=T)
time(presidents)
frequency(presidents)
deltat(presidents)

searchpaths()

#cut converts numeric to factor by dividing the input range into interval.
#sort of like bin'ing for a histogram
tapply(airquality$Solar.R , cut(airquality$Wind,10), mean,na.rm=TRUE)

a <- table(airquality$Ozone>80, airquality$Month)
a
a2 <- addmargins(a)
a2

#return proportions wrt the given dimension
a3 <-prop.table(a,2)
a3
# a3 <-prop.table(a,1)
# a3

addmargins(a3)
round(addmargins(a3)*100)

xtabs(Freq ~ Gender + Admit + Dept, UCBAdmissions)
admissions <- xtabs(Freq ~ Gender + Admit + Dept, UCBAdmissions)
ftable(xtabs(Freq ~ Gender + Admit + Dept, UCBAdmissions))
margin.table(xtabs(Freq ~ Gender + Admit + Dept, UCBAdmissions), 1:2)
prop.table(margin.table(xtabs(Freq ~ Gender + Admit + Dept, UCBAdmissions), 2:3),1)
prop.table(margin.table(xtabs(Freq ~ Gender + Admit + Dept, UCBAdmissions), c(1:3)),1)
ftable(prop.table(margin.table(xtabs(Freq ~ Gender + Admit + Dept, UCBAdmissions), c(1:3)),1))

depA <- admissions[,,1]
ftable(depA)
prop.table(depA,1)

summary(swiss) # tapply(airquality$Solar.R , cut(airquality$Wind,10), mean,na.rm=TRUE)
my.cut2<-cut(swiss$Agriculture,breaks=10*(0:10))
my.cut3<-cut(swiss$Catholic,breaks=10*(0:10))

tapply(swiss$Fertility,list(my.cut2,my.cut3),mean, na.rm=TRUE)


#You would like to show how many times in a month the temperature goes below 65 degree. 
with(airquality, table(Month, Temp < 65)) 
table(LowTemp = airquality$Temp < 65, airquality$Month) 

#You would like to show the percentage in each month the temperature goes above 80 degree.
prop.table(with(airquality, table(Month, Temp > 80)),1) 
prop.table(with(airquality, table(Temp > 80, Month)),2) 

my.data<-data.frame(Treatment=c(rep("A",4),rep("B",4)),
                    Stone=rep(rep(c("Small","Large"),c(2,2)),2),
                    Success=rep(c(1,0),4),
                    Count=c(81,6,192,71,234,36,55,25))
(my.data)
my.table<-xtabs(Count~Treatment+Success+Stone,data=my.data) 
ftable(my.table)

#You want to compare the success of treatment A vs treatment B using margin.table() and prop.table().
prop.table(margin.table(my.table, 1:2),1) 


#success rate of each treatment for the small stone type.
prop.table(my.table[,,2],1) 

#success rate of each treatment for the large stone type
prop.table(my.table[,,1],1) 

#success rate of each stone types for Treatment A.
prop.table(my.table[1,,],2) 

#success rate of each stone types for Treatment A.
prop.table(my.table[2,,],2) 

#Are the stones randomly allocated to the treatments? Let us look at the marginal table for Stone and Treatment.
prop.table(margin.table(my.table, c(3,1)),1) 
prop.table(margin.table(my.table, c(1,3)),2) 


plot(rnorm(10000))

par(mfrow=c(2,2))
set.seed(779)

# plot 4 histograms of 12 simulated standard normally distributed numbers. Add 
#the density of the standard normal to the plot with the curve() function.
#as we increase the sample size, the distribution is visibly more normal. Not to
#say that the low sample sizes did not contains numbers from a normal distribution.
for(i in 1:4){
  hist(rnorm(12), probability=TRUE,main=paste("Histogram",i))
  curve(dnorm,add=TRUE,col="red",lwd=3)
}

for(i in 1:4){
  hist(rnorm(25), probability=TRUE,main=paste("Histogram",i))
  curve(dnorm,add=TRUE,col="red",lwd=3)
}

for(i in 1:4){
  hist(rnorm(10000), probability=TRUE,main=paste("Histogram",i))
  curve(dnorm,add=TRUE,col="red",lwd=3)
}


#integration using random numbers. What we are doing here is trying to integrate a
#function between 0 and 2*pi. So we do a monte-carlo simulation to generate uniform(gaussian)
#distribution between the given integration limits and then find proportion of numbers 
#from that distribution that fall under the curve we are trying to integrate. As we 
#increase the number of points in the distribution, we get closer to
#expected answer.
#what i do not understand here is the that why was mean used. We wanted 
#the (count of numbers below the curve)/N where as what we are doing 
#is (sum of numbers below the curve)/N
#that is because of the way mean is being used...it is being used on a logical
#condtion(1 or 0). And is thus summing up the 1 and 0 s which will actually give
#the count of 1's. And thus the mean of that is what we want. Example case 
#showing the similarities is given below
N <- 100
x <- runif(N,0,2*pi)
y <- runif(N,0,8)

#do not plot as the points in distribution increase as it would just hang.
plot(x,y)
plot(x,exp(2*cos(x-pi)))
#asas <- c(1,2,3)
#mean(asas>1)
phat.under <- mean(y<exp(2*cos(x-pi)))
phat.under * 16 * pi

#asas <- c(1,2,3)
#length(asas[asas>1])/N
phat.under <- length(y[y<exp(2*cos(x-pi))])/N
phat.under * 16 * pi

#integration using numerical integrate function
f <- function(x){exp(2*cos(x-pi))}
integrate(f,0,2*pi)


doone <- function(){
  x <- rbinom(1,size = 50,prob = 1/6)  
  p <- x/50
  p
}

p.sim <- replicate(1000,doone())
hist(p.sim,breaks = 15)


x <- rnorm(1000)
plot(density(x), xlim = c(-8,16))

# mean shifted 8 units to the right..but the curve is similar other than that
x <- rnorm(1000, mean = 8)
lines(density(x), col="blue")

#variation increase by increasing the sd value...so a bigger spread
lines(density(rnorm(1000,sd=2)), col="red")
lines(density(rnorm(1000,mean=8,sd=2)), col="green")

#still more variation due to bigger sd value
lines(density(rnorm(1000,sd=4)), col="purple")
lines(density(rnorm(1000,mean=8,sd=4)), col="cyan")


#Cleanse the ozone variable in the airquality data from missing values
my.ozone<-airquality$Ozone[!is.na(airquality$Ozone) & airquality$Ozone>1]

#my.ozone should be normally distributed, the best guess of the mean and 
#standard deviation would be as follows
mean.1<-mean(my.ozone)
sd.1<-sd(my.ozone)


#Simulate a number of normally distributed numbers with mean mean.1 and standard 
#deviation sd.1, equal to the amount of data in my.ozone
length(my.ozone)

set.seed(55789)
simulated.1<-rnorm(115,mean=mean.1,sd=sd.1)

#Compare the simulated values with my.ozone through qqplot() 
#so what we trying to do here is that we have generated/simulated normal distribution 
#values and now are  trying to find plot the values on graph and to see if there
#is linear correlation between the values(straight line). Shouldn't we have sorted the 
#values in the 2 distributions before plotting??
#after sorting I got the same plot but i do not know why i got the same result without
#sorting..Ans: Q-Q plots take your sample data, sort it in ascending order, and then 
#plot them 
#qqplot(c(23,22,25),c(16,8,10))
#qqplot(sort(c(22,23,25)),sort(c(16,8,10)))
qqplot(simulated.1,my.ozone)
lines(0:200,0:200,type="l",lwd=3,col="red")

qqplot(sort(simulated.1),sort(my.ozone))
lines(0:200,0:200,type="l",lwd=3,col="red")

#use ggplot to display same. Remember ggplot would not sort the data being plotted. 
#So to get the same plot as qqplot, sort the data
#df <- data.frame(sim=c(23,22,25), nonsim=c(16,8,10))
#df <- data.frame(sim=sort(c(23,22,25)), nonsim=sort(c(16,8,10)))
#ggplot(data = df, aes(x= sim, y= nonsim)) + geom_point(aes(colour=nonsim)) + geom_text(aes(x=sim-0.05, y=nonsim-0.15, label=df$sim))
df <- data.frame(sim=sort(simulated.1), nonsim=sort(my.ozone))
ggplot(data = df, aes(x=sim, y=nonsim))+ geom_point(aes(colour=nonsim)) 
                + geom_text(aes(x=sim-0.25, y=nonsim-2.45, 
                label=sprintf("%0.2f", round(df$sim, digits = 2)), size=1))+
          geom_text(aes(x=sim-0.25, y=nonsim+2.45, 
                    label=sprintf("%0.2f", round(df$nonsim, digits = 2)), size=1))


#consider a log-transform of the data. It might be that data are at different scales. 
#So take log of both arrays and then plot to see correlation. If my.ozone data was normally
#distributed, it's plot against data from a nomally distribution would fit around a 
#straight line. Wrong!!!that is now what we are trying to do here as we take 
#exponential for plot!
#If the log-transformed data should be normally distributed, a best guess on 
#mean and standard deviation would be as follows:
mean.2<-mean(log(my.ozone))
sd.2<-sd(log(my.ozone))

#now simulate the new values:
set.seed(8942)
simulated.2<-rnorm(115,mean=mean.2,sd=sd.2)

#Compare the exponential to the simulated points with my.ozone in a qqplot
#why did we take the exponential here?? because we generated simulated data
#from which scaled down(using log) ozone data. So now to compare them, we have to take
#exp of simulated data to cancel out the log (i think).
qqplot(exp(simulated.2),my.ozone)
lines(0:200,0:200,type="l",lwd=3,col="red")


doone <- function(){
  x <- sum(sample(1:6,2,replace=TRUE))
  y<-sum(sample(1:6,x,replace=TRUE))
  y
}

#simulate 1000 times
set.seed(457778) 
y.values<-replicate(1000,doone())
hist(y.values) 
box()
summary(y.values)
boxplot(y.values)


a<-rnorm(3, mean=2, sd=1)
b1<- c(3.373546, 4.183643, 3.164371)
summary(b1)
sd(b1)
b2<- c(1.373546, 2.183643, 1.164371)
summary(b2)
sd(b2)
b3<- c(5.373546, 6.183643, 5.164371 )
summary(b3)
sd(b3)
b4<- c(-2.626454, -1.816357, -2.835629)
summary(b4)
sd(b4)
#qqplot(a,b)
#lines(0:50,0:50,type="l",lwd=3,col="red")




n<-10000 
doone <- function(){ 
  x<-rbinom(1,50,1/6) 
  p<-x/50 
  p 
} 
p.sim<-replicate(n,doone())
length(p.sim)
summary(p.sim)

n<-100000
doone <- function(){ 
  x<-rbinom(1,50,1/6) 
  p<-x/50 
  p 
} 
p.sim<-replicate(n,doone()) 
hist(p.sim,breaks=20)


a<-rnorm(1000, mean=5, sd=1) 
hist(a)
plot(density(a))

my.data <- read.csv("Lab10.csv")
summary(my.data)
nrow(my.data)
ncol(my.data)

data1<-my.data$systolic.bp[my.data$Genotype=="BA"] 
data2<-my.data$systolic.bp[my.data$Genotype=="BB"] 

testResult <- t.test(data1,data2)
testResult

plot(density(data1))
lines(density(data2))

n1 <- length(my.data[my.data$Genotype=="BA",c("Genotype")])

testResult <- t.test(runif(26,min=1, max=2),data2)
testResult

set.seed(1234)
my.new.data<-my.data
my.new.data$Genotype<-"BB"
index.temp<-sample(1:50,n1)  
my.new.data$Genotype[index.temp]<-"BA"
new.data1<-my.new.data$systolic.bp[my.new.data$Genotype=="BA"]
new.data2<-my.new.data$systolic.bp[my.new.data$Genotype=="BB"]
t.test(new.data1,new.data2)$statistic


my.new.data<-my.data
my.new.data$Genotype<-"BB"
doone<-function(){
  index.temp<-sample(1:50,n1)  
  my.new.data$Genotype[index.temp]<-"BA"
  new.data1<-my.new.data$systolic.bp[my.new.data$Genotype=="BA"]
  new.data2<-my.new.data$systolic.bp[my.new.data$Genotype=="BB"]
  return(t.test(new.data1,new.data2)$statistic)
}
set.seed(554)
my.t.values<-replicate(100000,doone())
mean((1*(my.t.values)^2>2.027021^2))

df<-replicate(n=10,rnorm(2))
mx_df <-apply(df,2,max,drop=FALSE)
mean(mx_df)


#heatmap of a matrix color coded according to a condition
library('plot.matrix') 
a<-data.frame(rnorm(10),rnorm(10),rnorm(10),rnorm(10))
a<-as.matrix(a)
edit(a)
plt<-ifelse(a>.5,1,0)
# par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(plt)

#can also use image to produce a similar heatmap to visualize the values in a matrix
mat1 <- apply(plt, 2, rev)#this is visualizing the condition result on matrix values
image(t(mat1))

mat1 <- apply(a, 2, rev)#this is visualizing the matrix values
image(t(a))
