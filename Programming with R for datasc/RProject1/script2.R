c(1,2,3) #vector
1:5 #sequence which is internally a vector?
? typeof
? sample

? rnorm #normal distribution using random no.s rnorm(5,mean=2, sd=1)
? rgamma #gamma distribution using random no.s rgamma(5,shape=2, rate=2)
? rbinom #binomial distribution using random no.s rbinom(5,size=100, prob=.3)
? rmultinom #multi nominal distribution using random no.s rmultinom(5,size=100, prob=c(.3,.2,.5))
? runif #uniform distribution using random no.s runif(5,min=1, max=2)

? c #this is the concatenate function to produce a vector, not cat
? cat # catalogue function like print. Can also redirect output to a file
? paste # is for concatenating cahracter vectors. Like a '+' for string
? lapply #output is list. also tapply given below
? sapply #output is vector or matrix. sapply(df,class). also tapply given below
#data.frame(5:7, c(14, 12, 13))  note that what u specify in dataframe constructor is equal length columns. Matrix can also
#be used to instantiate data.frame if the constructor of matrix is more useful for a use case.
? seq
? rep
? rbind # this takes inputs which becomes rows
? cbind # this takes inputs which becomes columns
? matrix # this takes one long vector and then you supply shape to transform the vector into the matrix
? mode
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
? dput #dput and dget also do something similar as dump and source but with some differences. But this is ascii serialization.
? dget #dput and dget also do something similar as dump and source but with some differences. But this is ascii serialization.
? save # serialize in r-binary format
? load # de-serialize in r-binary format
? ls #list r objects from workspace
? rm # remove r objects from workspace
? library(RODBC) # sql data connections
? merge #merge dataframes..sort of like inner join
? names #names(df) gets col names
? class #class of df column, i.e., numeric or factor, etc.
? str #str(df) gets object structure 
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
? tapply #tapply(df$height, df$gender, mean). could also be done by aggregate(height~gender, df, mean) AND by(df$height, df$gender, mean)
#tapply(df$height, list(df$gender, df$tmt), mean). could also be done by aggregate(height~gender+tmt, df, mean) AND by(df$height, list(df$gender, df$tmt), mean).
#but the way the output is returned is different for tapply, aggregate and by. But for large datasets, tapply is the fastest and aggregate the slowest
#note the difference between these 3 functions with table. Also similar is xtabs 
? aggregate(height~gender, df, mean) #see tapply
? by(df$height, df$gender, mean)  #see tapply
? attach #attached datasets so that you do have to refer them by name. But it has it's issues around hiding names. So avoid
? detach
? with #temporary attach
? searchpaths
? table #table(airquality$Ozone>80, airquality$Month). compare it with xtabs, tapply, aggregate and by
? xtabs # compare it with table, tapply, aggregate and by
? ftable
? dplyr  
? par # par(mfrow=c(3,3))

x <- rnorm(10000)
head(x)
mean(x)
var(x)
#p<-ggplot2::ggplo
t(data=x) + ggplot2::geom_histogram()

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
histogram(x)

set.seed(1234)
x <- rnorm(100, mean = .5, sd = .3)
mean(x)
sd(x)
histogram(x,axes=FALSE)
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

plot(fofx, -7, 7)
matrix(c(5, 4, 3, 2, 1, 0) + 2, nrow = 2) < 5
sin
? cat
? density
? hist
? "if"
? plot
? summary

summary(a)
sd(a)
?? skewed

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
sum(A[row <= column]) # A, row and column are 3 different matrices. c how other matrices themselves are used to index into a separate matrix(that too in a condition!)


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
my.data.frame <- read.table("data.exercise5.2.dat", skip = 2, header = TRUE, sep = ";", dec = ",")

my.data.frame <- read.csv("Exercise 5.3.csv", na.strings = c("","NA"), skip = 2)[, -1]
is.na(my.data.frame[4, 8])
is.na(my.data.frame[4, 7])
is.na(my.data.frame[3, 2])

my.data.frame <- read.table("Exercise 5.4a.txt", header = TRUE, skip = 2)
dimnames(my.data.frame)[[2]] <- c("", "", "")

l <- levels(unlist(read.table("Exercise 5.4a.txt", header = F, skip = 1, nrows = 1)))
dimnames(my.data.frame)[[2]] <- levels(unlist(read.table("Exercise 5.4a.txt", header = F, skip = 1, nrows = 1)))[2:length(l)]

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
write.table(additional.data, "Exercise 6.1.txt", row.names = FALSE, col.names = FALSE, append = T)

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

library(RODBC)
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

data.frame.x<-data.frame(names=c("Gretha","Robert","John","Heather"),
                         age=c(30,18,25,70),
                         nickname=c("Quicksilver","The Man","Nifty","Starlight"))

data.frame.y<-data.frame("Person_name"=c("William","Nancy","Charlotte","Henry"),
                         age=c(15,75,32,51),
                         "pet_dog"=c("King","Whity","Captain Vom","Doggie"))

data.frame.z<-merge(data.frame.x,data.frame.y)
data.frame.z<-merge(data.frame.y,data.frame.x,
                    by.x=c("Person_name","age"),
                    by.y=c("names","age"),all=TRUE)

s<-subset(airquality,!is.na(Ozone) & !is.na(Solar.R))#just keep in mind complete rows are being removed here
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
cycle(presidents) 
tapply(presidents,cycle(presidents),mean,na.rm=T)
time(presidents)
frequency(presidents)
deltat(presidents)

searchpaths()


tapply(airquality$Solar.R , cut(airquality$Wind,10), mean,na.rm=TRUE)

a <- table(airquality$Ozone>80, airquality$Month)
a
a2 <- addmargins(a)
a2

a3 <-prop.table(a,2)

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

par(mfrow=c(3,3))
set.seed(779)

# plot 9 histograms of 12 simulated standard normally distributed numbers. Add the density of the standard normal to the plot with the curve() function.
#as we increase the sample size, the distribution is visibly more normal. Not to say that the low sample sizes did not contains numbers from a normal distribution.
for(i in 1:9){
  hist(rnorm(12), probability=TRUE,main=paste("Histogram",i))
  curve(dnorm,add=TRUE,col="red",lwd=3)
}

for(i in 1:9){
  hist(rnorm(25), probability=TRUE,main=paste("Histogram",i))
  curve(dnorm,add=TRUE,col="red",lwd=3)
}

for(i in 1:9){
  hist(rnorm(10000), probability=TRUE,main=paste("Histogram",i))
  curve(dnorm,add=TRUE,col="red",lwd=3)
}


