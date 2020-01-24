Sys.getlocale()
# Sys.setlocale("LC_ALL", "C")
getwd()
setwd("C:/Users/sin17h/Downloads")
ls()

#ggplots drag and drop plot builder
#esquisse::esquisser()
#tidylog
#heatmaply

w = read.csv('who.csv')
str(w)
mn = mean(w$Over60)
mn
min(w$Over60)
w$Country[w$Over60==min(w$Over60)]
w$Country[w$LiteracyRate==max(w$LiteracyRate,na.rm = TRUE )]
which.max(w$LiteracyRate)

tapply(w$ChildMortality, w$Region, mean)

boxplot(c(1,2,3,4,5))

getwd()

a = read.csv("USDA.csv")
str(a)
summary(a)
boxplot(a$Sodium)
a$Description[a$Sodium==max(a$Sodium,na.rm = TRUE)]
table

cr = read.csv("mvtWeek1.csv")
summary(cr)
str(cr)
mvt =cr
max(cr$ID)
min(cr$Beat)
cr$Date[0]
DateConvert = as.Date(strptime(mvt$Date, "%m/%d/%y %H:%M"))
summary(DateConvert)

mvt$Month = months(DateConvert)
mvt$Weekday = weekdays(DateConvert)
mvt$Date = DateConvert

tapply(mvt$Month,  mvt$Month, length)
table(mvt$Month)
table(mvt$Month)==min(table(mvt$Month))

tapply(mvt$Weekday,  mvt$Weekday, length)
table(mvt$Weekday)
table(mvt$Weekday)==max(table(mvt$Weekday))

tapply(mvt$Arrest,  mvt$Month, sum, na.rm = TRUE)
tapply(mvt$Arrest,  mvt$Month, sum, na.rm = FALSE)
tapply(mvt$Arrest,  mvt$Month, sum)
max(tapply(mvt$Arrest,  mvt$Month, sum, na.rm = TRUE))
table(mvt$Arrest,mvt$Month)

hist(mvt$Date, breaks=100)

boxplot(mvt$Date)
boxplot(mvt$Date ~ mvt$Arrest)

table(mvt$Arrest,mvt$Year)
#15536/176105

sort(table(mvt$LocationDescription))

Top5 = subset(mvt, LocationDescription=="STREET" | LocationDescription=="PARKING LOT/GARAGE(NON.RESID.)" | LocationDescription=="ALLEY" | LocationDescription=="GAS STATION" | LocationDescription=="DRIVEWAY - RESIDENTIAL")

TopLocations = c("STREET", "PARKING LOT/GARAGE(NON.RESID.)", "ALLEY", "GAS STATION", "DRIVEWAY - RESIDENTIAL")
Top5 = subset(mvt, LocationDescription %in% TopLocations)

str(Top5)

unique(factor(Top5$LocationDescription))
unique(Top5$LocationDescription)
Top5$LocationDescription = factor(Top5$LocationDescription)
str(Top5)

table(Top5$LocationDescription)
table(Top5$LocationDescription, Top5$Arrest)

table(Top5$LocationDescription, Top5$Weekday)

IBM = read.csv('IBMStock.csv')
GE = read.csv('GEStock.csv')
ProcterGamble = read.csv('ProcterGambleStock.csv')
CocaCola = read.csv('CocaColaStock.csv')
Boeing = read.csv('BoeingStock.csv')
str(Boeing)

IBM$Date = as.Date(IBM$Date, "%m/%d/%y")
GE$Date = as.Date(GE$Date, "%m/%d/%y")
CocaCola$Date = as.Date(CocaCola$Date, "%m/%d/%y")
ProcterGamble$Date = as.Date(ProcterGamble$Date, "%m/%d/%y")
Boeing$Date = as.Date(Boeing$Date, "%m/%d/%y")

str(IBM)
head(sort(IBM$Date),n=1)
summary(IBM)
sd(ProcterGamble$StockPrice)

plot(CocaCola$Date, CocaCola$StockPrice)
plot(CocaCola$Date, CocaCola$StockPrice, type="l", col='red')
CocaCola$Date[CocaCola$StockPrice==max(CocaCola$StockPrice)]
CocaCola$Date[CocaCola$StockPrice==min(CocaCola$StockPrice)]
lines(ProcterGamble$Date, ProcterGamble$StockPrice, col='blue')
abline(v=as.Date(c("2000-03-01")), lwd=2)

plot(CocaCola$Date[301:432], CocaCola$StockPrice[301:432], type="l", col="red", ylim=c(0,210))
lines(ProcterGamble$Date[301:432], ProcterGamble$StockPrice[301:432], col="blue")
lines(IBM$Date[301:432], IBM$StockPrice[301:432], col="green")
lines(GE$Date[301:432], GE$StockPrice[301:432], col="purple")
lines(Boeing$Date[301:432], Boeing$StockPrice[301:432], col="orange")

abline(v=as.Date(c("1997-09-01")), lwd=2)
abline(v=as.Date(c("1997-11-01")), lwd=2)

abline(v=as.Date(c("2004-01-01")), lwd=2)
abline(v=as.Date(c("2005-12-01")), lwd=2)

summary(IBM$StockPrice)
tapply(IBM$StockPrice, months(IBM$Date), mean)>mean(IBM$StockPrice)

max(tapply(GE$StockPrice, months(GE$Date), mean))
max(tapply(CocaCola$StockPrice, months(CocaCola$Date), mean))

CPS = read.csv('CPSData.csv')
str(CPS)
summary(CPS)
sort(table(CPS$Industry))
sort(table(CPS$Citizenship)) 

table(CPS$Race, CPS$Hispanic)
h = subset(CPS, Hispanic=="1")
sort(tapply(h$Age,h$Race,length))>250

summary(CPS)

is.na(CPS$Married) 
table(CPS$Region, is.na(CPS$Married))
table(CPS$Sex, is.na(CPS$Married))
table(CPS$Age, is.na(CPS$Married))

is.na(CPS$MetroAreaCode) 
table(CPS$State,is.na(CPS$MetroAreaCode))
table(CPS$State,is.na(CPS$MetroAreaCode))==0


table(CPS$Region,is.na(CPS$MetroAreaCode))

sort(tapply(is.na(CPS$MetroAreaCode), CPS$State, mean))

MTC = read.csv('MetroAreaCodes.csv')
CC = read.csv('CountryCodes.csv')
str(MTC); nrow(MTC)
str(CC); str(CC)
summary(MTC)

CPS = merge(CPS, MTC, by.x="MetroAreaCode", by.y="Code", all.x=TRUE)
summary(CPS); str(CPS)
sort(table(CPS$MetroArea))

sort(tapply(!is.na(CPS$Hispanic),CPS$MetroArea,length))
sort(tapply(CPS$Hispanic,CPS$MetroArea,mean))

table(CPS$MetroArea, CPS$Race == "Asian")
sort(tapply(CPS$Race == "Asian", CPS$MetroArea, mean))>.20

sort(tapply(CPS$Education == "No high school diploma", CPS$MetroArea, mean,na.rm = TRUE))


CPS = merge(CPS, CC, by.x="CountryOfBirthCode", by.y="Code", all.x=TRUE)
summary(CPS$Country); str(CPS)
sort(tapply(CPS$Country != 'United States', CPS$Country, length))

h = subset(CPS, CPS$Country != 'United States' & CPS$Country != 'Canada' & CPS$Country != 'Mexico')
summary(h)

table(CPS$MetroArea == "New York-Northern New Jersey-Long Island, NY-NJ-PA", CPS$Country != "United States")
table(CPS$MetroArea, CPS$Country == "India")
sort(tapply(CPS$Country == "India", CPS$MetroArea, sum, na.rm=TRUE))

poll = read.csv('AnonymityPoll.csv')
summary(poll)
str(poll)
table(poll$Smartphone)
table(poll$Sex, poll$Region)

table(poll$Region, poll$State)

summary(poll$Internet.Use)
summary(poll$Smartphone)
table(poll$Smartphone, poll$Internet.Use) 

limited = subset(poll, poll$Internet.Use == 1 | poll$Smartphone == 1)
str(limited)
summary(limited)

summary(poll$Info.On.Internet)

table(poll$Info.On.Internet)
summary(poll$Info.On.Internet == 0)
summary(poll$Info.On.Internet == 11)

summary(poll$Worry.About.Info)
summary(!is.na(poll$Worry.About.Info))#790
summary(poll$Worry.About.Info==1)#386
table(limited$Worry.About.Info)#note that table does not give u data about na values

table(poll$Anonymity.Possible)

hist(limited$Age)
plot(limited$Age, limited$Info.On.Internet)
table(limited$Age, limited$Info.On.Internet)
table(limited$Age, limited$Info.On.Internet) == max(table(limited$Age, limited$Info.On.Internet))

jitter(c(1, 2, 3))

plot(limited$Age, limited$Info.On.Internet)
plot(jitter(limited$Age), jitter(limited$Info.On.Internet))

tapply(poll$Info.On.Internet, poll$Smartphone, summary)
tapply(poll$Tried.Masking.Identity, poll$Smartphone, table)


#linear regression
wine =read.csv("wine.csv")
str(wine);summary(wine)

model1 =lm(Price~AGST, data = wine)
str(model1); 
summary(model1)

model2 = lm(Price ~ AGST + HarvestRain, data=wine)
summary(model2)

ml1 = lm(Price ~ HarvestRain + WinterRain, data=wine)
summary(ml1)

cor(wine)
cor(wine$Price, wine$Year)
table(wine$Price, wine$Year)
plot(wine$Price, wine$Year)
plot(wine$Price, wine$Year, type="l", col='red')

cor(wine$HarvestRain, wine$WinterRain)

bball = read.csv("baseball.csv")
mball = subset(bball, Year <2002)
mball$RD = mball$RS - mball$RA
winsReg = lm(W~RD, data=mball)
wins = 80.8814 + 99*.1058

# VIDEO 1

# Read in the data
NBA = read.csv("NBA_train.csv")
str(NBA)


# VIDEO 2

# How many wins to make the playoffs?
table(NBA$W, NBA$Playoffs)

# Compute Points Difference
NBA$PTSdiff = NBA$PTS - NBA$oppPTS

# Check for linear relationship
plot(NBA$PTSdiff, NBA$W)

# Linear regression model for wins
WinsReg = lm(W ~ PTSdiff, data=NBA)
summary(WinsReg)


# VIDEO 3

# Linear regression model for points scored
PointsReg = lm(PTS ~ X2PA + X3PA + FTA + AST + ORB + DRB + TOV + STL + BLK, data=NBA)
summary(PointsReg)

# Sum of Squared Errors
PointsReg$residuals
SSE = sum(PointsReg$residuals^2)
SSE

# Root mean squared error
RMSE = sqrt(SSE/nrow(NBA))
RMSE

# Average number of points in a season
mean(NBA$PTS)

# Remove insignifcant variables
summary(PointsReg)

PointsReg2 = lm(PTS ~ X2PA + X3PA + FTA + AST + ORB + DRB + STL + BLK, data=NBA)
summary(PointsReg2)

PointsReg3 = lm(PTS ~ X2PA + X3PA + FTA + AST + ORB + STL + BLK, data=NBA)
summary(PointsReg3)

PointsReg4 = lm(PTS ~ X2PA + X3PA + FTA + AST + ORB + STL, data=NBA)
summary(PointsReg4)

# Compute SSE and RMSE for new model
SSE_4 = sum(PointsReg4$residuals^2)
RMSE_4 = sqrt(SSE_4/nrow(NBA))
SSE_4
RMSE_4




# VIDEO 4

# Read in test set
NBA_test = read.csv("NBA_test.csv")

# Make predictions on test set
PointsPredictions = predict(PointsReg4, newdata=NBA_test)

# Compute out-of-sample R^2
SSE = sum((PointsPredictions - NBA_test$PTS)^2)
SST = sum((mean(NBA$PTS) - NBA_test$PTS)^2)
R2 = 1 - SSE/SST
R2

# Compute the RMSE
RMSE = sqrt(SSE/nrow(NBA_test))
RMSE