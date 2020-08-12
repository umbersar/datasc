library(DBI)
con <- dbConnect(odbc::odbc(), Driver = "SQL Server", Server = "localhost\\sql2019", 
                 Database = "NatSoil", Trusted_Connection = "True")

# dbReadTable()
# a <- dbExecute(con,"select * from Test")
Lab_results_pivoted <- dbGetQuery(con,"select * from Test")
# write.csv(Lab_results_pivoted,"c:/temp/Lab_results_pivoted.csv")
summary(Lab_results_pivoted)
str(Lab_results_pivoted)
dim(Lab_results_pivoted)

#find if there are any empty string in the df
blank_cols = matrix(nrow = ncol(Lab_results_pivoted), ncol = 5)
totalCols = ncol(Lab_results_pivoted)
totalRows = nrow(Lab_results_pivoted)

#col1         col2                  col3        col4            col5
#'labm_code'  'empty strings count' 'NA count'  'Not NA count'  '% of NAs'
for (i in 1:ncol(Lab_results_pivoted)) {
  blank_cols[i,1] = colnames(Lab_results_pivoted)[i]
  blank_cols[i,2] = sum(Lab_results_pivoted[,i][Lab_results_pivoted[,i] == ''],
                        na.rm = TRUE)  
  blank_cols[i,3] = sum(is.na(Lab_results_pivoted[,i]) == TRUE)  
  blank_cols[i,4] = sum(!is.na(Lab_results_pivoted[,i]) == TRUE)  
  blank_cols[i,5] = sum(is.na(Lab_results_pivoted[,i]) == TRUE)/totalRows
}

# pdf(file = "C:/Temp/mat.pdf",width = 4, height = 4)
png(file = "C:/Temp/heatmap.png")

library('plot.matrix') 
plt<-ifelse(!is.na(Lab_results_pivoted),1,0)
plt<-as.matrix(plt)
dim(plt)
image(x=1:nrow(Lab_results_pivoted), y=1:ncol(Lab_results_pivoted),plt)

plot(plt)

dev.off()

write.csv(plt, "C:/Temp/heatmap.csv")
heatmap(plt)

sum(Lab_results_pivoted[,11] >0, na.rm = TRUE)
#this should not return True if there is not a single empty string which is what 
#we want
any(blank_cols[,2]>0,na.rm = TRUE)

#set all empty strings to NA if there are any
# Lab_results_pivoted==''
# Lab_results_pivoted[Lab_results_pivoted==''] = NA

#sort the matrix by the number of 'Not NA' values
blank_cols[order(as.integer(blank_cols[,4])),]

#there is not single row in df which has values for all columns. 
complete <- Lab_results_pivoted[complete.cases(Lab_results_pivoted),]
dim(complete)

#create correlation matrix. We still haven't dealt with NA values as yet.
cor_mat <- cor(Lab_results_pivoted)

#cor function handling of NAs is explained in 'learningstatisticswithR' book
cor_mat <- cor(Lab_results_pivoted,use="complete.obs")#this will remove the complete 
#row from feature matrix even if one feature had a NA value
cor_mat <- cor(Lab_results_pivoted,use="pairwise.complete.obs")#this will drop the
#the feature observation if and only if it is being used for pairwise correlation 
any(cor_mat!=1&& !is.na(cor_mat))

#synthetic data to show the problem. Few Random NAs make it difficult to make a 
# correlation matrix
set.seed(123)
x<-matrix(rnorm(9),3,3)
x[sample(9,1)]<-NA#increase the number of random NAs to notice the diff
x
#fix(x)

cor_x <- cor(x)
any(cor_x!=1 & !is.na(cor_x))#&& did not work!
which( cor_x != 1 & !is.na(cor_x)) 

cor(x,use="everything")#i think that is the default

cor_x <- cor(x,use="complete.obs")
any(cor_x!=1 & !is.na(cor_x))#&& did not work!
which( cor_x != 1 & !is.na(cor_x)) 

cor_x <- cor(x,use="pairwise.complete.obs")
any(cor_x!=1 & !is.na(cor_x))#&& did not work!
which( cor_x != 1 & !is.na(cor_x)) 

summary(x)



Lab_methods_description <- dbGetQuery(con,"select distinct lm.labm_code,lm.LABM_NAME,lm.LABM_SHORT_NAME, lm.LABM_UNITS
from HORIZONS as hr inner join [NatSoil].[dbo].[LAB_RESULTS] as lr
on hr.agency_code = lr.agency_code and hr.proj_code = lr.proj_code and hr.s_id = lr.s_id
and hr.o_id = lr.o_id and hr.h_no = lr.h_no
inner join LAB_METHODS as lm on lm.LABM_CODE = lr.labm_code")
write.csv(Lab_methods_description,"c:/temp/Lab_methods_description.csv")
