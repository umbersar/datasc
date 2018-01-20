library(ggplot2)
df <- read.csv('C:/temp/LabResults.csv', header = TRUE,sep =',' )
#df <- read.csv('LabResults_incomplete.csv', header = TRUE,sep =',' )
#head(df)
df <- subset(df, select = c(labm_code,labr_value))

for(labmCode in unique(df$labm_code)){
  filtered_df <- subset(df,labm_code==labmCode)  
  filtered_df <- droplevels(filtered_df)
  
  # png(filename=paste("C:/Users/sin17h/Documents/NatSoilDataSc/test/",labmCode,".png"))
  # boxplot(labr_value~labm_code,data=filtered_df,	xlab="labm_code", ylab="labr_value")
  # dev.off()
  
  pp <- ggplot2::ggplot(filtered_df,aes(x=labm_code, y=labr_value)) + ggplot2::geom_boxplot(outlier.shape = NA) + geom_jitter()
  # pp
  #ppp <-pp + geom_jitter()
  # ppp
  ggsave(paste("C:/Users/sin17h/Documents/NatSoilDataSc/test/",labmCode,".png"))
  #break;
}
  

# filtered_df <- df[df$labm_code=="4A1",]
# filtered_df <- subset(df,labm_code=="4A1")
# levels(droplevels(df$labm_code))
# df1 <- droplevels(filtered_df)
# pp <- ggplot2::ggplot(df1,aes(x=labm_code, y=labr_value)) + ggplot2::geom_boxplot()
# pp
# 
# dt <- read.table('C:/temp/LabResults.csv', header = TRUE,sep =',' )
# library(data.table)
# dt <- as.data.table(dt)
# dt <- subset(dt, select = c(labm_code,labr_value))
# dt <- dt[, mean.score := mean(labr_value), by = labm_code]
# dt <- dt[, sd.score := sd(labr_value), by = labm_code]
# dt <- dt[, max.score := max(labr_value), by = labm_code]
# dt <- dt[, min.score := min(labr_value), by = labm_code]
# dt <- dt[, length.score := length(labr_value), by = labm_code]
# dt <- dt[, outlier := abs(labr_value-mean.score) > 3 * sd.score, by = labm_code]
# dt <- dt[, sum.outlier := sum(outlier), by = labm_code]
# dt_out <- subset(dt, outlier==TRUE)  
# write.csv(dt, "C:/temp/stats.csv")
# write.csv(dt_out, "C:/temp/possible_outlier.csv")
# dt_fil <- subset(dt,labm_code=="4A1" & outlier==TRUE)  

