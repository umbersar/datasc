library(ggplot2)
df <- read.csv('C:/temp/LabResults.csv', header = TRUE, sep = ',')
#head(df)
df <- subset(df, select = c(labm_code, labr_value))

for (labmCode in unique(df$labm_code)) {
    filtered_df <- subset(df, labm_code == labmCode)
    filtered_df <- droplevels(filtered_df)

    # png(filename=paste("C:/Users/sin17h/Documents/NatSoilDataSc/test/",labmCode,".png"))
    # boxplot(labr_value~labm_code,data=filtered_df,	xlab="labm_code", ylab="labr_value")
    # dev.off()

    pp <- ggplot2::ggplot(filtered_df, aes(x = labm_code, y = labr_value)) + ggplot2::geom_boxplot()
    pp
    ppp <-pp + geom_jitter()
    ppp
    #ggsave(paste("C:/Users/sin17h/Documents/NatSoilDataSc/test/", labmCode, ".png"))
    break;
}


# filtered_df <- df[df$labm_code=="18F1_CU",]
# filtered_df <- subset(df,labm_code=="18F1_CU")
# levels(droplevels(df$labm_code))
# df1 <- droplevels(filtered_df)
# pp <- ggplot2::ggplot(df1,aes(x=labm_code, y=labr_value)) + ggplot2::geom_boxplot()
# pp
