
library(quantmod)
library(dplyr)
library(Hmisc)
library(lubridate)
library(plotrix)
library(caroline)
library(pander)
library(plotly)

dummy_creation <- function(df, column, day){
  tmp <- c()
  for (elem in df[,column]){
    if (elem == day){
      tmp <- c(tmp, 1)
    }else{
      tmp <- c(tmp, 0)
    }
  }
  name <- paste0("dum.", day)
  df[, name] <- tmp
  return(df)
}

# LSTM MAPE weekdays values -----------------------------------------------

#data <- read.csv("D:/GIT/enercibiddding/PyScripts/dumps/time_s_2015_2016.csv")
dates <- read.csv("D:/GIT/enercibiddding/PyScripts/dumps/2015_2016_dates_unique_sorted.csv", header = FALSE)
fitted_values <- read.csv("C:/Users/fabri/Downloads/Telegram Desktop/daily_pred_rank_30_row_330.csv")
data <- read.csv("C:/Users/fabri/Downloads/Telegram Desktop/y_real_daily_pred_rank_30_row_330.csv")

plot(x = seq(1:51), y = fitted_values[nrow(fitted_values)-5,], type = "l")


dim(data)
data$X <- NULL
dim(fitted_values)
fitted_values$X <- NULL
# dim(fitted_values)
# data[1:2,1:8]
# fitted_values[1:2,1:8]
# data[1,1]

names(fitted_values)
dates$weekdays <- wday(as.Date(dates$V1), label = TRUE)
dates <- dummy_creation(dates, "weekdays", "lun")
dates <- dummy_creation(dates, "weekdays", "sab")
dates <- dummy_creation(dates, "weekdays", "dom")
dates[50,]

# dates <- dates[-(1:168),]
# dim(dates)
# data <- data[-(1:168),]
# data <- data[(nrow(data)-499):nrow(data),]
# dim(data)
# data <- data[1:332,]
# dim(data)

test <- data
# dim(test)
# dim(data)
#test[1:2,1:8]
test <- as.matrix(test)
fitted_values <- as.matrix(fitted_values)
# dim(fitted_values)
# data[1:2,1:8]
# fitted_values[1:2,1:8]
res <- test - fitted_values
# res_a <- abs(res/test)
# mean(res_a)*100


res_w <- data.frame(abs(res/test))
rownames(res_w) <- NULL
res_w <- cbind(res_w, dates[(nrow(dates)-329):nrow(dates),2,drop=FALSE])
res_w <- res_w %>% 
  group_by(weekdays) %>% summarise_all("mean")
columns <- ncol(res_w)
res_w <- data.frame(res_w)
for (i in 1:nrow(res_w)){
  res_w[i,"mape"] <- sum(res_w[i,2:columns])/(columns-1)*100
}

mape_week <- res_w[,c(1,columns+1)]
pander(mape_week)
pander(mean(mape_week[,2]))
