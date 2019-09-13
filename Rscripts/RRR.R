
# packages ----------------------------------------------------------------

library(quantmod)
library(dplyr)
library(Hmisc)
library(lubridate)
library(plotrix)
library(caroline)
library(pander)
library(plotly)

# import data -------------------------------------------------------------

rm(list = ls())
data <- read.csv("D:/GIT/enercibiddding/PyScripts/dumps/q_p_list_2015_2016.csv")
#data <- read.csv("D:/GIT/enercibiddding/PyScripts/dumps/time_s_2015_2016.csv")
dates <- read.csv("D:/GIT/enercibiddding/PyScripts/dumps/2015_2016_dates_unique_sorted.csv", header = FALSE)


#data_max <- max(data)
#data_min <- min(data)

# functions ---------------------------------------------------------------

#norm <- function(x){(x-min(x))/(max(x)-min(x))}

#back_norm <- function(x, max, min){x*(max-min)+min}

# data <- as.matrix(data)
# 
# p <- plot_ly(z = ~data[1:5000,]) %>% add_surface(
#   contours = list(
#     z = list(
#       show=TRUE,
#       usecolormap=TRUE,
#       highlightcolor="#ff0000",
#       project=list(z=TRUE)
#     )
#   )
# ) %>%
#   layout(
#     scene = list(
#       camera=list(
#         eye = list(x=1.87, y=0.88, z=-0.64)
#       )
#     )
#   )
# 
# p
# rm(p)
# data <- range01(data)
# data <- data.frame(data)

#plot(x = seq(1,51), y = data[1,], type = "l")

correct_zeros <- function(prediction_matrix){
  for (i in 1:nrow(prediction_matrix)){
    for (j in 1:ncol(prediction_matrix)){
      if (prediction_matrix[i,j] < 0){
        prediction_matrix[i,j] <- 0
      }
    }
  }
  return(prediction_matrix)
}

back_transformation <- function(transposed_matrix, real = FALSE){
  if (real == FALSE){
    for (i in 1:ncol(transposed_matrix)){
      if(i == 1){
        transposed_matrix[,i] <- exp(transposed_matrix[,i] + (std.error(transposed_matrix[,i])^2)/2)
      }else{
        transposed_matrix[,i] <- exp(transposed_matrix[,i] + (std.error(transposed_matrix[,i])^2)/2) + transposed_matrix[,i-1] - 1
      }
    }
  }else if (real == TRUE){
    for (i in 1:ncol(transposed_matrix)){
      if(i == 1){
        transposed_matrix[,i] <- exp(transposed_matrix[,i])
      }else{
        transposed_matrix[,i] <- exp(transposed_matrix[,i]) + transposed_matrix[,i-1] - 1
      }
    }
  }
  return(transposed_matrix)
}

dataframe.Lag <- function(dat, lag) data.frame(unclass(dat[c(rep(NA, lag), 1:(nrow(dat)-lag)),]))

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

rrrcv <- function(Y, X, Z=matrix(,0,0), step=10, ranks=1:ncol(Y), const=TRUE){
  
  n <- nrow(Y)
  k <- n %/% step
  mse <- matrix(0,length(ranks),step)
  cat("\nCross validation for RRR\n")
  h = 0
  best = c(Inf,NA)
  for (i in ranks)
  {
    h = h+1
    cat("Estimating rank =",i)
    for (j in 1:step)
    {
      excl <- (k*(j-1)+1):(k*j)
      out <- rrr(Y[-excl,], X[-excl,], Z[-excl,,drop=F], rank=i, const=const)
      res <- Y[excl,] - outer(rep(1,length(excl)), out$coef[1,]) -
        X[excl,]%*%out$B%*%out$coef[2:(i+1),] -
        Z[excl,,drop=F]%*%out$coef[-(1:(i+1)),]
      
      mse[h,j] <- mean(res^2)
    }
    rmse = sqrt(mean(mse[h,]))
    if (rmse < best[1])
    {
      best[1] = rmse
      best[2] = i
    }
    cat("  RMSE =", rmse,"\n")
  }
  cat("\nBest rank =",best[2],"( RMSE =", best[1],")\n")
  rownames(mse) <- ranks
  return(list(mse = mse,
              best_rank = best[1],
              best_rank_index = best[2]))
}


rrr <- function(Y, X, Z=matrix(,0,0), rank=1, const=TRUE)
{
  n <- nrow(Y)
  k <- ncol(Y)
  if (nrow(Z) == n)
  {
    YY <- if (const) lm.fit(cbind(1,Z),Y)$residuals else lm.fit(Z,Y)$residuals
    B <- cancor(X, YY, xcenter = F, ycenter = F)$xcoef[,1:rank,drop=F]
    W <- X %*% B
    regr <- if (const) lm.fit(cbind(1,W,Z),Y) else lm.fit(cbind(W,Z),Y)
  }
  else
  {
    B <- cancor(X, Y, xcenter = const, ycenter = const)$xcoef[,1:rank,drop=F]
    W <- X %*% B
    regr <- if (const) lm.fit(cbind(1,W),Y) else lm.fit(W,Y)
  }
  append(regr,list(B=B))
}

make_Z <- function(options, df){
  if (options == "one"){
    Z <- as.matrix(select(df, contains("sinu365")))
    rownames(Z) <- NULL
  }else if (options == "two"){
    Z <- as.matrix(select(df, contains("sinu365")))
    Z <- cbind(Z, as.matrix(select(df, contains("dum"))))
    rownames(Z) <- NULL
  }else if (options == "three"){
    Z <- as.matrix(select(df, contains("sinu")))
    Z <- cbind(Z, as.matrix(select(df, contains("dum"))))
    rownames(Z) <- NULL
  }
  return(Z)
}

make_model <- function(prediction, Z_parameter, test_size, cross_validation = TRUE, matrix_rank = 0){

  if (prediction == "1-hour"){
    dframe <- cbind(data, Lag1, Lag24, Lag168, sinu365, sinu24, dates[,-c(1,2)])
    dframe <- dframe[complete.cases(dframe),]
  }else if (prediction == "24-hours"){
    dframe <- cbind(data, Lag24, Lag168, sinu365, sinu24, dates[,-c(1,2)])
    dframe <- dframe[complete.cases(dframe),]
  }else if (prediction == "168-hours"){
    dframe <- cbind(data, Lag168, sinu365, sinu24, dates[,-c(1,2)])
    dframe <- dframe[complete.cases(dframe),]
  }
  dates <- dates[-(1:168),,drop=FALSE]
  Y <- as.matrix(select(dframe, contains("Y.")))
  X <- as.matrix(select(dframe, contains("Lag")))
  rownames(X) <- NULL
  Z <- make_Z(Z_parameter, dframe)
  
  if (cross_validation == TRUE){
    rank <- rrrcv(Y, X, Z)$best_rank_index
  }else{
    rank <- matrix_rank
  }
  
  test <- (nrow(Y)-(test_size-1)):nrow(Y)
  
  out <- rrr(Y[-test,], X[-test,], Z[-test,,drop=F], rank=rank)
  
  #Y_t <- t(Y[-test,])
  #Y_t <- back_norm(Y[-test,], data_max, data_min)
  Y_t <- back_transformation(Y[-test,]
                             #,data_max
                             #,data_min
                             ,real = TRUE)
  
  Yp_t <- out$fitted.values
  #Yp_t <- back_norm(Yp_t, data_max, data_min)
  Yp_t <- back_transformation(Yp_t
                              #,data_max
                              #,data_min
                              #,real = TRUE
                              )
  rmse_insample <- sqrt(mean((Y_t - Yp_t)^2))
  mape_insample <- mean(abs((Y_t-Yp_t)/Y_t))*100
  #rmse_log_insample <- sqrt(mean((Y[-test,]-out$fitted.values)^2))
  rmse_log_insample <- sqrt(mean(out$residuals^2))
  cat("\nRMSE_log_insample: ", rmse_log_insample)
  cat("\nRMSE_insample: ", rmse_insample)
  cat("\nMAPE_insample: ", mape_insample)
  
  res <- Y[test,] - outer(rep(1,length(test)), out$coef[1,]) -
    X[test,]%*%out$B%*%out$coef[2:(rank+1),] -
    Z[test,,drop=F]%*%out$coef[-(1:(rank+1)),]
  
  Ypred <- outer(rep(1,length(test)), out$coef[1,]) +
    X[test,]%*%out$B%*%out$coef[2:(rank+1),] +
    Z[test,,drop=F]%*%out$coef[-(1:(rank+1)),]
  
  Ypred <- correct_zeros(Ypred)
  Yreal <- Y[test,]
  mse_log <- mean(res^2)
  rmse_log <- sqrt(mse_log)
  mape_log <- mean(abs(res/Yreal))*100
  cat("\nModel RMSE_log: ", rmse_log)
  cat("\nModel MAPE_log: ", mape_log)
  #cat("\nPrint example plot")
  
  par(mfrow = c(2,1))
  plot(y = Y[nrow(Y),],
       x = seq(1,ncol(Y)),
       type = "l", col = "red",
       lwd = 1.5
       #,ylim = c(min(Y[nrow(Y),]), (max(Y[nrow(Y),])-1000))
       ,main = "Example plot with last date")
  lines(y = Ypred[nrow(Ypred),],
        x = seq(1,ncol(Ypred)),
        col = "blue",
        lwd = 1.5, lty = 2)
  
  #Ypred_t <- t(Ypred)
  rownames(Ypred) <- dates[test,1]
  #Ypred <- back_norm(Ypred, data_max, data_min)
  Ypred_t <- back_transformation(Ypred
                                 #,data_max
                                 #,data_min
                                 ,real = TRUE)
  #Yreal_t <- t(Yreal)
  rownames(Yreal) <- dates[test,1]
  #Yreal <- back_norm(Yreal, data_max, data_min)
  Yreal_t <- back_transformation(Yreal
                                 #,data_max
                                 #,data_min
                                 #,real = TRUE
                                 )
  res_t <- Yreal_t - Ypred_t
  rmse <- sqrt(mean(res_t^2))
  mape <- mean(abs(res_t/Yreal_t))*100
  cat("\nRMSE: ", rmse)
  cat("\nMAPE: ", mape)
  res_t <- data.frame(abs(res_t/Yreal_t))
  res_w <- res_t
  rownames(res_w) <- NULL
  res_w <- cbind(res_w, dates[test,2,drop=FALSE])
  res_w <- res_w %>% 
    group_by(weekdays) %>% summarise_all("mean")
  columns <- ncol(res_w)
  res_w <- data.frame(res_w)
  for (i in 1:nrow(res_w)){
    res_w[i,"mape"] <- sum(res_w[i,2:columns])/(columns-1)*100
  }
  mape_week <- res_w[,c(1,columns+1)]
  cat("\nWeekdays MAPE:\n")
  pander(mape_week)
  plot(x = seq(1,ncol(Yreal_t)),
       y = Yreal_t[nrow(Yreal_t),],
       type = "l",
       col = "red",
       ylim = c((min(Yreal_t[nrow(Yreal_t),])-1000), (max(Yreal_t[nrow(Yreal_t),])+1000)),
       lwd = 1.5)
  lines(x = seq(1,ncol(Ypred_t)), y =  Ypred_t[nrow(Ypred_t),], col = "blue", lwd = 1.5, lty = 2)
    
  return(list(df = df,
              Y = Y,
              X = X,
              Z = Z,
              rank = rank,
              rmse_log_insample = rmse_log_insample,
              rmse_insample = rmse_insample,
              mape_insample = mape_insample,
              rmse_log = rmse_log,
              mape_log = mape_log,
              RMSE = rmse,
              MAPE = mape,
              WeekD_MAPE = mape_week))
}

# pre processing ----------------------------------------------------------


dates$weekdays <- wday(as.Date(dates$V1), label = TRUE)

dates <- dummy_creation(dates, "weekdays", "lun")
dates <- dummy_creation(dates, "weekdays", "sab")
dates <- dummy_creation(dates, "weekdays", "dom")
dates[50,]
#dim(dates[-(1:168),,drop=FALSE])
#data <- norm(data)
rownames(data) <- dates[,1]
names(data) <- paste0("Y.", seq(1,ncol(data)))
Lag1 <- dataframe.Lag(data, 1)
names(Lag1) <- paste0("Lag1.", seq(1,ncol(Lag1)))
Lag24 <- dataframe.Lag(data, 24)
names(Lag24) <- paste0("Lag24.", seq(1,ncol(Lag24)))
Lag168 <- dataframe.Lag(data, 168)
names(Lag168) <- paste0("Lag168.", seq(1,ncol(Lag168)))
freq  <- outer(1:nrow(data), 1:16) * 2 * pi / 365.25*24
sinu365  <- cbind(cos(freq), sin(freq))
colnames(sinu365) <- paste0("sinu365.", seq(1,ncol(sinu365)))
freq  <- outer(1:nrow(data), 1:6) * 2 * pi / 24
sinu24  <- cbind(cos(freq), sin(freq))
colnames(sinu24) <- paste0("sinu24.", seq(1,ncol(sinu24)))



# 1-hour predictions ------------------------------------------------------

#51
model_one_1 <- make_model(prediction = "1-hour",
                    Z_parameter = "one",
                    test_size = 1000,
                    cross_validation = FALSE
                    ,matrix_rank = 51
                    )

#51 
model_two_1 <- make_model(prediction = "1-hour",
                          Z_parameter = "two",
                          test_size = 1000,
                          cross_validation = FALSE
                          ,matrix_rank = 51
                          )

#51
model_three_1 <- make_model(prediction = "1-hour",
                          Z_parameter = "three",
                          test_size = 1000,
                          cross_validation = FALSE
                          ,matrix_rank = 51
                          )


# 24-hours predictions ----------------------------------------------------

#50
model_one_24 <- make_model(prediction = "24-hours",
                            Z_parameter = "one",
                            test_size = 1000,
                            cross_validation = FALSE
                            ,matrix_rank = 50
                            )

#50
model_two_24 <- make_model(prediction = "24-hours",
                            Z_parameter = "two",
                            test_size = 1000,
                            cross_validation = FALSE
                            ,matrix_rank = 50
                            )

#50
model_three_24 <- make_model(prediction = "24-hours",
                            Z_parameter = "three",
                            test_size = 1000,
                            cross_validation = FALSE
                            ,matrix_rank = 50
                            )


# 168-hours predictions ---------------------------------------------------


#51 
model_one_168 <- make_model(prediction = "168-hours",
                           Z_parameter = "one",
                           test_size = 1000,
                           cross_validation = FALSE
                           ,matrix_rank = 51
                            )

#51 
model_two_168 <- make_model(prediction = "168-hours",
                           Z_parameter = "two",
                           test_size = 1000,
                           cross_validation = FALSE
                           ,matrix_rank = 51
                            )

#50
model_three_168 <- make_model(prediction = "168-hours",
                             Z_parameter = "three",
                             test_size = 1000,
                             cross_validation = FALSE
                             ,matrix_rank = 50
                              )

# make tables -------------------------------------------------------------

model_one_1$RMSE
model_two_1$RMSE
model_three_1$RMSE
model_one_24$RMSE
model_two_24$RMSE
model_three_24$RMSE
model_one_168$RMSE
model_two_168$RMSE
model_three_168$RMSE


model_one_1$MAPE
model_two_1$MAPE
model_three_1$MAPE
model_one_24$MAPE
model_two_24$MAPE
model_three_24$MAPE
model_one_168$MAPE
model_two_168$MAPE
model_three_168$MAPE


model_one_1$WeekD_MAPE
model_two_1$WeekD_MAPE
mean(model_three_1$WeekD_MAPE[,2])
mean(model_one_24$WeekD_MAPE[,2])
mean(model_two_24$WeekD_MAPE[,2])
mean(model_three_24$WeekD_MAPE[,2])
mean(model_one_168$WeekD_MAPE[,2])
model_two_168$WeekD_MAPE
model_three_168$WeekD_MAPE



# plot supply function ----------------------------------------------------

df <- cbind(data, Lag1, Lag24, Lag168, sinu365, sinu24, dates[,-c(1,2)])
df_24 <- cbind(data, Lag24, Lag168, sinu365, sinu24, dates[,-c(1,2)])
df_168 <- cbind(data, Lag168, sinu365, sinu24, dates[,-c(1,2)])
df <- df[complete.cases(df),]
df_24 <- df_24[complete.cases(df_24),]
df_168 <- df_168[complete.cases(df_168),]

Y <- as.matrix(select(df, contains("Y.")))
X <- as.matrix(select(df, contains("Lag")))
Z <- make_Z("three", df)
X_24<- as.matrix(select(df_24, contains("Lag")))
Z_24 <- make_Z("three", df_24)
X_168 <- as.matrix(select(df_168, contains("Lag")))
Z_168 <- make_Z("three", df_168)

test <- (nrow(Y)-999):nrow(Y)

out <- rrr(Y[-test,], X[-test,], Z[-test,,drop=F], rank=51)
out_24 <- rrr(Y[-test,], X_24[-test,], Z_24[-test,,drop=F], rank=51)
out_168 <- rrr(Y[-test,], X_168[-test,], Z_168[-test,,drop=F], rank=51)

res <- Y[test,] - outer(rep(1,length(test)), out$coef[1,]) -
  X[test,]%*%out$B%*%out$coef[2:(51+1),] -
  Z[test,,drop=F]%*%out$coef[-(1:(51+1)),]
res_24 <- Y[test,] - outer(rep(1,length(test)), out_24$coef[1,]) -
  X_24[test,]%*%out_24$B%*%out_24$coef[2:(51+1),] -
  Z_24[test,,drop=F]%*%out_24$coef[-(1:(51+1)),]
res_168 <- Y[test,] - outer(rep(1,length(test)), out_168$coef[1,]) -
  X_168[test,]%*%out_168$B%*%out_168$coef[2:(51+1),] -
  Z_168[test,,drop=F]%*%out_168$coef[-(1:(51+1)),]

Ypred <- outer(rep(1,length(test)), out$coef[1,]) +
  X[test,]%*%out$B%*%out$coef[2:(51+1),] +
  Z[test,,drop=F]%*%out$coef[-(1:(51+1)),]
Ypred_24 <- outer(rep(1,length(test)), out_24$coef[1,]) +
  X_24[test,]%*%out_24$B%*%out_24$coef[2:(51+1),] +
  Z_24[test,,drop=F]%*%out_24$coef[-(1:(51+1)),]
Ypred_168 <- outer(rep(1,length(test)), out_168$coef[1,]) +
  X_168[test,]%*%out_168$B%*%out_168$coef[2:(51+1),] +
  Z_168[test,,drop=F]%*%out_168$coef[-(1:(51+1)),]


Ypred <- back_transformation(Ypred)
Ypred_24 <- back_transformation(Ypred_24)
Ypred_168 <- back_transformation(Ypred_168)

Yreal <- Y[test,]
Yreal <- back_transformation(Yreal, real = TRUE)

sqrt(mean((Yreal-Ypred)^2))
mean(abs((Yreal-Ypred_168)/Yreal))*100

res <- Yreal-Ypred
res_24 <- Yreal -Ypred_24
res168 <- Yreal -Ypred_168
best <- c(999999999,1)
for (i in 1:nrow(res)){
  rmse <- sqrt(mean(res[i,]^2))+sqrt(mean(res_24[i,]^2))+sqrt(mean(res_168[i,]^2))
  if (rmse < best[1]){
    best[1] <- rmse
    best[2] <- i
  }
}
best

try <- 909
#par(mfrow = c(1,1))
plot(x = seq(1:ncol(Yreal)), y = Yreal[try,], col = "red", type = "l", lwd = 1.5, lty = 1, ylim=c(16500, 21000))
lines(x = seq(1:ncol(Ypred)), y = Ypred[try,], col = "blue", type = "l", lwd = 1.5, lty = 2)
lines(x = seq(1:ncol(Ypred_24)), y = Ypred_24[try,], col = "green", type = "l", lwd = 1.5, lty = 3)
lines(x = seq(1:ncol(Ypred_168)), y = Ypred_168[try,], col = "orange", type = "l", lwd = 1.5, lty = 4)

library(ggplot2)
library(ggthemes)

viz <- data.frame(x = seq(1:ncol(Yreal)), y = Yreal[806,])
class(viz)
tiff("test.tiff", units="in", width=5, height=5, res=300)
ggplot(aes(x = x, y = y), data = viz) + 
  geom_line(color = "red", size = 1)+ 
  xlab("Prices (50-iles)") +
  ylab("Quantity (MWh)") +
  theme(panel.background = element_rect("black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background=element_rect(fill = "black"),
        axis.title = element_text(colour = "red"),
        axis.line = element_line(colour = "red"),
        axis.text = element_text(colour = "red"))
dev.off()

# deprecated --------------------------------------------------------------

df <- cbind(data, Lag1, Lag24, Lag168, sinu365, sinu24, dates[,-c(1,2)])
df <- df[complete.cases(df),]

Y <- as.matrix(select(df, contains("Y.")))
X <- as.matrix(select(df, contains("Lag")))
Z <- make_Z("three", df)

dim(Y)
dim(X)
dim(Z)

# for 24-hours predictions

df_24 <- cbind(data, Lag24, Lag168, sinu365, sinu24, dates[,-c(1,2)])
df_24 <- df_24[complete.cases(df_24),]
#head(df_24, 2)

X_24 <- as.matrix(select(df_24, contains("Lag")))
rownames(X_24) <- NULL

# models ------------------------------------------------------------------

# cross validation in order to find optimal rank

rrrcv(Y, X, Z) #51
rrrcv(Y, X_24, Z) #50

# define test set

test <- (nrow(Y)-999):nrow(Y)


# model 1-hour forecasting

out <- rrr(Y[-test,], X[-test,], Z[-test,,drop=F], rank=51)
#out$residuals

#dim(out$fitted.values)
res <- back_transformation(t(out$residuals))
sqrt(mean(res^2))
Yt <- back_transformation(t(Y[-test,]))
mean(abs(res/Yt))

sqrt(mean((Y[-test,]-out$fitted.values)^2))
Y_t <- t(Y[-test,])
Y_t <- back_transformation(Y_t, data_max, data_min)
Yp_t <- t(out$fitted.values)
Yp_t <- back_transformation(Yp_t, data_max, data_min)
rmse_insample <- sqrt(mean((Y_t - Yp_t)^2))
rmse_insample
mape_insample <- mean(abs((Y_t-Yp_t)/Y_t))*100
mape_insample



res <- Y[test,] - outer(rep(1,length(test)), out$coef[1,]) -
  X[test,]%*%out$B%*%out$coef[2:(51+1),] -
  Z[test,,drop=F]%*%out$coef[-(1:(51+1)),]

Ypred <- outer(rep(1,length(test)), out$coef[1,]) +
  X[test,]%*%out$B%*%out$coef[2:(51+1),] +
  Z[test,,drop=F]%*%out$coef[-(1:(51+1)),]

Yreal <- Y[test,]

dim(Ypred)
Ypred[1,1:10]
Yreal[1,1:10]-Ypred[1,1:10]
Yreal[1,1:10]
# 
# for (i in 1:nrow(Ypred)){
#   for (j in 1:ncol(Ypred)){
#     if (Ypred[i,j] < 0){
#       Ypred[i,j] <- 0
#     }
#   }
# }
# 
# count = 0
# for (i in 1:nrow(Ypred)){
#   for (j in 1:ncol(Ypred)){
#     if (Ypred[i,j] < 0){
#       count <- count + 1
#     }
#   }
# }
# count


Ypred[1,1]

Yreal[1,1]
dim(res)

# for (i in 1:ncol(res)){
#   if(i == 1){
#     res[,i] <- (exp(res[,i] + (var(res[,i])/2)))
#   }else{
#     res[,i] <- (exp(res[,i] + (var(res[,i])/2)) + res[,i] - 1)
#   }
# }

mse <- mean(res^2)
rmse <- sqrt(mse)
mse

par(mfrow = c(1,1))
plot(y = Y[nrow(Y),], x = seq(1,ncol(Y)), type = "l", col = "red", lwd = 1.5)
lines(y = Ypred[nrow(Ypred),], x = seq(1,ncol(Ypred)), col = "blue", lwd = 1.5, lty = 2)

mean(abs(res/Yreal))*100

#plot(y = Y[,20], x = seq(1,nrow(Y)))
# #dim(Y)
# Ypred_t <- t(Ypred)
# colnames(Ypred_t) <- dates[test,1]
# head(Ypred_t[,2:6],2)
#Ypred_t <- back_transformation(Ypred_t, data_max, data_min)

for (i in 1:ncol(Ypred)){
  if(i == 1){
    Ypred[,i] <- exp(Ypred[,i] + (std.error(Ypred[,i])^2)/2)
  }else{
    Ypred[,i] <- exp(Ypred[,i] + (std.error(Ypred[,i])^2)/2) + Ypred[,i-1] - 1
  }
}


Ypred[1:2, 1:8]

dim(Ypred)


# Yreal_t <- t(Yreal)
# colnames(Yreal_t) <- dates[test,1]
# Yreal_t <- back_transformation(Yreal_t, data_max, data_min)

for (i in 1:ncol(Yreal)){
  if(i == 1){
    Yreal[,i] <- exp(Yreal[,i]) #+ (std.error(Yreal[,i])^2)/2)
  }else{
    Yreal[,i] <- exp(Yreal[,i] + (std.error(Yreal[,i])^2)/2) + Yreal[,i-1] - 1
  }
}

std.error(Yreal[,1])^2/2
dim(Yreal)


log(35)
exp(3.555)


Yr <- data_real[test,]
Yr[1:3,1:8]
Ypred[1:2,1:8]
rownames(Yreal) <- NULL
Yreal[1:3,1:8]
max(Yr-Yreal)
mean(resr^2)

#data_real <- data_real[(nrow(data_real)-1000):nrow(data_real),1:10]
#data_real_t <- t(data_real)
#Yreal_t[,1, drop=FALSE]
#data_real_t[,1,drop=FALSE]
Yr <- as.matrix(Yr)

res_t <- Yr - Yreal
sqrt(mean(res_t^2))
mean(abs(res_t/Yreal))*100

Yreal_t[35,803]
Ypred_t[35,803]

# model 24-hours forecasting

out <- rrr(Y[-test,], X_24[-test,], Z[-test,,drop=F], rank=50)
res <- Y[test,] - outer(rep(1,length(test)), out$coef[1,]) -
  X_24[test,]%*%out$B%*%out$coef[2:(50+1),] -
  Z[test,,drop=F]%*%out$coef[-(1:(50+1)),]

Ypred_24 <- Y[test,]-res
Ypred_24[1,1]
Yreal_24 <- Y[test,]
Yreal_24[1,1]

mse <- mean(res^2)
rmse <- sqrt(mse)
rmse

for (i in 1:nrow(Ypred_24)){
  for (j in 1:ncol(Ypred_24)){
    if (Ypred_24[i,j] < 0){
      Ypred_24[i,j] <- 0
    }
  }
}


count = 0
for (i in 1:nrow(Ypred_24)){
  for (j in 1:ncol(Ypred_24)){
    if (Ypred_24[i,j] < 0){
      count <- count + 1
    }
  }
}
count

plot(y = Yreal_24[nrow(Yreal_24),], x = seq(1,51), type = "l", col = "red", lwd = 1.5)
lines(y = Ypred_24[nrow(Ypred_24),], x = seq(1,51), col = "blue", lwd = 1.5, lty = 2)

Ypred_t_24 <- t(Ypred_24)
colnames(Ypred_t_24) <- dates[test,1]
head(Ypred_t_24[,2:6],2)


for (i in 1:nrow(Ypred_t_24)){
  for (j in 1:ncol(Ypred_t_24)){
      
    if(i == 1){
      Ypred_t_24[i,j] <- exp(Ypred_t_24[i,j]+((std.error(Ypred_t_24[i,]))^2)/2)
    }else{
      Ypred_t_24[i,j] <- exp(Ypred_t_24[i,j]+((std.error(Ypred_t_24[i,]))^2)/2) + Ypred_t_24[i-1,j] - 1
    }
  }
}

exp(Ypred_t_24[1,]+((std.error(Ypred_t_24[1,]))^2)/2)

Ypred_t_24[1,2:6]

Yreal_t_24 <- t(Yreal_24)
colnames(Yreal_t_24) <- dates[test,1]

for (i in 1:nrow(Yreal_t_24)){
  for (j in 1:ncol(Yreal_t_24)){
      
    if(i == 1){
      Yreal_t_24[i,j] <- (exp(Yreal_t_24[i,j] + ((std.error(Yreal_t_24[i,j])^2)/2)))
    }else{
      Yreal_t_24[i,j] <- (exp(Yreal_t_24[i,j] + ((std.error(Yreal_t_24[i,j])^2)/2)) + Yreal_t_24[i-1,j] - 1)
    }
  }
}

Yreal_t_24[1,]
res_t_24 <- Yreal_t_24 - Ypred_t_24
sqrt(mean(res_t_24^2))

# plot supply function

rm(list = setdiff(ls(), c("test", "Ypred_t", "Yreal_t", "Ypred_t_24", "Yreal_t_24", "pb")))

pb <- txtProgressBar(min = 1, max  = length(test), style = 3)
for (i in 1:length(test)){
  setTxtProgressBar(pb, i)
  plot(x = seq(1:nrow(Yreal_t)), y = Yreal_t[,i],
       type = "l", col = "red", lwd = 1.5, main = paste0("Time: ", i))
  lines(x = seq(1:nrow(Ypred_t)), y = Ypred_t[,i], col = "blue", lwd = 1.5, lty = 2)
  lines(x = seq(1:nrow(Ypred_t_24)), y = Ypred_t_24[,i], col = "darkgreen", lwd = 1.5, lty = 3)
  Sys.sleep(0.5)
}
close(pb)

i <- 803
plot(x = seq(1:nrow(Yreal_t)), y = Yreal_t[,i],
     type = "l", col = "red", lwd = 1.5, main = paste0("Time: ", i))
lines(x = seq(1:nrow(Ypred_t)), y = Ypred_t[,i], col = "blue", lwd = 1.5, lty = 2)
lines(x = seq(1:nrow(Ypred_t_24)), y = Ypred_t_24[,i], col = "darkgreen", lwd = 1.5, lty = 3)

res

# correction
res_t <- Yreal_t - Ypred_t
res_t <- data.frame(abs(res_t/Yreal_t))
res_w <- t(res_t)
colnames(res_w) <- NULL
res_w <- cbind(res_w, dates[test,2,drop=FALSE])
res_w <- res_w %>% 
  group_by(weekdays) %>% summarise_all("mean")
columns <- ncol(res_w)
res_w <- data.frame(res_w)
for (i in 1:nrow(res_w)){
  res_w[i,"mape"] <- sum(res_w[i,2:columns])/(columns-1)*100
}
mape_week <- res_w[,c(1,columns+1)]
mape_week

res <- res^2
res <- cbind(res, dates[test,2,drop=FALSE])
rownames(res) <- NULL
res <- res %>% 
  group_by(weekdays) %>% summarise_all("mean")
res <- data.frame(res)
res[,1]
for (i in 1:nrow(res)){
  res[i,"mape"] <- sum(res[i,2:52])/51
}
res[,c(1,ncol(res))]


