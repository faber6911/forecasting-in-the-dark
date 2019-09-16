
# import data -------------------------------------------------------------

data <- read.csv("D:/GIT/enercibiddding/PyScripts/dumps/q_p_list.csv")
head(data)
View(data)

# functions ---------------------------------------------------------------

rrr <- function(Y, X, Z=matrix(,0,0), rank=1, const=TRUE)
{
  n <- nrow(Y)
  k <- ncol(Y)
  if (nrow(Z) == n)
  {
    YY <- if (const) lm.fit(cbind(1,Z),Y)$residuals else lm.fit(Z,Y)$residuals
    #return(YY)
    B <- cancor(X, YY, xcenter = F, ycenter = F)$xcoef[,1:rank,drop=F]
    #return(cancor(X, YY, xcenter = F, ycenter = F)$xcoef)
    W <- X %*% B
    regr <- if (const) lm.fit(cbind(1,W,Z),Y) else lm.fit(cbind(W,Z),Y)
  }
  else
  {
    B <- cancor(X, Y, xcenter = const, ycenter = const)$xcoef[,1:rank,drop=F]
    W <- X %*% B
    #     form <- if (const) Y~W else Y~0+W
    #     regr <- lm(form)
    regr <- if (const) lm.fit(cbind(1,W),Y) else lm.fit(W,Y)
  }
  append(regr,list(B=B))
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

myLag <- function(dat, lag) data.frame(unclass(dat[c(rep(NA, lag), 1:(nrow(dat)-lag)),]))


rrrcv <- function(Y, X, Z=matrix(,0,0), step=10, ranks=1:ncol(Y), const=TRUE)
{
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
      #print(c(i,k,j, excl))
      out <- rrr(Y[-excl,], X[-excl,], Z[-excl,,drop=F], rank=i, const=const)
      #print(c(dim(Y[-excl,]), dim(X[-excl,]), dim(Z[-excl,,drop=F])))
      #return(Y[excl,] - outer(rep(1,step), out$coef[1,]))
      res <- Y[excl,] - outer(rep(1,step), out$coef[1,]) -
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
  colnames(mse) <- paste("Subsample",1:k)
  mse
}

#rrrcv(as.matrix(Y), as.matrix(X), as.matrix(Z))

# pre-processing ----------------------------------------------------------

#creating sinusoids

freq  <- outer(1:nrow(data), 1:16) * 2 * pi / 365*24
sinu365  <- cbind(cos(freq), sin(freq))
class(sinu365)
# sinu365 <- data.frame(sinu365)
# names(sinu365) <- paste0("sinu365.",1:32)
# dim(sinu365)

freq  <- outer(1:nrow(data), 1:6) * 2 * pi / 24
sinu24  <- cbind(cos(freq), sin(freq))
class(sinu24)
# sinu24 <- data.frame(sinu24)
# names(sinu24) <- paste0("sinu24.",1:12)
# dim(sinu24)

sinusoid <- cbind(sinu365, sinu24)




Y <- data[-((nrow(data)-168):nrow(data)),]
Y <- as.matrix(Y)
class(Y)
# names(Y) <- paste0("Y.", 0:50)

X1 <- data[2:(nrow(Y)+1),]
X1 <- as.matrix(X1)
class(X1)
# names(X1) <- paste0("X1.", 0:50)


X24 <- data[24:(nrow(Y)+23),]
X24 <- as.matrix(X24)
class(X24)
# names(X24) <- paste0("X24.", 0:50)


X168 <- data[168:(nrow(Y)+167),]
X168 <- as.matrix(X168)
class(X168)
# names(X168) <- paste0("X168.", 0:50)

sinusoid <- cbind(sinu365, sinu24)
class(sinusoid)
# head(sinusoid)

Z <- sinusoid[-((nrow(sinusoid)-168):nrow(sinusoid)),]
class(Z)
# dim(Z)

X <- cbind(X1,X24,X168)
class(X)
#df <- cbind(Y, X1, X24, X168, Z)
#View(df)

rrr(Y, X, Z)

rrrcv(Y, X, Z)

lag(data, 1)

###############################################






#transpose data matrix
apply(Z, 2, function(x) any(is.na(x)))

data <- t(data)

#creating y X and Z
y <- data[,ncol(data),drop=FALSE]
X <- cbind(data[,ncol(data)-1,drop=FALSE], data[,ncol(data)-24,drop=FALSE], data[,ncol(data)-168,drop=FALSE])# data[nrow(data)-24,,drop=FALSE], data[nrow(data)-168,,drop=FALSE])

dim(X)
dim(y)

Z <- cbind(sinu365, sinu24)
Z <- Z[nrow(Z),, drop=FALSE]
dim(Z)

# 
# dim(rbind(rep(sinu365[nrow(sinu365),,drop=FALSE],51)))
# dim(sinu365[nrow(sinu365),,drop=FALSE])
# Z <- rep.row(sinu365[nrow(sinu365),,drop=FALSE], 51)

Z <- rep.row(Z, 51)
dim(Z)

# check

dim(y)
dim(X)
dim(Z)

# models ------------------------------------------------------------------

lm.fit(cbind(1,Z), y)


rrr(y, X, Z, rank=1)

rrrcv(y, X, Z, step = 5)



################# shitty
# 
# library(dplyr)
# data %>% mutate_all(lag) -> X1
# head(lagged)                 


X1 <- myLag(data,1)
X24 <- myLag(data,24)
X168 <- myLag(data,168)
df <- cbind(data, X1, X24, X168, sinusoid)
finished <- df[complete.cases(df),]

Y <- as.matrix(finished[,1:51])

plot(x = seq(2,50), y = X1[2:50,1],type="l",col="red")
lines(y = data[2:50, 1], x = seq(2,50),col="green")
lines(y = X24[24:50, 1], x = seq(24,50),col="blue")

X <- as.matrix(finished[,52:204])
Z <- as.matrix(finished[,205:ncol(finished)])

# rrrcv(Y,X,Z)
# rrr(Y,X,Z)$B
# head(Z)
# library(svMisc)
# 
# mse <- c()
# pb <- txtProgressBar(min = 0, max = 153, style = 3)
# for (i in 1:153){
  
#setTxtProgressBar(pb, i)
YY <- lm.fit(cbind(1,Z),Y)$residuals
B <- cancor(X, YY, xcenter = F, ycenter = F)$xcoef[,1:51]
W <- X%*%B
Ypred <- lm.fit(cbind(1,X),Y)$fitted.values
res <- lm.fit(cbind(1,W,Z),Y)$residuals
dim(res)
mse <- mean(res^2)
mse

plot(x = seq(1,51), y = Ypred[nrow(Ypred),], type="l", col="red", xlab = "real", ylab = "pred", main = "pred vs real last time series")
lines(x = seq(1,51), y = Y[nrow(Y),], col="blue", lty = 2)

  #if (i == 153) cat("Done \n")
#}
#close(pb)


length(mse)

plot(y = mse, x = seq(1:153), type = "l")
