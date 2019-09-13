library(dplyr)

grid <- c(0.0000,2.5164,5.0128,7.5092,10.0056,12.5020,14.9984,17.4948,19.9912, 
          22.4876,24.984,27.4804,29.976,32.4732,34.969,37.4660,39.962,42.4588, 
          44.9552,47.4516,49.9480,52.444,54.9408,57.437,59.9336,62.430,64.9264, 
          67.4228,69.919,72.4156,74.912,77.4084,79.904,82.4012,84.897,87.3940, 
          89.8904,92.3868,94.8832,97.4396,100.0060,102.8324,105.6688,108.6952,
          111.8316,115.2660,119.3044,124.6716,129.5072,138.3468,259.0300)


sgrid <- function(dt, grid) {
  PQ <- dt %>%
    group_by(ENERGY_PRICE_NO) %>%
    summarise(quantity = sum(QUANTITY_NO)) %>%
    mutate(quantity = cumsum(quantity)) %>%
    rename(price = ENERGY_PRICE_NO)
  out <- sapply(grid, function(x) {y <- PQ$quantity[PQ$price <= x]; y[length(y)]})
  names(out) <- grid
  out
}

df <- read.csv("try.csv")
head(df)
names(df)
df$X <- NULL
head(df)

new_vec <- sgrid(df, grid)
names(new_vec)
new_df <- data.frame(new_vec)
new_df
names(new_df)

rownames(new_df)
new_df <- cbind(as.numeric(rownames(new_df)), data.frame(new_df, row.names=NULL))
new_df
colnames(new_df) <- c("price", "qnt")
head(new_df)
new_df[order(new_df$price),]
typeof(new_df$price)
arrange(new_df, price)
new_df
library(ggplot2)

ggplot(data = new_df, aes(x = price, y = qnt, group  = 1)) + geom_line()
