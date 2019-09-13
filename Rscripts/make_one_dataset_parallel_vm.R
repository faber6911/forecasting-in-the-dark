#!/usr/bin/env Rscript

# ---- unzip files in the current folder

# zipfiles <- dir(pattern = "*.zip")
# files <- zipfiles
# substr(files, 29, 31) <- "xml"
#  
# for (i in 1:length(files)) unzip(zipfiles[i], files[i])

start_t <- Sys.time()

# ---- build data.frame from XML
if(!require(xml2)) install.packages("xml2")
library(xml2)
setwd("/home/studente/dslabproject/Rscripts/")

#msd <- NULL # will contain full dataframe
dataset_list <- NULL
fnames <- dir(path = "../raw_data/", pattern = "*.xml")
#fnames
if(!require(foreach)) install.packages("foreach")
library(foreach)
if(!require(doParallel)) install.packages("doParallel")
library(doParallel)
#if(!require(doMC)) install.packages("doMC")
#library(doMC)
#registerDoMC(2)  #change the 2 to your number of CPU cores  

#this way, dataset_list is a list containing each dataset of the GME dir
#TODO: 
#iterate over the list, 
#convert each element to a dataset ()
#join the datasets
#save the final total dataset as a csv

#library(parallel)
if(!require(doSNOW)) install.packages("doSNOW")
library(doSNOW)

numCores<-detectCores()
cat("using", numCores, "cores" )
#numCores=2
cl <- makeCluster(numCores)
registerDoSNOW(cl)

# progress bar ------------------------------------------------------------
library(progress)

iterations <- length(fnames)                               # used for the foreach loop  
#iterations=4
pb <- progress_bar$new(
  format = "iteration = :letter [:bar] :elapsed | eta: :eta",
  total = iterations,    # 100 
  width = 60)

#progress_letter <- rep(LETTERS[1:4], 1)  # token reported in progress bar
prog <- 1:iterations
# allowing progress bar to be used in foreach -----------------------------
progress <- function(n){
  #pb$tick(tokens = list(letter = progress_letter[n]))
  pb$tick(tokens = list(letter = prog[n]))
  
} 

opts <- list(progress = progress)


#dataset_list <- foreach (i=1:4) %dopar% {
dataset_list <- foreach (i=1:iterations, .combine = rbind, .options.snow = opts) %dopar% {
  library(xml2)
  setwd("/home/studente/dslabproject/Rscripts/")
  fname= fnames[i]
  file_path=paste0("../raw_data/", fname)
  print(i)
  #cat("File", fname, date(), "\n")
  #cat("blah-blah-blah\n", file=stdout())
  mb <- read_xml(file_path)
  #mb
  # extract variable names and types
  schema <- as_list(xml_child(mb, 1))[[1]][[1]][[1]][[1]][[1]][[1]]
  var_keep  <- c("PURPOSE_CD", "STATUS_CD", "UNIT_REFERENCE_NO", "INTERVAL_NO", "BID_OFFER_DATE_DT",
                 "QUANTITY_NO", "AWARDED_QUANTITY_NO", "ENERGY_PRICE_NO",
                 "MERIT_ORDER_NO", "PARTIAL_QTY_ACCEPTED_IN", "ADJ_QUANTITY_NO", "ZONE_CD",
                 "AWARDED_PRICE_NO", "OPERATORE")
  var_names <-sapply(schema, function(x) attr(x, "name"))
  var_types <- sapply(schema, function(x) attr(x, "type"))
  
  nrecords <- length(xml_children(mb)) - 1
  df <- data.frame(PURPOSE_CD = character(nrecords),
                   STATUS_CD = character(nrecords),
                   UNIT_REFERENCE_NO = character(nrecords),
                   INTERVAL_NO = integer(nrecords),
                   BID_OFFER_DATE_DT = integer(nrecords),
                   QUANTITY_NO = numeric(nrecords),
                   AWARDED_QUANTITY_NO = numeric(nrecords),
                   ENERGY_PRICE_NO = numeric(nrecords),
                   MERIT_ORDER_NO = integer(nrecords),
                   PARTIAL_QTY_ACCEPTED_IN = character(nrecords),
                   ADJ_QUANTITY_NO = numeric(nrecords),
                   ZONE_CD = character(nrecords),
                   AWARDED_PRICE_NO = numeric(nrecords),
                   OPERATORE = character(nrecords),
                   stringsAsFactors = FALSE
                   )
  #cat(length(df))  
  for (i in 1:length(var_keep)) {
    #cat("Working on ", var_keep[i])
    var_content <- xml_find_all(mb, paste0("//", var_keep[i]))
    nrectmp <- length(var_content)
    #cat("  n =", nrectmp, "\n")
    df[[i]] <- as(xml_text(var_content), typeof(df[[i]]))
  }
  #msd <- if (fname == fnames[1]) df else rbind(msd, df)
  
  df
}

write.csv(dataset_list, "../data/final_dataset_vm.csv")


#SHOULD BE NO LONGER NEEDED

#something like
# final_dataset <- data.frame()
# for(i in 1:length(dataset_list)){
#   dataset=dataset_list[i]
#   dat <- data.frame(dataset)
#   #join dat and final_dataset
#   final_dataset <- if (i == 1) dat else rbind(final_dataset, dat)
# }
# final_dataset


#NB: the whole parallelization is NOT TESTED THO!!

#save(msd, file = "MSDOffertePubbliche.Rdata")
#object.size(msd)
#length(fnames)
#5*133


end_t <- Sys.time()

time <- end_t -start_t

cat("exec time", time)



