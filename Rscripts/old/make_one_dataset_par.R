# ---- unzip files in the current folder

# zipfiles <- dir(pattern = "*.zip")
# files <- zipfiles
# substr(files, 29, 31) <- "xml"
#  
# for (i in 1:length(files)) unzip(zipfiles[i], files[i])

# ---- build data.frame from XML
library(xml2)

#msd <- NULL # will contain full dataframe
dataset_list <- NULL
fnames <- dir(pattern = "*.xml")
library(foreach)
#install.packages("doParallel")
library(doParallel)
#install.packages("doMC")
library(doMC)
registerDoMC(2)  #change the 2 to your number of CPU cores  

#this way, dataset_list is a list containing each dataset of the GME dir
#TODO: 
#iterate over the list, 
#convert each element to a dataset ()
#join the datasets
#save the final total dataset as a csv


dataset_list <- foreach (i=1:4) %dopar% {
  fname= fnames[i]
  print(i)
  cat("File", fname, date(), "\n")
  mb <- read_xml(fname)
  
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
  
  for (i in 1:length(var_keep)) {
    cat("Working on ", var_keep[i])
    var_content <- xml_find_all(mb, paste0("//", var_keep[i]))
    nrectmp <- length(var_content)
    cat("  n =", nrectmp, "\n")
    df[[i]] <- as(xml_text(var_content), typeof(df[[i]]))
  }
  #msd <- if (fname == fnames[1]) df else rbind(msd, df)
}

#something like
final_dataset <- data.frame()
for(i in 1:length(dataset_list)){
  dataset=dataset_list[i]
  dat <- data.frame(dataset)
  #join dat and final_dataset
  final_dataset <- if (i == 1) dat else rbind(final_dataset, dat)
}
final_dataset
#the whole parallelization is NOT TESTED THO!!

#save(msd, file = "MSDOffertePubbliche.Rdata")
#object.size(msd)
#length(fnames)
#5*133

write.csv(final_dataset, "final_dataset.csv")



