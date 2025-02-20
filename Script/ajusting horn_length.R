# Benjamin Larue
# Version 1          date: 2020-10-05
#  script pour obtenir les masse ajuste en juin et en septembre



###############################################################################.
#####     loading required libraries   ######
library(lme4)
library(plyr)
source("function.R")

# adjust POP --------------------------------------------------------------
bigRes <- data.frame()             # make empty data frame to house result
nsim=1000                         # set number of iteration for bootstrap used to get SE of mixed model
DJJ.eqal=1                         # if  observation id <= to DJJ.equal and missing, should use actual method
for(year in 1971:2023){           #  for all these years
  rm(list = c("d1", "data", "dt1", "EwtfMM", "EwtlbMM", "Ewtlbpop", "EwtMM",
              "Ewtoldfpop", "EwtoldMpop",
              "EwtyongMpop", "EwtyrlMM",
              "Ewtyrlpop",   "jun15", "jun5", "l.id",  "sept15",
              "WtfeMales", "Wtlb", "WtMales", "Wtyrl"))    # remove dis-used object
  
  data <- read.csv("Data/data_clean.csv", sep=";")               # load daily data
  data$ID <- as.factor(data$ID)
  data$sex <- as.factor(data$sex)
  data <- data[data$JJ>100,]                            # remove id with strange ID (cadomin transplant measuements?)
  
  # year=1978
  data <- data[!is.na(data$pheno),]
  
  # getting Marco style julian. days
  d1 <- as.Date(paste(year,"/05/25",sep=""),format="%Y/%m/%d")       #
  # d1 <- as.Date(paste(year,"/01/01",sep=""),format="%Y/%m/%d")       # un-comment to use real JJ
  data$JJ <- data$JJ + as.numeric(as.Date(paste(year,"/01/01",sep=""),format="%Y/%m/%d")-d1)+1
  jun5 <- as.numeric(as.Date(paste(year,"/06/05",sep=""),format="%Y/%m/%d")-d1)+1
  jun5 <- as.numeric(as.Date(paste(year,"/06/05",sep=""),format="%Y/%m/%d")-d1)+1
  jun15 <- as.numeric(as.Date(paste(year,"/06/15",sep=""),format="%Y/%m/%d")-d1)+1
  sept15 <- as.numeric(as.Date(paste(year,"/09/15",sep=""),format="%Y/%m/%d")-d1)+1
  
  #
  dt1 <- data[data$yr==year,]            #  subset for year to analyse
  
  dt1$JJ.sqrt <- sqrt(dt1$JJ)           #  calculate sqrt for adult models
  dt1 <- dt1[!is.nan(dt1$JJ.sqrt),]     # remove NAN; necessary when using Marco's date starting on 25 of may
  
  # run males
  EwtMM <- try(Male.mixed(dt1,nsim))    #  try fonction to estimate male with mixed model
  if(!is.data.frame(EwtMM)){EwtMM <- as.data.frame( # if  result is not a data frame, make empty data.frame of same form
    matrix(NA,nrow = 0,ncol = 7, dimnames = list(NULL,c("ID" , "JJ", "n", "age"  ,  "adjwt", "adjwtSE" ,"method"))))
  print(paste("failure for mixed model for males of",year ))    # and print an error message
  }
  
  EwtoldMpop <- try(Male.pop.old(dt1,nsim,max.DeltaJ = 50))    #  repeat for all necessary fonctions
  if(!is.data.frame(EwtoldMpop)){EwtoldMpop <- as.data.frame(
    matrix(NA,nrow = 0,ncol = 7, dimnames = list(NULL,c("ID" , "JJ", "n", "age"  ,  "adjwt", "adjwtSE" ,"method"))))
  print(paste("failure for pop model for old males of",year ))
  }
  
  EwtyongMpop <- try(Male.pop.young(dt1,nsim,max.DeltaJ = 50))
  if(!is.data.frame(EwtyongMpop)){EwtyongMpop <- as.data.frame(
    matrix(NA,nrow = 0,ncol = 7, dimnames = list(NULL,c("ID" , "JJ", "n", "age"  ,  "adjwt", "adjwtSE" ,"method"))))
  print(paste("failure for pop model for young males of",year ))
  }
  
  # compile male Res
  l.id <- droplevels(unique( dt1[dt1$sex=="male"& !is.na(dt1$pheno),"ID"]))  # find all ID
  WtMales <- data.frame(ID=l.id)                                             # put in dataframe
  WtMales <- merge(WtMales,EwtMM,all.x = T)                               #  add estimated pheno from mixed model
  l.id <- droplevels(unique( WtMales[is.na(WtMales$adjwt),"ID"]))  #  find id which where not adjusted
  WtMales <- rbind(WtMales,EwtoldMpop[EwtoldMpop$ID %in%l.id,])     # fill with adjusted pheno from pop models
  WtMales <- rbind(WtMales,EwtyongMpop[EwtyongMpop$ID %in%l.id,])
  WtMales <- WtMales[!is.na(WtMales$method),]                       # remove ID which could not be adjusted by anything
  #   WtMales$sex="male"
  if(nrow(WtMales)>0) WtMales$sex="male"                            #  add sex because I forgot to do it in the function
  
  # run female
  EwtfMM <- try(feMale.mixed(dt1,nsim))
  if(!is.data.frame(EwtfMM)){EwtfMM <- as.data.frame(
    matrix(NA,nrow = 0,ncol = 7, dimnames = list(NULL,c("ID" , "JJ", "n", "age"  ,  "adjwt", "adjwtSE" ,"method"))))
  print(paste("failure for mixed model for females of",year ))
  }
  
  Ewtoldfpop <- try(feMale.pop.old(dt1,nsim,max.DeltaJ = 50))
  if(!is.data.frame(Ewtoldfpop)){Ewtoldfpop <- as.data.frame(
    matrix(NA,nrow = 0,ncol = 7, dimnames = list(NULL,c("ID" , "JJ", "n", "age"  ,  "adjwt", "adjwtSE" ,"method"))))
  print(paste("failure for pop model for old females of",year ))
  }
  
  
  # compile female Res
  l.id <- droplevels(unique( dt1[dt1$sex=="female"& !is.na(dt1$pheno),"ID"]))
  WtfeMales <- data.frame(ID=l.id)
  WtfeMales <- merge(WtfeMales,EwtfMM,all.x = T)
  l.id <- droplevels(unique( WtfeMales[is.na(WtfeMales$adjwt),"ID"]))
  WtfeMales <- rbind(WtfeMales,Ewtoldfpop[Ewtoldfpop$ID %in%l.id,])
  WtfeMales <- WtfeMales[!is.na(WtfeMales$method),]
  if(nrow(WtfeMales)>0) WtfeMales$sex="female"
  
  # run lb
  EwtlbMM <- try(lb.mixed(dt1,nsim))
  if(!is.data.frame(EwtlbMM)){EwtlbMM <- as.data.frame(
    matrix(NA,nrow = 0,ncol = 8, dimnames = list(NULL,c("ID" , "JJ", "n", "age"  ,  "adjwt", "adjwtSE" ,"method","sex"))))
  print(paste("failure for mixed model for lb of",year ))
  }
  
  Ewtlbpop <- try(lb.pop(dt1,nsim,max.DeltaJ = 50))
  if(!is.data.frame(Ewtlbpop)){Ewtlbpop <- as.data.frame(
    matrix(NA,nrow = 0,ncol = 8, dimnames = list(NULL,c("ID" , "JJ", "n", "age"  ,  "adjwt", "adjwtSE" ,"method", "sex"))))
  print(paste("failure for pop model for lb of",year ))
  }
  
  # compile lb Res
  l.id <- droplevels(unique( dt1[dt1$age==0 & !is.na(dt1$pheno),"ID"]))
  Wtlb <- data.frame(ID=l.id)
  Wtlb <- merge(Wtlb,EwtlbMM,all.x = T)
  l.id <- droplevels(unique( Wtlb[is.na(Wtlb$adjwt),"ID"]))
  Wtlb <- rbind(Wtlb,Ewtlbpop[Ewtlbpop$ID %in% l.id,])
  Wtlb <- Wtlb[!is.na(Wtlb$method),]
  
  ## run yrl
  EwtyrlMM <- try(yrl.mixed(dt1,nsim))
  if(!is.data.frame(EwtyrlMM)){EwtyrlMM <- as.data.frame(
    matrix(NA,nrow = 0,ncol = 8, dimnames = list(NULL,c("ID" , "JJ", "n", "age"  ,  "adjwt", "adjwtSE" ,"method","sex"))))
  print(paste("failure for mixed model for yrl of",year ))
  }
  
  Ewtyrlpop <- try(yrl.pop(dt1,nsim,max.DeltaJ = 50))
  if(!is.data.frame(Ewtyrlpop)){Ewtyrlpop <- as.data.frame(
    matrix(NA,nrow = 0,ncol = 8, dimnames = list(NULL,c("ID" , "JJ", "n", "age"  ,  "adjwt", "adjwtSE" ,"method","sex"))))
  print(paste("failure for pop model for yrl of",year ))
  }
  
  # compile yrl Res
  l.id <- droplevels(unique( dt1[dt1$age==1 & !is.na(dt1$pheno),"ID"]))
  Wtyrl <- data.frame(ID=l.id)
  Wtyrl <- merge(Wtyrl,EwtyrlMM,all.x = T)
  l.id <- droplevels(unique( Wtyrl[is.na(Wtyrl$adjwt),"ID"]))
  Wtyrl <- rbind(Wtyrl,Ewtyrlpop[Ewtyrlpop$ID %in%l.id,])
  Wtyrl <- Wtyrl[!is.na(Wtyrl$method),]
  
  
  # compile all age groupe sex of the pop
  Pop.adjWT <- rbind(WtMales,WtfeMales,Wtyrl,Wtlb)
  Pop.adjWT <- Pop.adjWT[!is.na(Pop.adjWT$adjwt),]
  Pop.adjWT$yr=year              # add year
  bigRes <- rbind(bigRes,Pop.adjWT)      #   put after other reasults
}


#  fix no estimation points    #######
library(lubridate)
data <- read.csv("Data/dailyMass_18_22.csv", sep =";", stringsAsFactors = F)
daily.E= 0.125/3  # specify nb of kg of error per days off (ref: Marco 1996)
id.yrz <- unique(data[,c("ID","yr")])
for(i in 1:nrow(data)){
  jj.15sep <- yday(ymd(paste0(data[i,"yr"],"-09-15")))
  jj.5jn <-  yday(ymd(paste0(data[i,"yr"],"-06-05")))
  jj.15jun <-  yday(ymd(paste0(data[i,"yr"],"-06-15")))
  
  if(abs(data[i,"JJ"]- jj.15sep )<=DJJ.eqal){   # was the measurement taken less than DJJ.eqal days around 15 sep
    wid <- which(bigRes$ID==data[i,"ID"] & bigRes$yr==data[i,"yr"]  & bigRes$JJ==114) #  was its mass adjusted
    if(length(wid)==0){bigRes <- rbind(bigRes,data.frame(ID=data[i,"ID"],JJ=114,      # if not,  use actual method
                                                         n=1,age=data[i,"age"],       # standard stuff
                                                         adjwt=data[i,"pheno"],      # use observed mass
                                                         adjwtSE=abs(data[i,"JJ"]- jj.15sep )*daily.E,method="actual",                # calculate SE using possible change in mass by days
                                                         sex=data[i,"sex"],yr=data[i,"yr"])) }
  }
  
  if(abs(data[i,"JJ"]- jj.5jn )<=DJJ.eqal & data[i,"age"]>0 ){            #  repeat for 5 jun
    wid <- which(bigRes$ID==data[i,"ID"] & bigRes$yr==data[i,"yr"]  & bigRes$JJ==12)
    if(length(wid)==0){bigRes <- rbind(bigRes,data.frame(ID=data[i,"ID"],JJ=12,
                                                         n=1,age=data[i,"age"],
                                                         adjwt=data[i,"pheno"],
                                                         adjwtSE=abs(data[i,"JJ"]- jj.5jn )*daily.E,method="actual",
                                                         sex=data[i,"sex"],yr=data[i,"yr"])) }
  }
  
  if(abs(data[i,"JJ"]- jj.15jun )<=DJJ.eqal & data[i,"age"]==0 ){       #  repeat for june15 for lambs
    wid <- which(bigRes$ID==data[i,"ID"] & bigRes$yr==data[i,"yr"]  & bigRes$JJ==22)
    if(length(wid)==0){bigRes <- rbind(bigRes,data.frame(ID=data[i,"ID"],JJ=22,
                                                         n=1,age=data[i,"age"],
                                                         adjwt=data[i,"pheno"],
                                                         adjwtSE=abs(data[i,"JJ"]- jj.15jun )*daily.E,method="actual",
                                                         sex=data[i,"sex"],yr=data[i,"yr"])) }
  }
  
}



###   print result
write_csv(bigRes, "Data/adjWT_2018_2022.csv")