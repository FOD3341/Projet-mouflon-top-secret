# Benjamin Larue
# Version 1          date: 2020-10-05
#  script pour obtenir les masse ajuste en juin et en septembre



###############################################################################.
#####     loading required libraries   ######
library(lme4)
library(plyr)
library(dplyr)
source("script/function.R")

# adjust POP --------------------------------------------------------------
bigRes <- data.frame()             # make empty data frame to house result
nsim=1000                         # set number of iteration for bootstrap used to get SE of mixed model
DJJ.eqal=1                         # if  observation id <= to DJJ.equal and missing, should use actual method
for(year in 1976:2023){           #  for all these years
  rm(list = c("d1", "data", "dt1", "EwtfMM", "EwtlbMM", "Ewtlbpop", "EwtMM",
              "Ewtoldfpop", "EwtoldMpop",
              "EwtyongMpop", "EwtyrlMM",
              "Ewtyrlpop",   "jun15", "jun5", "l.id",  "sept15",
              "WtfeMales", "Wtlb", "WtMales", "Wtyrl"))    # remove dis-used object
  
  data <- read.csv("Cleaned_data/cleaned_data.csv", sep=",")
  
  # load daily data
  data$ID <- as.factor(data$ID)
  data$Sex <- as.factor(data$Sex)
  
  # mutate hgl for pheno
  data <- data %>%
    rename(pheno = hlg)
  
  # mutate Year for year
  data <- data %>%
    rename(yr = Year)
  
  # mutate Year for year
  data <- data %>%
    rename(age = Age)
  
  # mutate Year for year
  data <- data %>%
    rename(sex = Sex)
  
  # year=1978
  data <- data[!is.na(data$pheno),]
  #
  dt1 <- data[data$yr == year,]            #  subset for year to analyse
  
  dt1$JJ.sqrt <- dt1$sq.date           #  calculate sqrt for adult models
  dt1 <- dt1[!is.nan(dt1$JJ.sqrt),]     # remove NAN; necessary when using Marco's date starting on 25 of may
  
  # getting Marco style julian. days
  d1 <- as.Date(paste(year,"/05/25",sep=""),format="%Y/%m/%d")       #
  # d1 <- as.Date(paste(year,"/01/01",sep=""),format="%Y/%m/%d")       # un-comment to use real JJ
  data$JJ <- data$date
  jun5 <- as.numeric(as.Date(paste(year,"/06/05",sep=""),format="%Y/%m/%d")-d1)+1
  jun5 <- as.numeric(as.Date(paste(year,"/06/05",sep=""),format="%Y/%m/%d")-d1)+1
  jun15 <- as.numeric(as.Date(paste(year,"/06/15",sep=""),format="%Y/%m/%d")-d1)+1
  sept15 <- as.numeric(as.Date(paste(year,"/09/15",sep=""),format="%Y/%m/%d")-d1)+1
  
  # run female
  EwtfMM <- try(feMale.mixed(dt1,nsim))
  if(!is.data.frame(EwtfMM)){EwtfMM <- as.data.frame(
    matrix(NA,nrow = 0,ncol = 7, dimnames = list(NULL,c("ID" , "JJ", "n", "age"  ,  "adjhlg", "adjhlgSE" ,"method"))))
  print(paste("failure for mixed model for females of",year ))
  }
  
  Ewtoldfpop <- try(feMale.pop.old(dt1,nsim,max.DeltaJ = 50))
  if(!is.data.frame(Ewtoldfpop)){Ewtoldfpop <- as.data.frame(
    matrix(NA,nrow = 0,ncol = 7, dimnames = list(NULL,c("ID" , "JJ", "n", "age"  ,  "adjhlg", "adjhlgSE" ,"method"))))
  print(paste("failure for pop model for old females of",year ))
  }
  
  
  # compile female Res
  l.id <- droplevels(unique( dt1[dt1$sex=="2"& !is.na(dt1$pheno),"ID"]))
  hlgfeMales <- data.frame(ID=l.id)
  hlgfeMales <- merge(hlgfeMales,EwtfMM,all.x = T)
  l.id <- droplevels(unique( hlgfeMales[is.na(hlgfeMales$adjhlg),"ID"]))
  hlgfeMales <- rbind(hlgfeMales,Ewtoldfpop[Ewtoldfpop$ID %in%l.id,])
  hlgfeMales <- hlgfeMales[!is.na(hlgfeMales$method),]
  if(nrow(hlgfeMales)>0) hlgfeMales$sex="female"
  
  # compile all age groupe sex of the pop
  Pop.adjhlg <- hlgfeMales
  Pop.adjhlg <- Pop.adjhlg[!is.na(Pop.adjhlg$adjhlg),]
  Pop.adjhlg$yr=year              # add year
  bigRes <- rbind(bigRes,Pop.adjhlg)      #   put after other reasults
}



#  fix no estimation points    #######
library(lubridate)
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
  
}



###   print result
write_csv(bigRes, "Cleaned_data/adjhlg.csv")
