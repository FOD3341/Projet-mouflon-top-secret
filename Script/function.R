library(lme4)
library(plyr)
library(readr)

###  males ####
# fct to estimate pheno with mixed model
Male.mixed <- function(dt1,nsim=nsim){
  tmp <- subset(dt1, dt1$sex=="male" & dt1$age>1 & !is.na(pheno))  # only use adult males
  l.id <- table(droplevels(tmp$ID))                                # count observation per id
  M <- lmer( pheno ~ JJ.sqrt + (1 + JJ.sqrt | ID), #+ (1 | yr),
             data = tmp[tmp$ID %in% names(l.id[l.id>1]),], na.action = na.omit, #     only used ID with > 1 capture
             control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))   # run model longuer to reduce non-convergeance
  #
  newdat <-  data.frame(JJ=rep(c((jun5),(sept15)),each=length(l.id[l.id>1])),      #   create a fake data.frame to generate prediction
                        JJ.sqrt=rep(c(sqrt(jun5),sqrt(sept15)),each=length(l.id[l.id>1])),
                        ID=c(names(l.id[l.id>1]),names(l.id[l.id>1])),
                        n=c(l.id[l.id>1],n=l.id[l.id>1]) )
  newdat <- merge(newdat,unique(tmp[,c("ID","age")]))

  # mySumm <- function(.) { getME(., "theta")}
  mySumm3 <- function(.) {predict(.,newdata=newdat,type="response")}             #  fonction to male prediction per ID
  boo01 <- bootMer(M, mySumm3, nsim = nsim)              # bootstrap the prediction fonction  nsim times
  # boo02 <- bootMer(M, mySumm, nsim = 100)
  PM <- cbind(newdat,adjwt=boo01$t0,adjwtSE=apply(boo01$t,2,sd),method="mixed")     #generate result data.frame from bootstrap result
  return(PM[,-3])
}


# old male pop
Male.pop.old <- function(dt1,nsim=nsim,max.DeltaJ=20){
  tmp <- subset(dt1, dt1$sex=="male" & dt1$age>2 & !is.na(pheno))    # for males 3 and older, old estimation methode with linear model
  tmp$ID <- droplevels(tmp$ID)
  l.id <- table(tmp$ID)                                             #   count observation per id
  l.fittable.age <- table(tmp$age,tmp$ID)                           #   count observation by id per age
  pop.V <- data.frame()
  for(i.age in as.numeric(names(rowSums(l.fittable.age>=2)[rowSums(l.fittable.age>=2)>=1]))){  # for age for which at least 1 ID as more than 2 observation
    lm <- ldply(colnames(l.fittable.age)[l.fittable.age[rownames(l.fittable.age)==i.age ,] >=2 ],function(i.id){  # for id with > than 2 capture
      lmt <- lm(pheno ~ JJ.sqrt,data = tmp[tmp$ID==i.id,], na.action = na.omit) #                                  fit lm
      pred <- predict(lmt,newdata = data.frame(JJ.sqrt=c(sqrt(jun5),sqrt(sept15))),se.fit = T)     #              make prediction
      res <- data.frame(i.id,summary(lmt)$coefficients[2,1],summary(lmt)$coefficients[2,2],length(lmt$residuals),
                        wt12=pred$fit[1],wt12.se=pred$se.fit[1],wt114=pred$fit[2],wt114.se=pred$se.fit[2])     # generate prediction dataframe
      names(res) <- c("id","Beta","SE","N","wt12","wt12.se","wt114","wt114.se")
      return(res)
    })

    # calculate mean individual growth  and error on growth
    bhat <- list(bhat=mean(lm$Beta))
    if(nrow(lm)==1){bhat$SE <- lm[,"SE"]}   # if only one regresion, use it's SE
    if(nrow(lm)>1){   # if more than 1 regression, use mean SE of the lot; if no SE is available, use sd of slopes as SE
      lm[is.nan(lm[,"SE"]),"SE"] <- NA
      SE <- c(0,mean(lm[,"SE"],na.rm = T))
      SE <- SE[!is.nan(SE)]
      bhat$SE <- max(SE)
      if(bhat$SE==0){bhat$SE <- sd(lm$Beta)}
    }


    #   calculate adjusted pheno or extract from lm if present
    age.res <- data.frame()
    for(i.id in colnames(l.fittable.age)[l.fittable.age[rownames(l.fittable.age)==i.age,]>0]){   # for Id of age==i.age
      if(l.fittable.age[as.character(i.age),i.id]>1){   #                                           if there is more than 1 point
        tmp.res <- data.frame(ID=i.id,JJ=c(12,114),n=l.fittable.age[as.character(i.age),i.id],age=i.age,
                              adjwt= as.numeric(lm[lm$id==i.id,c("wt12","wt114")]),       #           extract estimate from lm prediction object
                              adjwtSE= as.numeric(lm[lm$id==i.id,c("wt12.se","wt114.se")] ),      #    extract SE
                              method="pop")
        laply(c(12,114),function(x)min(abs(x-tmp[tmp$ID==i.id,"JJ"])))<50
        tmp.res$adjwt <- tmp.res$adjwt/as.numeric( laply(c(12,114),function(x)min(abs(x-tmp[tmp$ID==i.id,"JJ"])))<50)  # devides by 0 ( results in infinite) if measurement date is more than 50 days from adjustement date
        tmp.res[is.infinite(tmp.res$adjwt),c("adjwt","adjwtSE")] <- NA                  #  changes previously created infinites to NA
        age.res <- rbind(age.res, tmp.res)
      }
      if(l.fittable.age[as.character(i.age),i.id]==1){ #                                           if there is only 1 point
        DeltaJJ <- sqrt(c(12,114))-tmp[tmp$ID==i.id,"JJ.sqrt"]      #                              calculate diff between measurement and estimation date (DeltaJJ)
        tmp.res <- data.frame(ID=i.id,JJ=c(12,114),n=l.fittable.age[as.character(i.age),i.id],age=i.age,
                              adjwt= tmp[tmp$ID==i.id,"pheno"]+bhat$bhat*DeltaJJ,         #        calculate adjusted wt using DeltaJJ and population mean age specific growth rate
                              adjwtSE= abs(bhat$SE*DeltaJJ) ,method="pop"   )             #       calculate SE by multiplying DeltaJJ and SE; the farther adhusted date are from measurement, the higher the SE is
        tmp.res$adjwt <- tmp.res$adjwt/as.numeric( laply(c(12,114),function(x)min(abs(x-tmp[tmp$ID==i.id,"JJ"])))<max.DeltaJ)  # devides by 0 ( results in infinite) if measurement date is more than maxdelta.J days from adjustement date
        tmp.res[is.infinite(tmp.res$adjwt),c("adjwt","adjwtSE")] <- NA                  #  changes previously created infinites to NA
        age.res <- rbind(age.res, tmp.res)
      }
    }

    pop.V <- rbind(pop.V,age.res)
  }
  return(pop.V)
}

# young male pop
Male.pop.young <- function(dt1,nsim=nsim,max.DeltaJ=6){
  tmp <- subset(dt1, dt1$sex=="male" & dt1$age==2 & !is.na(pheno))
  tmp$ID <- droplevels(tmp$ID)
  l.id <- table(tmp$ID)
  l.fittable.age <- table(tmp$age,tmp$ID)
  i.age=2
  lm <- ldply(colnames(l.fittable.age)[l.fittable.age[rownames(l.fittable.age)==i.age ,] >=2 ],function(i.id){
    lmt <- lm(pheno ~ JJ.sqrt,data = tmp[tmp$ID==i.id,], na.action = na.omit)
    pred <- predict(lmt,newdata = data.frame(JJ.sqrt=c(sqrt(jun5),sqrt(sept15))),se.fit = T)
    res <- data.frame(i.id,summary(lmt)$coefficients[2,1],summary(lmt)$coefficients[2,2],length(lmt$residuals),
                      wt12=pred$fit[1],wt12.se=pred$se.fit[1],wt114=pred$fit[2],wt114.se=pred$se.fit[2])
    names(res) <- c("id","Beta","SE","N","wt12","wt12.se","wt114","wt114.se")
    return(res)
  })
  # calculate mean individual growth
  bhat <- list(bhat=mean(lm$Beta))
  if(nrow(lm)==1){bhat$SE <- lm[,"SE"]}   # if only one regresion, use it's SE
  if(nrow(lm)>1){   # if more than 1 regression, use max SE of the lot; if no SE is available, use slope as SE
    lm[is.nan(lm[,"SE"]),"SE"] <- NA
    SE <- c(0,mean(lm[,"SE"],na.rm = T))
    SE <- SE[!is.nan(SE)]
    bhat$SE <- max(SE)
    if(bhat$SE==0){bhat$SE <- sd(lm$Beta)}
  }
  #   calculate adjusted pheno or extract from lm if present
  age.res <- data.frame()#ID=NA,JJ=NA,n=NA,"age,"adjwt","adjwtSE",method="pop")
  for(i.id in colnames(l.fittable.age)[l.fittable.age[rownames(l.fittable.age)==i.age,]>0]){
    if(l.fittable.age[as.character(i.age),i.id]>1){
      tmp.res <- data.frame(ID=i.id,JJ=c(12,114),n=l.fittable.age[as.character(i.age),i.id],age=i.age,
                            adjwt= as.numeric(lm[lm$id==i.id,c("wt12","wt114")]),
                            adjwtSE= as.numeric(lm[lm$id==i.id,c("wt12.se","wt114.se")] ),
                            method="pop")
      tmp.res$adjwt <- tmp.res$adjwt/as.numeric( laply(c(12,114),function(x)min(abs(x-tmp[tmp$ID==i.id,"JJ"])))<50)  # devides by 0 ( results in infinite) if measurement date is more than 50 days from adjustement date
      tmp.res[is.infinite(tmp.res$adjwt),c("adjwt","adjwtSE")] <- NA
      age.res <- rbind(age.res, tmp.res)
    }
    if(l.fittable.age[as.character(i.age),i.id]==1){
      DeltaJJ <- sqrt(c(12,114))-tmp[tmp$ID==i.id,"JJ.sqrt"]
      tmp.res <- data.frame(ID=i.id,JJ=c(12,114),n=l.fittable.age[as.character(i.age),i.id],age=i.age,
                            adjwt= tmp[tmp$ID==i.id,"pheno"]+bhat$bhat*DeltaJJ,
                            adjwtSE= abs(bhat$SE*DeltaJJ) ,method="pop"   )
      tmp.res$adjwt <- tmp.res$adjwt/as.numeric( laply(c(12,114),function(x)min(abs(x-tmp[tmp$ID==i.id,"JJ"])))<max.DeltaJ)  # devides by 0 ( results in infinite) if measurement date is more than 50 days from adjustement date
      tmp.res[is.infinite(tmp.res$adjwt),c("adjwt","adjwtSE")] <- NA
      age.res <- rbind(age.res, tmp.res)
    }
  }

  #     pop.V <- rbind(pop.V,age.res)
  #   }
  return(age.res)
}



# females   ######
feMale.mixed <- function(dt1,nsim=nsim){
  tmp <- subset(dt1, dt1$sex=="female" & dt1$age>1 & !is.na(pheno))
  l.id <- table(droplevels(tmp$ID))
  M <- lmer( pheno ~ JJ.sqrt + (1 + JJ.sqrt | ID), #+ (1 | yr),
             data = tmp[tmp$ID %in% names(l.id[l.id>1]),], na.action = na.omit,
             control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
  # Mz$males <- M
  newdat <-  data.frame(JJ=rep(c((jun5),(sept15)),each=length(l.id[l.id>1])),
                        JJ.sqrt=rep(c(sqrt(jun5),sqrt(sept15)),each=length(l.id[l.id>1])),
                        ID=c(names(l.id[l.id>1]),names(l.id[l.id>1])),
                        n=c(l.id[l.id>1],n=l.id[l.id>1]) )
  newdat <- merge(newdat,unique(tmp[,c("ID","age")]))

  # mySumm <- function(.) { getME(., "theta")}
  mySumm3 <- function(.) {predict(.,newdata=newdat,type="response")}
  boo01 <- bootMer(M, mySumm3, nsim = nsim)
  # boo02 <- bootMer(M, mySumm, nsim = 100)
  PM <- cbind(newdat,adjwt=boo01$t0,adjwtSE=apply(boo01$t,2,sd),method="mixed")
  return(PM[,-3])
}

# old female pop    same as male.yuong, but for female
feMale.pop.old <- function(dt1,nsim=nsim,max.DeltaJ=6){
  tmp <- subset(dt1, dt1$sex=="female" & dt1$age>=2 & !is.na(pheno))
  tmp$ID <- droplevels(tmp$ID)
  l.id <- table(tmp$ID)
  l.fittable.age <- table(tmp$age,tmp$ID)
  pop.V <- data.frame()
  for(i.age in as.numeric(names(rowSums(l.fittable.age>=2)[rowSums(l.fittable.age>=2)>=1]))){
    lm <- ldply(colnames(l.fittable.age)[l.fittable.age[rownames(l.fittable.age)==i.age ,] >=2 ],function(i.id){
      lmt <- lm(pheno ~ JJ.sqrt,data = tmp[tmp$ID==i.id,], na.action = na.omit)
      pred <- predict(lmt,newdata = data.frame(JJ.sqrt=c(sqrt(jun5),sqrt(sept15))),se.fit = T)
      res <- data.frame(i.id,summary(lmt)$coefficients[2,1],summary(lmt)$coefficients[2,2],length(lmt$residuals),
                        wt12=pred$fit[1],wt12.se=pred$se.fit[1],wt114=pred$fit[2],wt114.se=pred$se.fit[2])
      names(res) <- c("id","Beta","SE","N","wt12","wt12.se","wt114","wt114.se")
      return(res)
    })
    # calculate mean individual growth
    bhat <- list(bhat=mean(lm$Beta))
    if(nrow(lm)==1){bhat$SE <- lm[,"SE"]}   # if only one regresion, use it's SE
    if(nrow(lm)>1){   # if more than 1 regression, use max SE of the lot; if no SE is available, use slope as SE
      lm[is.nan(lm[,"SE"]),"SE"] <- NA
      SE <- c(0,mean(lm[,"SE"],na.rm = T))
      SE <- SE[!is.nan(SE)]
      bhat$SE <- max(SE)
      if(bhat$SE==0){bhat$SE <- sd(lm$Beta)}
    }
    #   calculate adjusted pheno or extract from lm if present
    age.res <- data.frame()#ID=NA,JJ=NA,n=NA,"age,"adjwt","adjwtSE",method="pop")
    for(i.id in colnames(l.fittable.age)[l.fittable.age[rownames(l.fittable.age)==i.age,]>0]){
      if(l.fittable.age[as.character(i.age),i.id]>1){
        tmp.res <- data.frame(ID=i.id,JJ=c(12,114),n=l.fittable.age[as.character(i.age),i.id],age=i.age,
                              adjwt= as.numeric(lm[lm$id==i.id,c("wt12","wt114")]),
                              adjwtSE= as.numeric(lm[lm$id==i.id,c("wt12.se","wt114.se")] ),
                              method="pop")
        tmp.res$adjwt <- tmp.res$adjwt/as.numeric( laply(c(12,114),function(x)min(abs(x-tmp[tmp$ID==i.id,"JJ"])))<50)  # devides by 0 ( results in infinite) if measurement date is more than 50 days from adjustement date
        tmp.res[is.infinite(tmp.res$adjwt),c("adjwt","adjwtSE")] <- NA
        age.res <- rbind(age.res, tmp.res)
      }
      if(l.fittable.age[as.character(i.age),i.id]==1){
        DeltaJJ <- sqrt(c(12,114))-tmp[tmp$ID==i.id,"JJ.sqrt"]
        tmp.res <- data.frame(ID=i.id,JJ=c(12,114),n=l.fittable.age[as.character(i.age),i.id],age=i.age,
                              adjwt= tmp[tmp$ID==i.id,"pheno"]+bhat$bhat*DeltaJJ,
                              adjwtSE= abs(bhat$SE*DeltaJJ) ,method="pop"   )
        tmp.res$adjwt <- tmp.res$adjwt/as.numeric( laply(c(12,114),function(x)min(abs(x-tmp[tmp$ID==i.id,"JJ"])))<max.DeltaJ)  # devides by 0 ( results in infinite) if measurement date is more than 50 days from adjustement date
        tmp.res[is.infinite(tmp.res$adjwt),c("adjwt","adjwtSE")] <- NA
        age.res <- rbind(age.res, tmp.res)
      }
    }

    pop.V <- rbind(pop.V,age.res)
  }
  return(pop.V)
}



###lambs   #####
lb.mixed <- function(dt1,nsim=nsim){
  tmp <- subset(dt1,  dt1$age==0 & !is.na(pheno))    #      select lambs and remove lines with no phenotipe info
  l.id <- table(droplevels(tmp$ID))                  #      counts observation per id
  M <- lmer( pheno ~ JJ *sex + (1 + JJ | ID), #+ (1 | yr),    3 fit mixed model with linear effect of JJ and interaction with sex; random growth linear slope
             data = tmp[tmp$ID %in% names(l.id[l.id>1]),], na.action = na.omit,       #  only use id with more than one observation
             control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))    #  reduce chances of  non-convergeance
  # Mz$males <- M
  newdat <-  data.frame(JJ=rep(c((jun15),(sept15)),each=length(l.id[l.id>1])),     #  make new data for prediction
                        ID=c(names(l.id[l.id>1]),names(l.id[l.id>1])),
                        n=c(l.id[l.id>1],n=l.id[l.id>1]) )
  newdat <- merge(newdat,unique(tmp[,c("ID","sex")]))

  # mySumm <- function(.) { getME(., "theta")}
  mySumm3 <- function(.) {predict(.,newdata=newdat,type="response")}
  boo01 <- bootMer(M, mySumm3, nsim = nsim)                                  # bootstrap prediction fonction
  # boo02 <- bootMer(M, mySumm, nsim = 100)
  PM <- cbind(newdat,adjwt=boo01$t0,adjwtSE=apply(boo01$t,2,sd),method="mixed")   # extract predicted adjusted value and error
  PM$age <- 0
  return(PM[,c(1,2,3,8,5:7,4)])                                                #  return reasults
}


lb.pop <- function(dt1,nsim=nsim,max.DeltaJ=3){
  tmp <- subset(dt1, dt1$age==0 & !is.na(pheno))                                #  select lambs
  tmp$ID <- droplevels(tmp$ID)                                                  # drops un-used factors level
  l.id <- table(tmp$ID)                                                         # count observation per id
  l.fittable.age <- table(tmp$sex,tmp$ID)                                       # count observation per id by sex
  pop.V <- data.frame()                                                         # make empty data.frame that will house result
  for(i.sex in names(rowSums(l.fittable.age>=2)[rowSums(l.fittable.age>=2)>=1])){      #  for each sex
    lm <- ldply(colnames(l.fittable.age)[l.fittable.age[rownames(l.fittable.age)==i.sex ,] >=2 ],function(i.id){   #for id of sex==i.sex with 2 or more observation
      lmt <- lm(pheno ~ JJ,data = tmp[tmp$ID==i.id,], na.action = na.omit)                      # fit model
      pred <- predict(lmt,newdata = data.frame(JJ=c(jun15,sept15)),se.fit = T)                 #  make prediction
      res <- data.frame(i.id,summary(lmt)$coefficients[2,1],summary(lmt)$coefficients[2,2],length(lmt$residuals), #  extract prediction from models for each i.id tested
                        wt12=pred$fit[1],wt12.se=pred$se.fit[1],wt114=pred$fit[2],wt114.se=pred$se.fit[2])
      names(res) <- c("id","Beta","SE","N","wt12","wt12.se","wt114","wt114.se")
      return(res)
    })

    # calculate mean of  individual growths
    bhat <- list(bhat=mean(lm$Beta))
    if(nrow(lm)==1){bhat$SE <- lm[,"SE"]}   # if only one regresion, use it's SE
    if(nrow(lm)>1){   # if more than 1 regression, use max SE of the lot; if no SE is available, use slope as SE
      lm[is.nan(lm[,"SE"]),"SE"] <- NA
      SE <- c(0,mean(lm[,"SE"],na.rm = T))
      SE <- SE[!is.nan(SE)]
      bhat$SE <- max(SE)
      if(bhat$SE==0){bhat$SE <- sd(lm$Beta)}
    }
    #   calculate adjusted pheno or extract from lm if present
    age.res <- data.frame()#ID=NA,JJ=NA,n=NA,"age,"adjwt","adjwtSE",method="pop")
    for(i.id in colnames(l.fittable.age)[l.fittable.age[rownames(l.fittable.age)==i.sex,]>0]){  # for i.id of a certain sex that where observed
      if(l.fittable.age[as.character(i.sex),i.id]>1){                    # if the id was observed more than once, get prediction from lm object
        tmp.res <- data.frame(ID=i.id,JJ=c(jun15,114),n=l.fittable.age[as.character(i.sex),i.id],age=0,
                              adjwt= as.numeric(lm[lm$id==i.id,c("wt12","wt114")]),
                              adjwtSE= as.numeric(lm[lm$id==i.id,c("wt12.se","wt114.se")] ),
                              method="pop",sex=i.sex)
        tmp.res$adjwt <- tmp.res$adjwt/as.numeric( laply(c(12,114),function(x)min(abs(x-tmp[tmp$ID==i.id,"JJ"])))<50)  # devides by 0 ( results in infinite) if measurement date is more than 50 days from adjustement date
        tmp.res[is.infinite(tmp.res$adjwt),c("adjwt","adjwtSE")] <- NA             #  remove infinite and put NA
        age.res <- rbind(age.res, tmp.res)
      }
      if(l.fittable.age[as.character(i.sex),i.id]==1){                # if i.id was observed only once
        DeltaJJ <- c(jun15,114)-tmp[tmp$ID==i.id,"JJ"]              #  calculate difference with 15 june and 15 septembre
        tmp.res <- data.frame(ID=i.id,JJ=c(jun15,114),n=l.fittable.age[as.character(i.sex),i.id],age=0,
                              adjwt= tmp[tmp$ID==i.id,"pheno"]+bhat$bhat*DeltaJJ,         # observed measure + time diffenrece multiplied by growth
                              adjwtSE= abs(bhat$SE*DeltaJJ) ,method="pop"  ,sex=i.sex )   # error multiplied by time difference
        tmp.res$adjwt <- tmp.res$adjwt/as.numeric( laply(c(12,114),function(x)min(abs(x-tmp[tmp$ID==i.id,"JJ"])))<max.DeltaJ)  # devides by 0 ( results in infinite) if measurement date is more than max.DeltaJ days from adjustement date
        tmp.res[is.infinite(tmp.res$adjwt),c("adjwt","adjwtSE")] <- NA
        age.res <- rbind(age.res, tmp.res)
      }
    }

    pop.V <- rbind(pop.V,age.res)
  }
  return(pop.V)
}



###   yearlingg  #####
yrl.mixed <- function(dt1,nsim=nsim){
  tmp <- subset(dt1,  dt1$age==1 & !is.na(pheno))
  l.id <- table(droplevels(tmp$ID))
  M <- lmer( pheno ~ JJ *sex + (1 + JJ | ID), #+ (1 | yr),
             data = tmp[tmp$ID %in% names(l.id[l.id>1]),], na.action = na.omit,
             control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
  # Mz$males <- M
  newdat <-  data.frame(JJ=rep(c((jun5),(sept15)),each=length(l.id[l.id>1])),
                        ID=c(names(l.id[l.id>1]),names(l.id[l.id>1])),
                        n=c(l.id[l.id>1],n=l.id[l.id>1]) )
  newdat <- merge(newdat,unique(tmp[,c("ID","sex")]))

  # mySumm <- function(.) { getME(., "theta")}
  mySumm3 <- function(.) {predict(.,newdata=newdat,type="response")}
  boo01 <- bootMer(M, mySumm3, nsim = nsim)
  # boo02 <- bootMer(M, mySumm, nsim = 100)
  PM <- cbind(newdat,adjwt=boo01$t0,adjwtSE=apply(boo01$t,2,sd),method="mixed")
  PM$age <- 1
  return(PM[,c(1,2,3,8,5:7,4)])
}


yrl.pop <- function(dt1,nsim=nsim,max.DeltaJ=3){
  tmp <- subset(dt1, dt1$age==1 & !is.na(pheno))
  tmp$ID <- droplevels(tmp$ID)
  l.id <- table(tmp$ID)
  l.fittable.age <- table(tmp$sex,tmp$ID)
  pop.V <- data.frame()
  for(i.sex in names(rowSums(l.fittable.age>=2)[rowSums(l.fittable.age>=2)>=1])){
    lm <- ldply(colnames(l.fittable.age)[l.fittable.age[rownames(l.fittable.age)==i.sex ,] >=2 ],function(i.id){
      lmt <- lm(pheno ~ JJ,data = tmp[tmp$ID==i.id,], na.action = na.omit)
      pred <- predict(lmt,newdata = data.frame(JJ=c(jun5,sept15)),se.fit = T)
      res <- data.frame(i.id,summary(lmt)$coefficients[2,1],summary(lmt)$coefficients[2,2],length(lmt$residuals),
                        wt12=pred$fit[1],wt12.se=pred$se.fit[1],wt114=pred$fit[2],wt114.se=pred$se.fit[2])
      names(res) <- c("id","Beta","SE","N","wt12","wt12.se","wt114","wt114.se")
      return(res)
    })
    # calculate mean individual growth
    bhat <- list(bhat=mean(lm$Beta))
    if(nrow(lm)==1){bhat$SE <- lm[,"SE"]}   # if only one regresion, use it's SE
    if(nrow(lm)>1){   # if more than 1 regression, use max SE of the lot; if no SE is available, use slope as SE
      lm[is.nan(lm[,"SE"]),"SE"] <- NA
      SE <- c(0,mean(lm[,"SE"],na.rm = T))
      SE <- SE[!is.nan(SE)]
      bhat$SE <- max(SE)
      if(bhat$SE==0){bhat$SE <- sd(lm$Beta)}
    }
    #   calculate adjusted pheno or extract from lm if present
    age.res <- data.frame()#ID=NA,JJ=NA,n=NA,"age,"adjwt","adjwtSE",method="pop")
    for(i.id in colnames(l.fittable.age)[l.fittable.age[rownames(l.fittable.age)==i.sex,]>0]){
      if(l.fittable.age[as.character(i.sex),i.id]>1){
        tmp.res <- data.frame(ID=i.id,JJ=c(12,114),n=l.fittable.age[as.character(i.sex),i.id],age=1,
                              adjwt= as.numeric(lm[lm$id==i.id,c("wt12","wt114")]),
                              adjwtSE= as.numeric(lm[lm$id==i.id,c("wt12.se","wt114.se")] ),
                              method="pop",sex=i.sex)
        tmp.res$adjwt <- tmp.res$adjwt/as.numeric( laply(c(12,114),function(x)min(abs(x-tmp[tmp$ID==i.id,"JJ"])))<50)  # devides by 0 ( results in infinite) if measurement date is more than 50 days from adjustement date
        tmp.res[is.infinite(tmp.res$adjwt),c("adjwt","adjwtSE")] <- NA
        age.res <- rbind(age.res, tmp.res)
      }
      if(l.fittable.age[as.character(i.sex),i.id]==1){
        DeltaJJ <- c(12,114)-tmp[tmp$ID==i.id,"JJ"]
        tmp.res <- data.frame(ID=i.id,JJ=c(12,114),n=l.fittable.age[as.character(i.sex),i.id],age=1,
                              adjwt= tmp[tmp$ID==i.id,"pheno"]+bhat$bhat*DeltaJJ,
                              adjwtSE= abs(bhat$SE*DeltaJJ) ,method="pop"  ,sex=i.sex )
        tmp.res$adjwt <- tmp.res$adjwt/as.numeric( laply(c(12,114),function(x)min(abs(x-tmp[tmp$ID==i.id,"JJ"])))<max.DeltaJ)  # devides by 0 ( results in infinite) if measurement date is more than 50 days from adjustement date
        tmp.res[is.infinite(tmp.res$adjwt),c("adjwt","adjwtSE")] <- NA
        age.res <- rbind(age.res, tmp.res)
      }
    }

    pop.V <- rbind(pop.V,age.res)
  }
  return(pop.V)
}

