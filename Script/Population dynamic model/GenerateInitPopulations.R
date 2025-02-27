#######################################################################################################################################
########################################                                                ################################################
########################################                Simulation Model                #################################################    
########################################         Create initial metapopulations         ################################################
########################################                                                ################################################
#######################################################################################################################################

####### 
####### Note: Create metapopulations gathering initial protected and exploited populations
####### 

rm(list=ls()) 
## !!! NB. Below put your folder address instead of "[...]" !!!
setwd("[...]/SimulationModel/")

#######################################################################################################################################
########################################         Load of useful resources               #################################################    
#######################################################################################################################################

####### Librairies ########
if (!require("tidyr")) install.packages("tidyr") # file extraction
if (!require("RColorBrewer")) install.packages("RColorBrewer") # colors for plot
if (!require("latex2exp")) install.packages("latex2exp") # latex reading for greek letters
if (!require("openxlsx")) install.packages("openxlsx") # file extraction
if (!require("tidyverse")) install.packages("tidyverse") # data manipulation and plotting with ggplot2
if (!require("plyr")) install.packages("plyr") # data frames manipulation
if (!require("effects")) install.packages("effects") # plotting parameter effects
if (!require("dplyr")) install.packages("dplyr") # set of tools that solves a common set of problems
if (!require("ggthemes")) install.packages("ggthemes") # set of tools that solves a common set of problems
if (!require("modelr")) install.packages("modelr") # data manipulation and plotting with ggplot2
if (!require("tidybayes")) install.packages("tidybayes") # data manipulation and plotting with ggplot2
if (!require("lubridate")) install.packages("lubridate") # dates and time tool
if (!require("forecast")) install.packages("forecast") # dates and time tool
if (!require("janitor")) install.packages("janitor") # data frame proportion tables

######## Simulation Functions ########
source("./Scripts/Run_Model/SimulationFunctions.R")

#######################################################################################################################################
########################################         Creation of initial metapopulations         #################################################    
#######################################################################################################################################

Nb.Init.pops=1 # number of initial metapopulation populations to generate

for(p in 1:Nb.Init.pops){
  
  ######## Charcateristics to create ########
  #### Individuals variables
  # id:: unique identification nb.
  # a:: age (in time step)
  # s:: sex ("F": female, "M": male)
  # o:: origin ("P": protected, "H": harvested)
  # m:: migration status ("TRUE"/"FALSE": if migrate during current year)
  # r:: reproductive status ("TRUE"/"FALSE": repro. success during current mating season, so that females have 1 lamb/year)
  # n:: total record of reproductions
  # y:: horn size (in cm)
  # y_inf:: asymptotic horn length (in cm)
  # bv_f:: breeding value for female horn length (in cm)
  # bv_m:: breeding value for male horn length (in cm)
  
  
  #### Populational variables
  # Y:: sum of mature males horn length at each time step
  # y_mean_females/males:: populational mean horn length by sex at each time step
  # H:: nb. of males killed during time step t
  # q:: proportion of legal rams at each time step
  # N/NH/NP:: meta/harvested/protected population size at each time step
  # Nm/f:: number of mature males/females in population of origin at each time step
  
  ######## Initial metapopulation ########
  N_H=6500 # 6500/10
  N_P=4500 # 4500/10
  N=N_H+N_P
  
  ##### Estimated parameters for horn growth
  y_inf_females=25.69 # [25.23;26.15] # y_inf_females: female asymptotic horn length
  y_inf_males=102.51 # [99.89;105.25] # y_inf_females: male asymptotic horn length
  rho_females=0.60 # [0.57;0.63] # rho_females: female growth rate
  rho_males=0.22 # [0.21;0.23] # rho_males: male growth rate
  
  ##### Parameters for sex allocation
  sigma=0.5 # sigma: sex ratio
  
  ##### Parameters for quantitative genetics
  hf=0.223 # [0.090;0.446] # hf: female horn length heritability
  hm=0.397 # [0.203;0.534] # hm: male horn length heritability
  G_mm=17.884 # [9.82–25.881] # G_mm: intra-sexual genetic variance for male horn length
  G_mf=5.345 # [1.928–8.144] # G_mf: inter-sexual genetic variance horn length
  G_ff=1.622 # [0.748–3.963] # G_ff: intra-sexual genetic variance for female horn length
  E_ff=round(G_ff/hf-G_ff,digits=2) # E_ff: intra-sexual environmental variance for female horn length (VP=VA/h2;VE=VP-VA=VA/h2-VA)
  E_mm=round(G_mm/hm-G_mm,digits=2) # E_mm: intra-sexual environmental variance for male horn length (VP=VA/h2;VE=VP-VA=VA/h2-VA)
  
  ##### Function to compute initial horn length from initial age and asymptotic horn length
  # cf. von Bertalanffy growth equation with estimated parameters
  AsymptHornLgth.And.Age.To.HornLgth=function(a,s,y_inf){
    rho=ifelse((s=="F"),rho_females,rho_males)
    y=round(y_inf*(1-exp(-rho*a)),digits=2)
    return(y)
  }
  
  
  #################### Initial data frame #################### 
  
  ########## Harvested population
  ##### Creation of initial data frame
  init.id=rep(NA,N_H)
  init.a=rep(NA,N_H)
  init.s=rep(NA,N_H)
  init.o=rep(NA,N_H)
  init.m=rep(NA,N_H)
  init.r=rep(NA,N_H)
  init.y=rep(NA,N_H)
  init.y_inf=rep(NA,N_H)
  init.bv_f=rep(NA,N_H)
  init.bv_m=rep(NA,N_H)
  
  ##### Creation of initial data frame
  init_base_H=data.frame(id=init.id,a=init.a,s=init.s,o=init.o,m=init.m,r=init.r,y=init.y,y_inf=init.y_inf,bv_f=init.bv_m,bv_m=init.bv_m)
  
  ##### Assignation of initial variables
  # init_base_H$a=sample(1:12,N,replace=TRUE) # random pop.age in [1,12]
  init_base_H$a=round(rnorm(N_H,5,2),digits=0) # normal distribtuion of ages centered in 6 years, with 3 years of variance
  init_base_H$a=ifelse(init_base_H$a<1,1,init_base_H$a) # individuals of age at least 1 year
  init_base_H$s=ifelse(rbinom(N_H,1,sigma)==1,"F","M") # sex allocation with default sex ratio of 0.5 
  init_base_H$o="H" # attribution of pop. of origin with equal probability
  init_base_H$m=rep(FALSE,N_H) # no migrants
  init_base_H$r=rep(0,N_H) # no offspring
  init_base_H$y_inf[which(init_base_H$s=="F")]=round(rnorm(length(init_base_H$y_inf[which(init_base_H$s=="F")]),y_inf_females,2.5),digits=2) # random deviate of normal distribution with mean y_inf_females and standard deviation of 2.5 cm
  init_base_H$y_inf[which(init_base_H$s=="M")]=round(rnorm(length(init_base_H$y_inf[which(init_base_H$s=="M")]),y_inf_males,2.5),digits=2) # random deviate of normal distribution with mean y_inf_males and standard deviation of 2.5 cm
  
  ## Initialization of breeding values
  y_inf_mean_f=mean(init_base_H$y_inf[which(init_base_H$s=="F")])
  y_inf_mean_m=mean(init_base_H$y_inf[which(init_base_H$s=="M")])
  
  Nf=nrow(subset(init_base_H,init_base_H$s=="F"))
  Nm=nrow(subset(init_base_H,init_base_H$s=="M"))
  init_base_H$bv_f[which(init_base_H$s=="F")]=round((init_base_H$y_inf[which(init_base_H$s=="F")]-y_inf_mean_f)+rnorm(rep(0,Nf),0,2),digits=2)
  init_base_H$bv_m[which(init_base_H$s=="F")]=round((init_base_H$y_inf[which(init_base_H$s=="F")]-y_inf_mean_f)+rnorm(rep(0,Nf),0,10),digits=2)
  init_base_H$bv_f[which(init_base_H$s=="M")]=round((init_base_H$y_inf[which(init_base_H$s=="M")]-y_inf_mean_m)+rnorm(rep(0,Nm),0,10),digits=2)
  init_base_H$bv_m[which(init_base_H$s=="M")]=round((init_base_H$y_inf[which(init_base_H$s=="M")]-y_inf_mean_m)+rnorm(rep(0,Nm),0,2),digits=2)
  
  ##### Figures
  PLOT<-FALSE
  if (PLOT){
    ## Figures of bv for females
    # -> positive relation betwwen female bv and female asymptotic horn length
    # -> no relation betwwen female bv and male bv, nor between male bv and female asymptotic horn length
    plot(init_base_H$bv_f,init_base_H$bv_m)
    plot(init_base_H$bv_f[which(init_base_H$s=="F")],init_base_H$bv_m[which(init_base_H$s=="F")])
    plot(init_base_H$bv_f[which(init_base_H$s=="F")],init_base_H$y_inf[which(init_base_H$s=="F")])
    plot(init_base_H$bv_m[which(init_base_H$s=="F")],init_base_H$y_inf[which(init_base_H$s=="F")])
    
    ## Figures of bv for males
    # -> positive relation betwwen male bv and male asymptotic horn length
    # -> no relation betwwen male bv and female bv, nor between female bv and male asymptotic horn length
    plot(init_base_H$bv_f[which(init_base_H$s=="M")],init_base_H$bv_m[which(init_base_H$s=="M")])
    plot(init_base_H$bv_m[which(init_base_H$s=="M")],init_base_H$y_inf[which(init_base_H$s=="M")])
    plot(init_base_H$bv_f[which(init_base_H$s=="M")],init_base_H$y_inf[which(init_base_H$s=="M")])
  }
  
  ########## Protected population
  init.id=rep(NA,N_P)
  init.a=rep(NA,N_P)
  init.s=rep(NA,N_P)
  init.o=rep(NA,N_P)
  init.m=rep(NA,N_P)
  init.r=rep(NA,N_P)
  init.y=rep(NA,N_P)
  init.y_inf=rep(NA,N_P)
  init.bv_f=rep(NA,N_P)
  init.bv_m=rep(NA,N_P)
  
  ##### Creation of initial data frame
  init_base_P=data.frame(id=init.id,a=init.a,s=init.s,o=init.o,m=init.m,r=init.r,y=init.y,y_inf=init.y_inf,bv_f=init.bv_m,bv_m=init.bv_m)
  
  ##### Assignation of initial variables
  # init_base_P$a=sample(1:12,N,replace=TRUE) # random pop.age in [1,12]
  init_base_P$a=round(rnorm(N_P,5,2),digits=0) # normal distribtuion of ages centered in 6 years, with 3 years of variance
  init_base_P$a=ifelse(init_base_P$a<1,1,init_base_P$a) # individuals of age at least 1 year
  init_base_P$s=ifelse(rbinom(N_P,1,sigma)==1,"F","M") # sex allocation with default sex ratio of 0.5 
  init_base_P$o="P" # attribution of pop. of origin with equal probability
  init_base_P$m=rep(FALSE,N_P) # no migrants
  init_base_P$r=rep(FALSE,N_P) # no offspring
  init_base_P$y_inf[which(init_base_P$s=="F")]=round(rnorm(length(init_base_P$y_inf[which(init_base_P$s=="F")]),y_inf_females,2.5),digits=2) # random deviate of normal distribution with mean y_inf_females and standard deviation of 2.5 cm
  init_base_P$y_inf[which(init_base_P$s=="M")]=round(rnorm(length(init_base_P$y_inf[which(init_base_P$s=="M")]),y_inf_males,2.5),digits=2) # random deviate of normal distribution with mean y_inf_males and standard deviation of 2.5 cm
  
  ## Initialization of breeding values
  y_inf_mean_f=mean(init_base_P$y_inf[which(init_base_P$s=="F")])
  y_inf_mean_m=mean(init_base_P$y_inf[which(init_base_P$s=="M")])
  
  Nf=nrow(subset(init_base_P,init_base_P$s=="F"))
  Nm=nrow(subset(init_base_P,init_base_P$s=="M"))
  init_base_P$bv_f[which(init_base_P$s=="F")]=round((init_base_P$y_inf[which(init_base_P$s=="F")]-y_inf_mean_f)+rnorm(rep(0,Nf),0,2),digits=2)
  init_base_P$bv_m[which(init_base_P$s=="F")]=round((init_base_P$y_inf[which(init_base_P$s=="F")]-y_inf_mean_f)+rnorm(rep(0,Nf),0,10),digits=2)
  init_base_P$bv_f[which(init_base_P$s=="M")]=round((init_base_P$y_inf[which(init_base_P$s=="M")]-y_inf_mean_m)+rnorm(rep(0,Nm),0,10),digits=2)
  init_base_P$bv_m[which(init_base_P$s=="M")]=round((init_base_P$y_inf[which(init_base_P$s=="M")]-y_inf_mean_m)+rnorm(rep(0,Nm),0,2),digits=2)
  
  ##### Figures
  PLOT<-FALSE
  if (PLOT){
    ## Figures of bv for females
    # -> positive relation betwwen female bv and female asymptotic horn length
    # -> no relation betwwen female bv and male bv, nor between male bv and female asymptotic horn length
    plot(init_base_P$bv_f[which(init_base_P$s=="F")],init_base_P$bv_m[which(init_base_P$s=="F")])
    plot(init_base_P$bv_f[which(init_base_P$s=="F")],init_base_P$bv_m[which(init_base_P$s=="F")])
    plot(init_base_P$bv_f[which(init_base_P$s=="F")],init_base_P$y_inf[which(init_base_P$s=="F")])
    plot(init_base_P$bv_m[which(init_base_P$s=="F")],init_base_P$y_inf[which(init_base_P$s=="F")])
    
    ## Figures of bv for males
    # -> positive relation betwwen male bv and male asymptotic horn length
    # -> no relation betwwen male bv and female bv, nor between female bv and male asymptotic horn length
    plot(init_base_P$bv_f[which(init_base_P$s=="M")],init_base_P$bv_m[which(init_base_P$s=="M")])
    plot(init_base_P$bv_m[which(init_base_P$s=="M")],init_base_P$y_inf[which(init_base_P$s=="M")])
    plot(init_base_P$bv_f[which(init_base_P$s=="M")],init_base_P$y_inf[which(init_base_P$s=="M")])
  } 
  
  ########## Metapopulation
  ##### Creation of initial data frame
  init_base=rbind(init_base_H,init_base_P)
  init_base$id=seq(1,nrow(init_base))
  
  ## Verifications of initial distributions
  # table(init_base$a)
  # table(init_base$o)
  # table(init_base$s)
  # table(init_base$m)

  ##### Attribute asymptotic horn length and horn length to each individuals
  for (i in seq(1,nrow(init_base))){
    ## Step 1. Get asymptotic horn length from breeding values
    init_base$y_inf[i]=Asympt.Length(init_base$s[i],init_base$bv_f[i],init_base$bv_m[i],y_inf_mean_f,y_inf_mean_m,E_ff,E_mm)
    ## Step 2. Get horn length from age and asymptotic horn length
    init_base$y[i]=AsymptHornLgth.And.Age.To.HornLgth(init_base$a[i],init_base$s[i],init_base$y_inf[i])
  }
  
  ## grouping for visualisation
  init_base$y_inf_group<-as.factor(round(init_base$y_inf,digits=0))
  init_base$bv_f_group<-as.factor(round(init_base$bv_f,digits=0))
  init_base$bv_m_group<-as.factor(round(init_base$bv_m,digits=0))
  
  ##### Figures
  PLOT<-FALSE
  if (PLOT){
    
    subset(init_base,(init_base$s=="M"))%>%group_by(bv_m_group,s)%>%
      ggplot(aes(x=bv_m_group,y=y_inf,col=s))+
      labs(title="",x =TeX('$\\bv_{m}$'), y = TeX('$\\y_{\\infty,m}$'),color="Sex") +
      scale_color_brewer(palette = "Set2") +
      geom_point() +
      facet_wrap(~s) +
      theme_tufte() +
      theme(axis.line=element_line(color='black'),text=element_text(size=22))
  
    subset(init_base,(init_base$s=="M"))%>%group_by(bv_m_group,s)%>%
      ggplot(aes(x=bv_f_group,y=y_inf,col=s))+
      labs(title="",x =TeX('$\\bv_{f}$'), y = TeX('$\\y_{\\infty,m}$'),color="Sex") +
      scale_color_brewer(palette = "Set2") +
      geom_point() +
      facet_wrap(~s) +
      theme_tufte() +
      theme(axis.line=element_line(color='black'),text=element_text(size=22))
    
    ggplot(init_base,aes(x=a,y=y,col=s))+
      labs(title="",x =TeX('a'), y = TeX('y'),color="Sex") +
      scale_color_brewer(palette = "Set2") +
      geom_point() +
      facet_wrap(~s) +
      theme_tufte() +
      theme(axis.line=element_line(color='black'),text=element_text(size=22))
  }
  
  ##### Final base
  init_base_final=init_base
  init_base_final=subset(init_base_final,select=-c(y_inf_group,bv_f_group,bv_m_group))
 
  #######################################################################################################################################
  ########################################         Initial populations savings                       #################################################    
  #######################################################################################################################################
  init.pop.id=p
  init.pop_directory=file.path(getwd(),"/Inits/")
  
  ##### Saving of initial population in external files
  init.pop_file=paste("/Init_pop_",init.pop.id,".xlsx",sep="")
  write.xlsx(init_base_final,paste(init.pop_directory,init.pop_file,sep="/"),colNames=TRUE,rowNames=FALSE,startRow=1)

} # end loop of all initial meta-populations to generate 

