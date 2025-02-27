#######################################################################################################################################
########################################                                                ################################################
########################################                Simulation Model                #################################################    
########################################            Simulation Design and Runs          ################################################
########################################                                                ################################################
#######################################################################################################################################

rm(list=ls()) 
## !!! NB. Below put your folder address instead of "[...]" !!!
setwd("[...]/SimulationModel/")

#######################################################################################################################################
########################################         Load of useful resources               #################################################    
#######################################################################################################################################

####### Librairies ########
if (!require("openxlsx")) install.packages("openxlsx") # file extraction
if (!require("tidyr")) install.packages("tidyr") # create tidy data
if (!require("splines")) install.packages("splines")  #handling of model summaries
if (!require("mvtnorm")) install.packages("mvtnorm")  # multivariate normal distributions

######## Simulation Functions ########
source("./Scripts/Run_Model/SimulationFunctions.R")

#######################################################################################################################################
########################################         Selection of simulation design         #################################################    
#######################################################################################################################################

######## Definition of possible scenarios ########
#### Scenarios to simulate different:
#### - migration rates modulated through parameters epsilon_PH/HP
#### - harvest intensities modulated through parameters H_max, y_min, phi
## For final version

# -> To run
scen.eps=c(0.00,0.25,0.50,0.75)  #  flow scale for migration rates
scen.quota=c("NO_QUOTA") # c("NO_QUOTA","QUOTA") # quotas
scen.min.length=c("4_5Curl") # c("4_5Curl","Full_Curl") # min. legal horn size
scen.phi=c(0.00,0.35,0.70) # harvest rate of legal rams
scen.season=c(0.00,0.35,0.70) # scale hunt risk of migrant through harvest timing (30/09, 15/10, 31/10, 15/11)

# All combinations of scenario parameters
all_scen=crossing(scen.eps,scen.quota,scen.min.length,scen.phi,scen.season) 
all_scen_bis=subset(all_scen,!((all_scen$scen.eps==0)&(all_scen$scen.season!=0)|(all_scen$scen.phi==0)&(all_scen$scen.season!=0)))
all_scen=all_scen_bis

### Verification of possible scenarios
# all_scen
# length(scen.eps)
# length(scen.quota)
# length(scen.min.length)
# scen.eps
# scen.quota
# scen.min.length
# any((all_scen$scen.eps==0)&(all_scen$scen.season!=0))
# any((all_scen$scen.phi==0)&(all_scen$scen.season!=0))

######## Run over all selected Scenarios ########
for(Scen in (1:nrow(all_scen))){
  ## Update current scenario characteristics
  scen.eps=all_scen$scen.eps[Scen]
  scen.quota=all_scen$scen.quota[Scen]
  scen.min.length=all_scen$scen.min.length[Scen]
  scen.phi=all_scen$scen.phi[Scen]
  scen.season=all_scen$scen.season[Scen]
  scen.label=paste(as.character(scen.eps),scen.quota,scen.min.length,as.character(scen.phi),scen.season,sep="-")
  
  ## Print scenario characteristics
  cat("\n")
  print(paste("SCENARIO Nº: ",Scen,"/",nrow(all_scen),sep=""))
  print(paste("Migration rate:",as.character(scen.eps),"--- Quota:",scen.quota,"--- Legal size:",scen.min.length,"--- Harvest rate:",as.character(scen.phi),"--- Harvest season:",scen.season))

#######################################################################################################################################
########################################         Initialisations                        #################################################    
#######################################################################################################################################

  ######## Model variables ########
  #### Number of replicated runs of the same scenario
  # n:: replicate
  Nrep=10 # 10
  
  #### Simulation time
  # t:: time step 1t=1year
  TI=25 # 25 # TI: Initial simulation time to train and stabilize pop dynamics
  TF=TI+75 # TI+75 # TF: Final simulation time
  TH=TI+15 # TI+15 # TH: Start of harvest selection in time step
  
  ######## Model parameters ########
  #### 1) Aging and growing
  rho_females=0.60 # [0.57;0.63] # rho_females: female growth rate
  rho_males=0.22 # [0.21;0.23] # rho_males: male growth rate
  
  #### 2) Natural mortality 
  max_a_females=20 # max_a_females: maximal female longevity
  max_a_males=15 # max_a_males: maximal male longevity
  mu_females=c(2.059,-0.172,-0.001) # +/- c(0.143,0.028,0.006) # mu_females: female death rate coefficients of linear combination with age
  mu_males=c(1.028,-0.165,0.0006) # +/- c(0.183,0.096,0.017) # mu_males: male death rate coefficients of linear combination with age

  #### 3) Harvesting !!! -> selected with scenario !!! 
  phi=scen.phi # harvest rate of legal rams
  H_season=scen.season # harvest timing relative to breeding migrations
  # define H_max when quota=TRUE otherwise H_max/y_min=+infinity
  H_max=ifelse(scen.quota=="QUOTA", 10, Inf) # maximum legal harvest !!! -> selected with scenario !!! 
  # define y_min when min.length="4_5Curl"/"Full_Curl" otherwise y_min=0
  # y_min=ifelse(scen.min.length=="4_5Curl",76.4,0)+ifelse(scen.min.length=="Full_Curl",85.0,0) # minimum legal trophy size !!! -> selected with scenario !!! 
  y_min=ifelse(scen.min.length=="4_5Curl",76.8,0)+ifelse(scen.min.length=="Full_Curl",85.0,0) # minimum legal trophy size !!! -> selected with scenario !!! 
  
  #### 4) Migration !!! -> selected with scenario !!! 
  eps=scen.eps # factor applied to scale gene flow between both areas !!! -> selected with scenario !!! 
  lambda=c(-10.79+0.03*72,9.02,-0.61,-5.31,0.49) # c([-30.74;8.86],[4.00;15.35],[-4.58;3.19],[-11.37;0.39],[0.10;0.92]) # influence of age and nb. of competitors (intercept with correction for effect of sampling effort)

  #### 5) Reproduction 
  theta=c(-2.96,0.65,-0.04,-0.01) # c([-3.10;-1.91],[0.51;0.80],[-0.05;-0.03],[-0.01;0.00]) # theta: yearling production rate, as coefficients of linear combination with ewe's age and pop. size
  gamma=c(-6.21,0.08) # c([-7.60;-4.96],[0.06;0.09]) # gamma: coefficients of linear combination with rams's horn length influence for paternity assignation
  
  #### 6) Offspring 
  sigma=0.5 # sigma: sex ratio
  hf=0.223 # [0.090;0.446] # hf: female horn length heritability
  hm=0.397 # [0.203;0.534] # hm: male horn length heritability
  G_mm=17.884 # [9.82–25.881] # G_mm: intra-sexual genetic variance for male horn length
  G_mf=5.345 # [1.928–8.144] # G_mf: inter-sexual genetic variance horn length
  G_ff=1.622 # [0.748–3.963] # G_ff: intra-sexual genetic variance for female horn length
  E_ff=round(G_ff/hf-G_ff,digits=2) # E_ff: intra-sexual environmental variance for female horn length (VP=VA/h2;VE=VP-VA=VA/h2-VA)
  E_mm=round(G_mm/hm-G_mm,digits=2) # E_mm: intra-sexual environmental variance for male horn length (VP=VA/h2;VE=VP-VA=VA/h2-VA)
  
  ######## Initial population ########
  #### Random selection an initial population in external files
  init.pop.id=sample(1:1,1)
  init.pop_directory=file.path(getwd(),"/Inits/")
  init.pop_file=paste("/Init_pop_",init.pop.id,".xlsx",sep="")
  init_base=read.xlsx(paste(init.pop_directory,init.pop_file,sep=""),colNames=TRUE,rowNames=FALSE)
  
  ### Compute scalers for density effects on lambing and harvest
  No=nrow(init_base)

  #### Individuals variables
  # id:: unique identification nb.
  # a:: age (in time step)
  # s:: sex ("F": female, "M": male)
  # o:: origin ("P": protected, "H": harvested)
  # m:: migration status (TRUE/FALSE: if migrate during current year)
  # r:: total record of yearling produced
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
  
  #######################################################################################################################################
  ########################################         Model runs                             #################################################    
  #######################################################################################################################################
  
  ######## Runs ########
  ##### Loop for replication of scenarios
  for(n in 1:Nrep){
    stop=FALSE # flag stopping a replicate of a scenario in case of simulation error
    
    ##### Print simulation progress
    cat("\n")
    print(paste("--- SCENARIO Nº: ",Scen,"/",nrow(all_scen),", REPLICATE Nº:", n, "/", Nrep," ---", sep=""))
    
    ##### Creation of output variables
    ## -> Track recruitment
    recruitment_base=data.frame()
    ## -> Track migration
    migration_base=cbind(data.frame(t=integer(),NmP=integer(),NmH=integer(),init_base[0,]))
    ## -> Track harvest
    harvest_base=cbind(data.frame(t=integer(),init_base[0,]))
    ## -> Track average population characteristics over simulation time for all replicates of all scenarios
    out.t=rep(NA,TF+1) # time step
    ### Mean output values
    out.a=rep(NA,TF+1) # mean age at each time step
    out.aH=rep(NA,TF+1) # mean age in harvested population at each time step
    out.aP=rep(NA,TF+1) # mean age in protected population  at each time step
    out.am=rep(NA,TF+1) # mean male age at each time step
    out.amH=rep(NA,TF+1) # mean male age in harvested population at each time step
    out.amP=rep(NA,TF+1) # mean male age in protected population  at each time step
    out.yf=rep(NA,TF+1) # mean male horn size at each time step
    out.yfH=rep(NA,TF+1) # mean male horn size in harvested population at each time step
    out.yfP=rep(NA,TF+1) # mean male  horn size in protected population  at each time step
    out.ym=rep(NA,TF+1) # mean female horn size at each time step
    out.ymH=rep(NA,TF+1) # mean female horn size in harvested population at each time step
    out.ymP=rep(NA,TF+1) # mean female horn size in protected population  at each time step
    out.ym_six=rep(NA,TF+1) # mean 6y-male horn size at each time step
    out.ymH_six=rep(NA,TF+1) # mean 6y-male horn size in harvested population  at each time step
    out.ymP_six=rep(NA,TF+1) # mean 6y-male horn size in protected population  at each time step
    out.ym_sev=rep(NA,TF+1) # mean 7y-male horn size at each time step
    out.ymH_sev=rep(NA,TF+1) # mean 7y-male horn size in harvested population  at each time step
    out.ymP_sev=rep(NA,TF+1) # mean 7y-male horn size in protected population  at each time step
    out.ym_sevandolder=rep(NA,TF+1) # mean 7y-and-older-male horn size at each time step
    out.ymH_sevandolder=rep(NA,TF+1) # mean 7y-and-older-male  horn size in harvested population  at each time step
    out.ymP_sevandolder=rep(NA,TF+1)  # mean 7y-and-older-male  horn size in protected population  at each time step
    out.yf_inf=rep(NA,TF+1) # mean female asymptotic horn size at each time step
    out.yfH_inf=rep(NA,TF+1) # mean female asymptotic horn size in harvested population at each time step
    out.yfP_inf=rep(NA,TF+1) # mean female asymptotic horn size in protected population at each time step
    out.ym_inf=rep(NA,TF+1) # mean male asymptotic horn size at each time step
    out.ymH_inf=rep(NA,TF+1) # mean male asymptotic horn size in harvested population at each time step
    out.ymP_inf=rep(NA,TF+1) # mean male asymptotic horn size in protected population at each time step
    out.bvf_f=rep(NA,TF+1) # mean female breeding value for female horn length
    out.bvfH_f=rep(NA,TF+1)  # mean female breeding value for female horn length
    out.bvfP_f=rep(NA,TF+1) #mean female breeding value for female horn length
    out.bvf_m=rep(NA,TF+1)  # mean female breeding value for male horn length
    out.bvfH_m=rep(NA,TF+1) # mean female breeding value for male horn length
    out.bvfP_m=rep(NA,TF+1) #mean female breeding value for male horn length
    out.bvm_f=rep(NA,TF+1) # mean female breeding value for female horn length
    out.bvmH_f=rep(NA,TF+1) # mean female breeding value for female horn length
    out.bvmP_f=rep(NA,TF+1) #mean female breeding value for female horn length
    out.bvm_m=rep(NA,TF+1) # mean female breeding value for male horn length
    out.bvmH_m=rep(NA,TF+1) # mean female breeding value for male horn length
    out.bvmP_m=rep(NA,TF+1) # mean female breeding value for male horn length
    
    out.bvf_f_yearl=rep(NA,TF+1) # mean female breeding value for female horn length
    out.bvfH_f_yearl=rep(NA,TF+1)  # mean female breeding value for female horn length
    out.bvfP_f_yearl=rep(NA,TF+1) #mean female breeding value for female horn length
    out.bvf_m_yearl=rep(NA,TF+1)  # mean female breeding value for male horn length
    out.bvfH_m_yearl=rep(NA,TF+1) # mean female breeding value for male horn length
    out.bvfP_m_yearl=rep(NA,TF+1) #mean female breeding value for male horn length
    out.bvm_f_yearl=rep(NA,TF+1) # mean female breeding value for female horn length
    out.bvmH_f_yearl=rep(NA,TF+1) # mean female breeding value for female horn length
    out.bvmP_f_yearl=rep(NA,TF+1) #mean female breeding value for female horn length
    out.bvm_m_yearl=rep(NA,TF+1) # mean female breeding value for male horn length
    out.bvmH_m_yearl=rep(NA,TF+1) # mean female breeding value for male horn length
    out.bvmP_m_yearl=rep(NA,TF+1) # mean female breeding value for male horn length
    
    ### Variance output values
    out.sd.a=rep(NA,TF+1) # sd age at each time step
    out.sd.aH=rep(NA,TF+1) # sd age in harvested population at each time step
    out.sd.aP=rep(NA,TF+1) # sd age in protected population  at each time step
    out.sd.yf=rep(NA,TF+1) # sd female horn size at each time step
    out.sd.yfH=rep(NA,TF+1) # sd female horn size in harvested population at each time step
    out.sd.yfP=rep(NA,TF+1) # sd female  horn size in protected population  at each time step
    out.sd.ym=rep(NA,TF+1) # sd female horn size at each time step
    out.sd.ymH=rep(NA,TF+1) # sd female horn size in harvested population at each time step
    out.sd.ymP=rep(NA,TF+1) # sd female horn size in protected population  at each time step
    out.sd.yf_inf=rep(NA,TF+1) # sd female asymptotic horn size at each time step
    out.sd.yfH_inf=rep(NA,TF+1) # sd female asymptotic horn size in harvested population at each time step
    out.sd.yfP_inf=rep(NA,TF+1) # sd female asymptotic horn size in protected population at each time step
    out.sd.ym_inf=rep(NA,TF+1) # sd male asymptotic horn size at each time step
    out.sd.ymH_inf=rep(NA,TF+1) # sd male asymptotic horn size in harvested population at each time step
    out.sd.ymP_inf=rep(NA,TF+1) # sd male asymptotic horn size in protected population at each time step
    out.sd.bvf_f=rep(NA,TF+1) # sd female breeding value for female horn length
    out.sd.bvfH_f=rep(NA,TF+1)  # sd female breeding value for female horn length in harvested population at each time step
    out.sd.bvfP_f=rep(NA,TF+1) #sd female breeding value for female horn length in protected population at each time step
    out.sd.bvf_m=rep(NA,TF+1)  # sd female breeding value for male horn length at each time step
    out.sd.bvfH_m=rep(NA,TF+1) # sd female breeding value for male horn length in harvested population at each time step
    out.sd.bvfP_m=rep(NA,TF+1) #sd female breeding value for male horn length in protected population at each time step at each time step
    out.sd.bvm_f=rep(NA,TF+1) # sd female breeding value for female horn length at each time step
    out.sd.bvmH_f=rep(NA,TF+1) # sd female breeding value for female horn length in harvested population at each time step
    out.sd.bvmP_f=rep(NA,TF+1) #sd female breeding value for female horn length in protected population at each time step at each time step
    out.sd.bvm_m=rep(NA,TF+1) # sd female breeding value for male horn length at each time step
    out.sd.bvmH_m=rep(NA,TF+1) # sd female breeding value for male horn length in harvested population at each time step
    out.sd.bvmP_m=rep(NA,TF+1) # sd female breeding value for male horn length in protected population at each time step at each time step
    out.sd.bvf_f_yearl=rep(NA,TF+1) # sd female breeding value for female horn length
    out.sd.bvfH_f_yearl=rep(NA,TF+1)  # sd female breeding value for female horn length in harvested population at each time step
    out.sd.bvfP_f_yearl=rep(NA,TF+1) #sd female breeding value for female horn length in protected population at each time step
    out.sd.bvf_m_yearl=rep(NA,TF+1)  # sd female breeding value for male horn length at each time step
    out.sd.bvfH_m_yearl=rep(NA,TF+1) # sd female breeding value for male horn length in harvested population at each time step
    out.sd.bvfP_m_yearl=rep(NA,TF+1) # sd female breeding value for male horn length in protected population at each time step at each time step
    out.sd.bvm_f_yearl=rep(NA,TF+1) # sd female breeding value for female horn length at each time step
    out.sd.bvmH_f_yearl=rep(NA,TF+1) # sd female breeding value for female horn length in harvested population at each time step
    out.sd.bvmP_f_yearl=rep(NA,TF+1) # sd female breeding value for female horn length in protected population at each time step at each time step
    out.sd.bvm_m_yearl=rep(NA,TF+1) # sd female breeding value for male horn length at each time step
    out.sd.bvmH_m_yearl=rep(NA,TF+1) # sd female breeding value for male horn length in harvested population at each time step
    out.sd.bvmP_m_yearl=rep(NA,TF+1) # sd female breeding value for male horn length in protected population at each time step at each time step
    out.sd.ym_six=rep(NA,TF+1) # sd 6y-male horn size at each time step
    out.sd.ymH_six=rep(NA,TF+1) # sd 6y-male horn size in harvested population  at each time step
    out.sd.ymP_six=rep(NA,TF+1) # sd 6y-male horn size in protected population  at each time step
    out.sd.ym_sev=rep(NA,TF+1) # sd 7y-male horn size at each time step
    out.sd.ymH_sev=rep(NA,TF+1) # sd 7y-male horn size in harvested population  at each time step
    out.sd.ymP_sev=rep(NA,TF+1) # sd 7y-male horn size in protected population  at each time step
    
    ### Counting output values
    out.H=rep(NA,TF+1) # harvest at each time step
    out.prop_m_45=rep(NA,TF+1) # proportion of aged 4-5 year old males at each time step
    out.prop_mH_45=rep(NA,TF+1) # proportion of aged 4-5 year old males in harvested population at each time step
    out.prop_mP_45=rep(NA,TF+1) # proportion of aged 4-5 year old males in protected population at each time step
    out.struct_a=rep(NA,TF+1) # age structure at each time step
    out.struct_aH=rep(NA,TF+1) # age structure in harvested population at each time step
    out.struct_aP=rep(NA,TF+1) # age structure in protected population at each time step
    out.q=rep(NA,TF+1) # proportion of legal rams at each time step
    out.qH=rep(NA,TF+1) # proportion of legal rams in harvested population at each time step
    out.qP=rep(NA,TF+1) # proportion of legal rams in protected population at each time step
    out.sr=rep(NA,TF+1) # sex-ratio at each time step
    out.srH=rep(NA,TF+1) # sex-ratio in harvested population at each time step
    out.srP=rep(NA,TF+1) # sex-ratio in protected population at each time step
    out.N=rep(NA,TF+1) # metapopulation size at each time step
    out.NH=rep(NA,TF+1) # harvested population size at each time step
    out.NP=rep(NA,TF+1) # protected population size at each time step
    out.Nf=rep(NA,TF+1) # number of females at each time step
    out.NfH=rep(NA,TF+1) # number of females at each time step
    out.NfP=rep(NA,TF+1) # number of females at each time step
    out.Nm=rep(NA,TF+1) # number of males at each time step
    out.NmH=rep(NA,TF+1) # number of males in harvested population at each time step
    out.NmP=rep(NA,TF+1) # number of males in protected population  at each time step
    out.M_PtoH=rep(NA,TF+1) # number of migrants from P to H
    out.M_HtoP=rep(NA,TF+1) # number of migrants from H to P
    out.M_PtoH_Risk=rep(NA,TF+1) # number of migrants from P to H before harvest
    out.M_HtoP_NoRisk=rep(NA,TF+1) # number of migrants from H to P before harvest
    out.r_PtoH=rep(NA,TF+1) # migration rate from P to H
    out.r_HtoP=rep(NA,TF+1) # migration rate from H to P
    out.r_PtoH_Risk=rep(NA,TF+1) # migration rate from P to H before harvest
    out.r_HtoP_NoRisk=rep(NA,TF+1) # migration rate from H to P before harvest
    
    ##### Initialization of output variables
    out.t[1]=0 # time step
    data_base_out=subset(init_base,(init_base$a>0))
    ### Average output values
    out.a[1]=mean(data_base_out$a) 
    out.aH[1]=mean(data_base_out$a[which(data_base_out$o=="H")]) 
    out.aP[1]=mean(data_base_out$a[which(data_base_out$o=="P")]) 
    out.am[1]=mean(data_base_out$a[which(data_base_out$s=="M")]) 
    out.amH[1]=mean(data_base_out$a[which((data_base_out$o=="H")&(data_base_out$s=="M"))]) 
    out.amP[1]=mean(data_base_out$a[which((data_base_out$o=="P")&(data_base_out$s=="M"))]) 
    out.yf[1]=mean(data_base_out$y[which(data_base_out$s=="F")]) 
    out.yfH[1]=mean(data_base_out$y[which((data_base_out$o=="H")&(data_base_out$s=="F"))]) 
    out.yfP[1]=mean(data_base_out$y[which((data_base_out$o=="P")&(data_base_out$s=="F"))]) 
    out.ym[1]=mean(data_base_out$y[which(data_base_out$s=="M")]) 
    out.ymH[1]=mean(data_base_out$y[which((data_base_out$o=="H")&(data_base_out$s=="M"))]) 
    out.ymP[1]=mean(data_base_out$y[which((data_base_out$o=="P")&(data_base_out$s=="M"))]) 
    out.ym_six[1]=mean(data_base_out$y[which((data_base_out$a==6)&(data_base_out$s=="M"))])
    out.ymH_six[1]=mean(data_base_out$y[which((data_base_out$a==6)&(data_base_out$o=="H")&(data_base_out$s=="M"))]) 
    out.ymP_six[1]=mean(data_base_out$y[which((data_base_out$a==6)&(data_base_out$o=="P")&(data_base_out$s=="M"))])    
    out.ym_sev[1]=mean(data_base_out$y[which((data_base_out$a==7)&(data_base_out$s=="M"))])
    out.ymH_sev[1]=mean(data_base_out$y[which((data_base_out$a==7)&(data_base_out$o=="H")&(data_base_out$s=="M"))]) 
    out.ymP_sev[1]=mean(data_base_out$y[which((data_base_out$a==7)&(data_base_out$o=="P")&(data_base_out$s=="M"))])
    out.ym_sevandolder[1]=mean(data_base_out$y[which((data_base_out$a>=7)&(data_base_out$s=="M"))]) 
    out.ymH_sevandolder[1]=mean(data_base_out$y[which((data_base_out$a>=7)&(data_base_out$o=="H")&(data_base_out$s=="M"))])
    out.ymP_sevandolder[1]=mean(data_base_out$y[which((data_base_out$a>=7)&(data_base_out$o=="P")&(data_base_out$s=="M"))])
    out.yf_inf[1]=mean(data_base_out$y_inf[which(data_base_out$s=="F")]) 
    out.yfH_inf[1]=mean(data_base_out$y_inf[which((data_base_out$o=="H")&(data_base_out$s=="F"))]) 
    out.yfP_inf[1]=mean(data_base_out$y_inf[which((data_base_out$o=="P")&(data_base_out$s=="F"))]) 
    out.ym_inf[1]=mean(data_base_out$y_inf[which(data_base_out$s=="M")]) 
    out.ymH_inf[1]=mean(data_base_out$y_inf[which((data_base_out$o=="H")&(data_base_out$s=="M"))]) 
    out.ymP_inf[1]=mean(data_base_out$y_inf[which((data_base_out$o=="P")&(data_base_out$s=="M"))]) 
    out.bvf_f[1]=mean(data_base_out$bv_f[which((data_base_out$s=="F"))]) 
    out.bvfH_f[1]=mean(data_base_out$bv_f[which((data_base_out$o=="H")&(data_base_out$s=="F"))]) 
    out.bvfP_f[1]=mean(data_base_out$bv_f[which((data_base_out$o=="P")&(data_base_out$s=="F"))]) 
    out.bvf_m[1]=mean(data_base_out$bv_m[which((data_base_out$s=="F"))]) 
    out.bvfH_m[1]=mean(data_base_out$bv_m[which((data_base_out$o=="H")&(data_base_out$s=="F"))]) 
    out.bvfP_m[1]=mean(data_base_out$bv_m[which((data_base_out$o=="P")&(data_base_out$s=="F"))]) 
    out.bvm_f[1]=mean(data_base_out$bv_f[which((data_base_out$s=="M"))]) 
    out.bvmH_f[1]=mean(data_base_out$bv_f[which((data_base_out$o=="H")&(data_base_out$s=="M"))])
    out.bvmP_f[1]=mean(data_base_out$bv_f[which((data_base_out$o=="P")&(data_base_out$s=="M"))]) 
    out.bvm_m[1]=mean(data_base_out$bv_m[which((data_base_out$s=="M"))]) 
    out.bvmH_m[1]=mean(data_base_out$bv_m[which((data_base_out$o=="H")&(data_base_out$s=="M"))]) 
    out.bvmP_m[1]=mean(data_base_out$bv_m[which((data_base_out$o=="P")&(data_base_out$s=="M"))]) 
    out.bvf_f_yearl[1]=mean(data_base_out$bv_f[which((data_base_out$s=="F")&(data_base_out$a==1))]) 
    out.bvfH_f_yearl[1]=mean(data_base_out$bv_f[which((data_base_out$o=="H")&(data_base_out$s=="F")&(data_base_out$a==1))]) 
    out.bvfP_f_yearl[1]=mean(data_base_out$bv_f[which((data_base_out$o=="P")&(data_base_out$s=="F")&(data_base_out$a==1))]) 
    out.bvf_m_yearl[1]=mean(data_base_out$bv_m[which((data_base_out$s=="F")&(data_base_out$a==1))]) 
    out.bvfH_m_yearl[1]=mean(data_base_out$bv_m[which((data_base_out$o=="H")&(data_base_out$s=="F"))]) 
    out.bvfP_m_yearl[1]=mean(data_base_out$bv_m[which((data_base_out$o=="P")&(data_base_out$s=="F")&(data_base_out$a==1))]) 
    out.bvm_f_yearl[1]=mean(data_base_out$bv_f[which((data_base_out$s=="M")&(data_base_out$a==1))]) 
    out.bvmH_f_yearl[1]=mean(data_base_out$bv_f[which((data_base_out$o=="H")&(data_base_out$s=="M")&(data_base_out$a==1))])
    out.bvmP_f_yearl[1]=mean(data_base_out$bv_f[which((data_base_out$o=="P")&(data_base_out$s=="M")&(data_base_out$a==1))]) 
    out.bvm_m_yearl[1]=mean(data_base_out$bv_m[which((data_base_out$s=="M")&(data_base_out$a==1))]) 
    out.bvmH_m_yearl[1]=mean(data_base_out$bv_m[which((data_base_out$o=="H")&(data_base_out$s=="M")&(data_base_out$a==1))]) 
    out.bvmP_m_yearl[1]=mean(data_base_out$bv_m[which((data_base_out$o=="P")&(data_base_out$s=="M")&(data_base_out$a==1))]) 
    
    ### Variance output values
    out.sd.a[1]=sd(data_base_out$a) 
    out.sd.aH[1]=sd(data_base_out$a[which(data_base_out$o=="H")]) 
    out.sd.aP[1]=sd(data_base_out$a[which(data_base_out$o=="P")]) 
    out.sd.yf[1]=sd(data_base_out$y[which(data_base_out$s=="F")]) 
    out.sd.yfH[1]=sd(data_base_out$y[which((data_base_out$o=="H")&(data_base_out$s=="F"))]) 
    out.sd.yfP[1]=sd(data_base_out$y[which((data_base_out$o=="P")&(data_base_out$s=="F"))]) 
    out.sd.ym[1]=sd(data_base_out$y[which(data_base_out$s=="M")])
    out.sd.ymH[1]=sd(data_base_out$y[which((data_base_out$o=="H")&(data_base_out$s=="M"))]) 
    out.sd.ymP[1]=sd(data_base_out$y[which((data_base_out$o=="P")&(data_base_out$s=="M"))]) 
    out.sd.yf_inf[1]=sd(data_base_out$y_inf[which(data_base_out$s=="F")]) 
    out.sd.yfH_inf[1]=sd(data_base_out$y_inf[which((data_base_out$o=="H")&(data_base_out$s=="F"))]) 
    out.sd.yfP_inf[1]=sd(data_base_out$y_inf[which((data_base_out$o=="P")&(data_base_out$s=="F"))]) 
    out.sd.ym_inf[1]=sd(data_base_out$y_inf[which(data_base_out$s=="M")]) 
    out.sd.ymH_inf[1]=sd(data_base_out$y_inf[which((data_base_out$o=="H")&(data_base_out$s=="M"))]) 
    out.sd.ymP_inf[1]=sd(data_base_out$y_inf[which((data_base_out$o=="P")&(data_base_out$s=="M"))]) 
    out.sd.bvf_f[1]=sd(data_base_out$bv_f[which(data_base_out$s=="F")]) 
    out.sd.bvfH_f[1]=sd(data_base_out$bv_f[which((data_base_out$o=="H")&(data_base_out$s=="F"))])
    out.sd.bvfP_f[1]=sd(data_base_out$bv_f[which((data_base_out$o=="P")&(data_base_out$s=="F"))]) 
    out.sd.bvf_m[1]=sd(data_base_out$bv_m[which(data_base_out$s=="F")]) 
    out.sd.bvfH_m[1]=sd(data_base_out$bv_m[which((data_base_out$o=="H")&(data_base_out$s=="F"))]) 
    out.sd.bvfP_m[1]=sd(data_base_out$bv_m[which((data_base_out$o=="P")&(data_base_out$s=="F"))]) 
    out.sd.bvm_f[1]=sd(data_base_out$bv_f[which(data_base_out$s=="M")]) 
    out.sd.bvmH_f[1]=sd(data_base_out$bv_f[which((data_base_out$o=="H")&(data_base_out$s=="M"))]) 
    out.sd.bvmP_f[1]=sd(data_base_out$bv_f[which((data_base_out$o=="P")&(data_base_out$s=="M"))]) 
    out.sd.bvm_m[1]=sd(data_base_out$bv_m[which(data_base_out$s=="M")]) 
    out.sd.bvmH_m[1]=sd(data_base_out$bv_m[which((data_base_out$o=="H")&(data_base_out$s=="M"))])
    out.sd.bvmP_m[1]=sd(data_base_out$bv_m[which((data_base_out$o=="P")&(data_base_out$s=="M"))]) 

    out.sd.bvf_f_yearl[1]=sd(data_base_out$bv_f[which((data_base_out$s=="F")&(data_base_out$a==1))]) 
    out.sd.bvfH_f_yearl[1]=sd(data_base_out$bv_f[which((data_base_out$o=="H")&(data_base_out$s=="F")&(data_base_out$a==1))])
    out.sd.bvfP_f_yearl[1]=sd(data_base_out$bv_f[which((data_base_out$o=="P")&(data_base_out$s=="F")&(data_base_out$a==1))]) 
    out.sd.bvf_m_yearl[1]=sd(data_base_out$bv_m[which((data_base_out$s=="F")&(data_base_out$a==1))]) 
    out.sd.bvfH_m_yearl[1]=sd(data_base_out$bv_m[which((data_base_out$o=="H")&(data_base_out$s=="F")&(data_base_out$a==1))]) 
    out.sd.bvfP_m_yearl[1]=sd(data_base_out$bv_m[which((data_base_out$o=="P")&(data_base_out$s=="F")&(data_base_out$a==1))]) 
    out.sd.bvm_f_yearl[1]=sd(data_base_out$bv_f[which((data_base_out$s=="M")&(data_base_out$a==1)&(data_base_out$a==1))]) 
    out.sd.bvmH_f_yearl[1]=sd(data_base_out$bv_f[which((data_base_out$o=="H")&(data_base_out$s=="M")&(data_base_out$a==1))]) 
    out.sd.bvmP_f_yearl[1]=sd(data_base_out$bv_f[which((data_base_out$o=="P")&(data_base_out$s=="M")&(data_base_out$a==1))]) 
    out.sd.bvm_m_yearl[1]=sd(data_base_out$bv_m[which((data_base_out$s=="M")&(data_base_out$a==1))]) 
    out.sd.bvmH_m_yearl[1]=sd(data_base_out$bv_m[which((data_base_out$o=="H")&(data_base_out$s=="M")&(data_base_out$a==1))])
    out.sd.bvmP_m_yearl[1]=sd(data_base_out$bv_m[which((data_base_out$o=="P")&(data_base_out$s=="M")&(data_base_out$a==1))]) 
    
    out.sd.ym_six[1]=sd(data_base_out$y[which((data_base_out$a==6)&(data_base_out$s=="M"))])
    out.sd.ymH_six[1]=sd(data_base_out$y[which((data_base_out$a==6)&(data_base_out$o=="H")&(data_base_out$s=="M"))]) 
    out.sd.ymP_six[1]=sd(data_base_out$y[which((data_base_out$a==6)&(data_base_out$o=="P")&(data_base_out$s=="M"))]) 
    out.sd.ym_sev[1]=sd(data_base_out$y[which((data_base_out$a==7)&(data_base_out$s=="M"))])
    out.sd.ymH_sev[1]=sd(data_base_out$y[which((data_base_out$a==7)&(data_base_out$o=="H")&(data_base_out$s=="M"))]) 
    out.sd.ymP_sev[1]=sd(data_base_out$y[which((data_base_out$a==7)&(data_base_out$o=="P")&(data_base_out$s=="M"))]) 
    
    ### Counting output values
    out.H[1]=0 
    out.prop_m_45[1]=round(nrow(subset(data_base_out,(data_base_out$a%in%c(4,5))&(data_base_out$s=="M")))/nrow(subset(data_base_out,(data_base_out$s=="M"))),digits=2) 
    out.prop_mH_45[1]=round(nrow(subset(data_base_out,(data_base_out$o=="H")&(data_base_out$a%in%c(4,5))&(data_base_out$s=="M")))/nrow(subset(data_base_out,(data_base_out$o=="H")&(data_base_out$s=="M"))),digits=2) 
    out.prop_mP_45[1]=round(nrow(subset(data_base_out,(data_base_out$o=="P")&(data_base_out$a%in%c(4,5))&(data_base_out$s=="M")))/nrow(subset(data_base_out,(data_base_out$o=="P")&(data_base_out$s=="M"))),digits=2) 
    out.struct_a[1]=round(nrow(subset(data_base_out,data_base_out$a>=7))/nrow(data_base_out),digits=2) 
    out.struct_aH[1]=round(nrow(subset(data_base_out,(data_base_out$o=="H")&(data_base_out$a>=7)))/nrow(subset(data_base_out,(data_base_out$o=="H"))),digits=2) 
    out.struct_aP[1]=round(nrow(subset(data_base_out,(data_base_out$o=="P")&(data_base_out$a>=7)))/nrow(subset(data_base_out,(data_base_out$o=="P"))),digits=2) 
    out.q[1]=nrow(subset(data_base_out,(data_base_out$s=="M")&(data_base_out$y>y_min)))/nrow(subset(data_base_out,(data_base_out$s=="M"))) 
    out.qH[1]=nrow(subset(data_base_out,(data_base_out$o=="H")&(data_base_out$s=="M")&(data_base_out$y>y_min)))/nrow(subset(data_base_out,(data_base_out$o=="H")&(data_base_out$s=="M"))) 
    out.qP[1]=nrow(subset(data_base_out,(data_base_out$o=="P")&(data_base_out$s=="M")&(data_base_out$y>y_min)))/nrow(subset(data_base_out,(data_base_out$o=="P")&(data_base_out$s=="M"))) 
    out.sr[1]=nrow(subset(data_base_out,(data_base_out$s=="F")))/nrow(subset(data_base_out,(data_base_out$s=="M"))) 
    out.srH[1]=nrow(subset(data_base_out,(data_base_out$o=="H")&(data_base_out$s=="F")))/nrow(subset(data_base_out,(data_base_out$o=="H")&(data_base_out$s=="M"))) 
    out.srP[1]=nrow(subset(data_base_out,(data_base_out$o=="P")&(data_base_out$s=="F")))/nrow(subset(data_base_out,(data_base_out$o=="P")&(data_base_out$s=="M"))) 
    out.N[1]=nrow(data_base_out) 
    out.NH[1]=nrow(subset(data_base_out,(data_base_out$o=="H"))) 
    out.NP[1]=nrow(subset(data_base_out,(data_base_out$o=="P"))) 
    out.Nm[1]=nrow(subset(data_base_out,(data_base_out$s=="M"))) 
    out.NmH[1]=nrow(subset(data_base_out,(data_base_out$o=="H")&(data_base_out$s=="M")))
    out.NmP[1]=nrow(subset(data_base_out,(data_base_out$o=="P")&(data_base_out$s=="M"))) 
    out.Nf[1]=nrow(subset(data_base_out,(data_base_out$s=="F"))) 
    out.NfH[1]=nrow(subset(data_base_out,(data_base_out$o=="H")&(data_base_out$s=="F"))) 
    out.NfP[1]=nrow(subset(data_base_out,(data_base_out$o=="P")&(data_base_out$s=="F"))) 
    out.M_PtoH[1]=0 
    out.M_HtoP[1]=0 
    out.M_PtoH_Risk[1]=0 
    out.M_HtoP_NoRisk[1]=0 
    out.r_PtoH[1]=0 
    out.r_HtoP[1]=0 
    out.r_PtoH_Risk[1]=0 
    out.r_HtoP_NoRisk[1]=0 
    
    ##### Initialisation of record variables 
    data_base=init_base # reinitialise data_base at initial base for each replicate of each scenario
    id_count=nrow(data_base) # reinitialise counter for attribution of unique identification nb.
    
    siring_male_representation=c() # reinitialise male siring competition 
    siring_male_representation_P=c() # reinitialise male siring competition 
    siring_male_representation_H=c()  # reinitialise male siring competition 
    siring_male_phenotype=c() # reinitialise siring competitors' phenotypes
    siring_male_phenotype_P=c() # reinitialise siring competitors' phenotypes 
    siring_male_phenotype_H=c() # reinitialise siring competitors' phenotypes
    dead_fathers=c() # reinitialize record of paternities from previous siring competition of naturally dead fathers before lambing
    
    ## Compute first assignation of paternity (only for t=1, then computed according to previous siring competition)
    # One siring competition in both populations
    for (origins in c("P","H")){
      origin_male_data_base=subset(data_base,(data_base$s=="M")&(data_base$o==origins)&(data_base$a>1))
      Nm=nrow(origin_male_data_base) 
      if(Nm==0){
        stop=TRUE
        print("!!! FIRST SIRING ERROR: NO MALE !!!")
        break}
      rank=rep(0,Nm) # classification of mature males in same population of reproductive ewe according to horn length
      for(j in 1:Nm){
        p_siring=Male.Repro(origin_male_data_base$y[j],gamma)
        rank[j]=p_siring
      }
      unity=rep(1,length(rank))
      siring_weights=pmax(unity,round(rank*1.8)) # weigths of male mating opportunities with horn length dependent rank
      if (origins=="P"){
        siring_male_phenotype_P=subset(origin_male_data_base,select=c('id','y','y_inf','bv_f','bv_m','m')) # save phenotype of male competitors for siring
        siring_male_representation_P=rep(origin_male_data_base$id,siring_weights) # representation of male mating opportunities with horn length dependent rank
      }else{
        siring_male_phenotype_H=subset(origin_male_data_base,select=c('id','y','y_inf','bv_f','bv_m','m')) # save phenotype of male competitors for siring
        siring_male_representation_H=rep(origin_male_data_base$id,siring_weights) # representation of male mating opportunities with horn length dependent rank
      }
    }

    
    ##### Loop over time steps
    for(t in 1:(TF)){
      
      ##### Print simulation progress
      print(ifelse(t<=TI,paste("TRAINING TIME STEP:", t, "/", TI),paste("SIMULATION TIME STEP:", t, "/", TF)))
      
      if(nrow(data_base)==0){
        stop=TRUE
        print("!!! INITIAL ERROR: EMPTY POP. !!!")
        break}
      
      ##### Reinitialisation (1 time step = 1 year)
      # -> Reinitialisation of harvest by hunting season
      H=0 

      ##### Biological and artificial processes over all individual
      
      ## (1) Aging 
      # -> Increase age of all indiv. by 1 time step   
      data_base$a=data_base$a+1
        
      ## (2) Lambing 
      # -> Compute fate of all females for producing a yearling
      female_data_base=subset(data_base,(data_base$s=="F")&(data_base$a>1))
      N_updated=nrow(data_base)
      Nf=nrow(female_data_base) 
      birth_counter=0
      if(Nf==0){
          stop=TRUE
          print("!!! LAMBING ERROR: NO FEMALE !!!")
          break}
        for(i in 1:Nf){
          # -> Compute lambing proba. according to density in pop. of origin
          # NB. scale the absolute effect of pop size on yearling production (pop size in [50,200] in RM approx., while init pop size of 11000 indiv.)
          N_size=ifelse(female_data_base$o[i]=="H",nrow(subset(data_base,(data_base$o=="H"))),nrow(subset(data_base,(data_base$o=="P"))))
          N_scale_H=0.1*(round((200-50)/(700-300)*(N_size-300)+50,digits=2)) 
          N_scale_P=0.1*(round((200-50)/(400-200)*(N_size-200)+50,digits=2))
          N_scale=ifelse(female_data_base$o[i]=="H",N_scale_H,N_scale_P) 
          p_birth=Female.Repro(female_data_base$a[i],theta,N_scale)
          birth=rbinom(1,1,p_birth) # 1 random value in [0,1] with proba. of p_birth
          # -> Compute characteristics of new yearling indiv.
          if(birth==1){
            mother.id=female_data_base$id[i]
            birth_counter=birth_counter+1
            data_base$r[which(data_base$id==mother.id)]=data_base$r[which(data_base$id==mother.id)]+1 # update nb. of offspring of mother
            # -> Assign paternity according to previous siring competition
            origin_mother=female_data_base$o[i]
            if (origin_mother=="P"){
              siring_male_representation=siring_male_representation_P
              siring_male_phenotype=siring_male_phenotype_P
            }else{
              siring_male_representation=siring_male_representation_H
              siring_male_phenotype=siring_male_phenotype_H
            }
            # -> Select one father among siring competitors
            father.id=ifelse(length(siring_male_representation)==1,siring_male_representation,sample(siring_male_representation,1)) # assign random paternity in representation of male mating opportunities
            # -> Save male siring success and compute paternal characteristics of new yearling indiv. (for following computation of new.y_inf)
            if (any(data_base$id==father.id)){
              data_base$r[which(data_base$id==father.id)]=data_base$r[which(data_base$id==father.id)]+1 # update nb. of offspring of father
              BV_father_f=data_base$bv_f[which(data_base$id==father.id)]
              BV_father_m=data_base$bv_m[which(data_base$id==father.id)]
              y_inf_father=data_base$y_inf[which(data_base$id==father.id)]
              migrant_father=data_base$m[which(data_base$id==father.id)]
              male_data_base=subset(data_base,(data_base$s=="M")&(data_base$a>0)&(data_base$o==female_data_base$o[i]))
              y_inf_mean_males=mean(male_data_base$y_inf)
           }else{
              dead_fathers=c(dead_fathers,father.id,t)
              BV_father_f=siring_male_phenotype$bv_f[which(siring_male_phenotype$id==father.id)]
              BV_father_m=siring_male_phenotype$bv_m[which(siring_male_phenotype$id==father.id)]
              y_inf_father=siring_male_phenotype$y_inf[which(siring_male_phenotype$id==father.id)]
              migrant_father=siring_male_phenotype$m[which(siring_male_phenotype$id==father.id)]
              y_inf_mean_males=mean(siring_male_phenotype$y_inf)
            }

            # -> Compute maternal characteristics of new yearling indiv. (for following computation of new.y_inf)
            BV_mother_f=data_base$bv_f[which(data_base$id==mother.id)]
            BV_mother_m=data_base$bv_m[which(data_base$id==mother.id)]
            y_inf_mother=data_base$y_inf[which(data_base$id==mother.id)]
            y_inf_mean_females=mean(female_data_base$y_inf[which(female_data_base$o==female_data_base$o[i])])

            ## (3) Offspring
            # -> Compute characteristics of all new yearling indiv.
            id_count=id_count+1 # update counter for attribution of unique identification nb.
            new.id=id_count # unique identification nb.
            new.a=0 # production of yearling appearing next time step
            new.s=Sex.Assign(sigma) # sex assignation
            new.o=female_data_base$o[i] # population of origin according to mother characteristics
            new.bv=Breeding.Value(BV_mother_f,BV_mother_m,BV_father_f,BV_father_m,G_mm,G_mf,G_ff) # breeding value for sex-specific horn length
            new.bv_f=new.bv[1] # breeding value for female horn length
            new.bv_m=new.bv[2] # breeding value for male horn length
            new.y_inf=Asympt.Length(new.s,new.bv_f,new.bv_m,y_inf_mean_females,y_inf_mean_males,E_ff,E_mm) # new born lambs asymptotic horn length computation according to parents characteristics
            new.y=0 # null lambs horn length
            # -> Add new yearling indiv. to the metapop.
            new_yearling=data.frame(id=new.id,a=new.a,s=new.s,o=new.o,m=FALSE,r=0,y=new.y,y_inf=new.y_inf,bv_f=new.bv_f,bv_m=new.bv_m)
            data_base=rbind(data_base,new_yearling)
            recruitment_base=rbind(recruitment_base,cbind(t,new.bv_f,new.bv_m,new.y_inf,new.s,origin_mother,migrant_father,BV_mother_f,BV_mother_m,BV_father_f,BV_father_m,y_inf_mean_males,y_inf_mean_females))
          }
        }

      #### Print simulation progress
        print(paste("         --> Recruitment: ",birth_counter,"/",N_updated,"=",nrow(data_base),sep=""))
        
        ## (4) Growing
        # -> Compute new horn length
        N=nrow(data_base) # metapop. size
        for(i in 1:N){
          if (data_base$a[i]>0){
            data_base$y[i]=Length.Growth(data_base$s[i],data_base$y[i],data_base$y_inf[i],rho_females,rho_males)
          }
        }
        
        ## (5) Migration 
        # -> Compute fate of all males for temporal migration between protected and exploited pop.
        Nm_P_ante_migration=nrow(subset(data_base,(data_base$s=="M")&(data_base$o=="P")&(data_base$a>0)))
        Nm_H_ante_migration=nrow(subset(data_base,(data_base$s=="M")&(data_base$o=="H")&(data_base$a>0)))
        
        if (eps!=0){ 
          male_data_base=subset(data_base,(data_base$s=="M")&(data_base$a>0))
          Nm=nrow(male_data_base)
          if(Nm==0){
            stop=TRUE
            print("!!! MIGRATION ERROR: NO MALE !!!")
            break}
          for(i in 1:Nm){
            # -> Compute migration proba. according to density in pop. of origin
            if (male_data_base$o[i]=="H"){
              Nf_available=nrow(subset(data_base,(data_base$o=="H")&(data_base$s=="F")&(data_base$a>0))) # nb. of females in pop. of origin
              Nm_comp=nrow(subset(male_data_base,(male_data_base$o=="H")&(male_data_base$a>0))) # nb. of competitors in pop. of origin
              Nm_comp_scale=round((30-10)/(325-100)*(Nm_comp-100)+10,digits=2) # scale the absolute effect of the nb. of competitors on migration proba. (nb. of comp in [10,30] in RM, while init pop size of 11000 indiv. so approx. 2500 males in each pop.)
              older_m=ifelse(Nm_comp==0,12,max(male_data_base$a[which(male_data_base$o=="H")])) # age of older male in pop. of origin
             }else{
              Nf_available=nrow(subset(data_base,(data_base$o=="P")&(data_base$s=="F")&(data_base$a>0))) # nb. of females in pop. of origin
              Nm_comp=nrow(subset(male_data_base,(male_data_base$o=="P")&(male_data_base$a>0))) # nb. of competitors in pop. of origin
              Nm_comp_scale=round((30-10)/(175-0)*(Nm_comp-0)+10,digits=2) # scale the absolute effect of the nb. of competitors on migration proba. (nb. of comp in [10,30] in RM, while init pop size of 11000 indiv. so approx. 2500 males in each pop.)
              older_m=ifelse(Nm_comp==0,12,max(male_data_base$a[which(male_data_base$o=="P")])) # age of older male in pop. of origin
             }
            sexratio=Nf_available/Nm_comp
            p_migrate=Breed.Mig(male_data_base$a[i],lambda,sexratio,Nm_comp_scale,older_m)
            migrate=rbinom(1,1,p_migrate) # 1 random value in [0,1] with proba. of p_migrate
            # -> Update migrant characteristics
            if(migrate==1){
              migrant.id=male_data_base$id[i]
              # update migration status of migrants
              data_base$m[which(data_base$id==migrant.id)]=TRUE 
            }
          }
          
          # scale effective migration rate according to eps parameter
          potential_migrants=data_base$id[which(data_base$m==TRUE)]
          ids_no_migrants=sample(as.character(potential_migrants),round((1-eps)*length(potential_migrants)),replace=FALSE)
          data_base$m[which(data_base$id%in%ids_no_migrants)]=FALSE
          
          # update pop. of origin of migrants: transition of migrant males with risk to be hunted determined later
          data_base$o[which((data_base$o=="H")&(data_base$m==TRUE))]="To_P"
          data_base$o[which((data_base$o=="P")&(data_base$m==TRUE))]="To_H"

          ##### Print simulation progress
          print(paste("         --> Migration: ",nrow(subset(data_base,data_base$m==TRUE)),"/",nrow(data_base),sep=""))
        }
        
        ## (6) Temporality of migration and harvest
        ids_at_risk=c()
        ids_nomore_at_risk=c()
        # Non-null harvest rate with hunting season after breeding migrations
        if ((phi!=0)&(H_season!=0)&(eps!=0)){
            # -> Compute hunt risk status for migrant males in transition
            transit_toH=data_base$id[which(data_base$o=="To_H")]
            transit_toP=data_base$id[which(data_base$o=="To_P")]
            ids_at_risk=sample(as.character(transit_toH),round(H_season*length(transit_toH)),replace=FALSE)
            ids_nomore_at_risk=sample(as.character(transit_toP),round(H_season*length(transit_toP)),replace=FALSE)
            data_base$o[which(data_base$id%in%ids_at_risk)]="H"
            data_base$o[which(data_base$id%in%ids_nomore_at_risk)]="P"
        }
        # -> Save migrant indiv. within migration base
        migrant.info=subset(data_base,data_base$m==TRUE)
        migration_base=rbind(migration_base,cbind(t=rep(t,nrow(migrant.info)),NmP=rep(Nm_P_ante_migration,nrow(migrant.info)),NmH=rep(Nm_H_ante_migration,nrow(migrant.info)),migrant.info))
        
        ##### Update of output variables by time step for Migration
        data_base_out=subset(data_base,(data_base$a>0))
        ## absolute migration nb.
        out.M_PtoH[t+1]=nrow(subset(data_base_out,(data_base_out$o%in%c("H","To_H"))&(data_base_out$m==TRUE))) 
        out.M_HtoP[t+1]=nrow(subset(data_base_out,(data_base_out$o%in%c("P","To_P"))&(data_base_out$m==TRUE))) 
        out.M_PtoH_Risk[t+1]=ifelse(!((phi!=0)&(H_season!=0)&(eps!=0)),0,nrow(subset(data_base_out,(data_base_out$id%in%ids_at_risk)))) 
        out.M_HtoP_NoRisk[t+1]=ifelse(!((phi!=0)&(H_season!=0)&(eps!=0)),0,nrow(subset(data_base_out,(data_base_out$id%in%ids_nomore_at_risk)))) 
        ## migration rates
        out.r_PtoH[t+1]=nrow(subset(data_base_out,(data_base_out$o%in%c("H","To_H"))&(data_base_out$m==TRUE)))/Nm_P_ante_migration
        out.r_HtoP[t+1]=nrow(subset(data_base_out,(data_base_out$o%in%c("P","To_P"))&(data_base_out$m==TRUE)))/Nm_H_ante_migration 
        out.r_PtoH_Risk[t+1]=ifelse(!((phi!=0)&(H_season!=0)&(eps!=0)),0,nrow(subset(data_base_out,(data_base_out$id%in%ids_at_risk)))/Nm_P_ante_migration)  
        out.r_HtoP_NoRisk[t+1]=ifelse(!((phi!=0)&(H_season!=0)&(eps!=0)),0,nrow(subset(data_base_out,(data_base_out$id%in%ids_nomore_at_risk)))/Nm_H_ante_migration)
        
        ## (7) Harvesting 
        # Non-null harvest rate with hunting season after breeding migrations
        if (phi!=0){ 
          if (t>=TH){
            # -> Compute fate of all legal males for hunting after start of hunting period
            harv.deads=c() # current storage of id. of killed indiv.
            male_data_base=subset(data_base,(data_base$s=="M")&(data_base$o=="H")&(data_base$a>0))
            Nm=nrow(male_data_base)
            if(Nm==0){
              stop=TRUE
              print("!!! HARVEST ERROR: NO MALE !!!")
              break}
            for(i in 1:Nm){
              p_harvest=Harv.Mort(male_data_base$s[i],male_data_base$o[i],male_data_base$y[i],H,H_max,y_min,phi)
              harv.death=rbinom(1,1,p_harvest) # 1 random value in [0,1] with proba. of p_harvest
              # -> Update and save harvest
              if(harv.death==1){
                killed.id=male_data_base$id[i]
                harv.deads=c(harv.deads,killed.id)
                H=H+1 # storage of nb. of males killed by time step (= 1 year = 1 hunting season) to verify respect of quota when applied
              }
            }
            # -> Save killed indiv. within harvest base
            killed.info=subset(data_base,data_base$id%in%harv.deads)
            harvest_base=rbind(harvest_base,cbind(t=rep(t,nrow(killed.info)),killed.info))
            # -> Remove killed indiv. from metapop. 
            N_updated=nrow(data_base)
            data_base=subset(data_base,!(data_base$id%in%harv.deads)) # removal of killed indiv.
            
            ##### Print simulation progress
            print(paste("         --> Harvest: ",H,"/",N_updated,"=",nrow(data_base),sep=""))
          }
        }
        
        ## (8) Siring competition
        # -> Update origin status of males in transition so that they compete in appropriate siring competition
        if (eps!=0){
          data_base$o[which(data_base$o=="To_H")]="H"
          data_base$o[which(data_base$o=="To_P")]="P"
        }
        # -> Compute fate of all mature males for siring next yearling indiv.
        for (origins in c("P","H")){
          origin_male_data_base=subset(data_base,(data_base$s=="M")&(data_base$o==origins)&(data_base$a>1))
          Nm=nrow(origin_male_data_base) 
          if(Nm==0){
            stop=TRUE
            print("!!! SIRING ERROR: NO MALE !!!")
            break}
          rank=rep(0,Nm) # classification of mature males in same population of reproductive ewe according to horn length
          for(j in 1:Nm){
            p_siring=Male.Repro(origin_male_data_base$y[j],gamma)
            rank[j]=p_siring
          }
          unity=rep(1,length(rank))
          siring_weights=pmax(unity,round(rank*1.8)) # weights of male mating opportunities with horn length dependent rank
          if (origins=="P"){
            siring_male_phenotype_P=subset(origin_male_data_base,select=c('id','y','y_inf','bv_f','bv_m','m')) # save phenotype of male competitors for siring
            siring_male_representation_P=rep(origin_male_data_base$id,siring_weights) # representation of male mating opportunities with horn length dependent rank
          }else{
            siring_male_phenotype_H=subset(origin_male_data_base,select=c('id','y','y_inf','bv_f','bv_m','m')) # save phenotype of male competitors for siring
            siring_male_representation_H=rep(origin_male_data_base$id,siring_weights) # representation of male mating opportunities with horn length dependent rank
          }
        }
        
        ## (9) Natural mortality 
        # -> Compute fate of all indiv. for natural mortality
        N_updated=nrow(data_base) # metapop. size
        nat.deads=c() # current storage of id. of nat. deads 
        for(i in 1:N_updated){
          if (data_base$a[i]>0){
            p_dying=Nat.Mort(data_base$a[i],data_base$s[i],mu_females,mu_males,max_a_females,max_a_males)
            nat.death=rbinom(1,1,p_dying) # 1 random value in [0,1] with proba. of p_dying
            if(nat.death==1){
              deads.id=data_base$id[i]
              nat.deads=c(nat.deads,deads.id)
            }
          }
      }
      # -> Remove natural dead indiv. from metapop.
      N_updated=nrow(subset(data_base,(data_base$a>0)))
      data_base=subset(data_base,!(data_base$id%in%nat.deads))

      ##### Print simulation progress
      print(paste("         --> Natural mortality: ",length(nat.deads),"/",N_updated,"=",nrow(data_base),sep=""))
      
      # -> Reinitialisation of annual migration status
      if (eps!=0){
        data_base$o[which((data_base$m==TRUE)&(data_base$o=="P"))]="H_temp"
        data_base$o[which((data_base$m==TRUE)&(data_base$o=="H"))]="P"
        data_base$o[which(data_base$o=="H_temp")]="H"
        data_base$m=FALSE
      }
      
      ##### Update of output variables by time step
      out.t[t+1]=t
      data_base_out=subset(data_base,(data_base$a>0))
      
      ### Average output values
      out.a[t+1]=mean(data_base_out$a) 
      out.aH[t+1]=mean(data_base_out$a[which(data_base_out$o=="H")]) 
      out.aP[t+1]=mean(data_base_out$a[which(data_base_out$o=="P")]) 
      out.am[t+1]=mean(data_base_out$a[which(data_base_out$s=="M")]) 
      out.amH[t+1]=mean(data_base_out$a[which((data_base_out$o=="H")&(data_base_out$s=="M"))]) 
      out.amP[t+1]=mean(data_base_out$a[which((data_base_out$o=="P")&(data_base_out$s=="M"))]) 
      out.yf[t+1]=mean(data_base_out$y[which(data_base_out$s=="F")]) 
      out.yfH[t+1]=mean(data_base_out$y[which((data_base_out$o=="H")&(data_base_out$s=="F"))]) 
      out.yfP[t+1]=mean(data_base_out$y[which((data_base_out$o=="P")&(data_base_out$s=="F"))]) 
      out.ym[t+1]=mean(data_base_out$y[which(data_base_out$s=="M")]) 
      out.ymH[t+1]=mean(data_base_out$y[which((data_base_out$o=="H")&(data_base_out$s=="M"))]) 
      out.ymP[t+1]=mean(data_base_out$y[which((data_base_out$o=="P")&(data_base_out$s=="M"))]) 
      out.ym_six[t+1]=mean(data_base_out$y[which((data_base_out$a==6)&(data_base_out$s=="M"))]) 
      out.ymH_six[t+1]=mean(data_base_out$y[which((data_base_out$a==6)&(data_base_out$o=="H")&(data_base_out$s=="M"))]) 
      out.ymP_six[t+1]=mean(data_base_out$y[which((data_base_out$a==6)&(data_base_out$o=="P")&(data_base_out$s=="M"))]) 
      out.ym_sev[t+1]=mean(data_base_out$y[which((data_base_out$a==7)&(data_base_out$s=="M"))]) 
      out.ymH_sev[t+1]=mean(data_base_out$y[which((data_base_out$a==7)&(data_base_out$o=="H")&(data_base_out$s=="M"))]) 
      out.ymP_sev[t+1]=mean(data_base_out$y[which((data_base_out$a==7)&(data_base_out$o=="P")&(data_base_out$s=="M"))]) 
      out.ym_sevandolder[t+1]=mean(data_base_out$y[which((data_base_out$a>=7)&(data_base_out$s=="M"))]) 
      out.ymH_sevandolder[t+1]=mean(data_base_out$y[which((data_base_out$a>=7)&(data_base_out$o=="H")&(data_base_out$s=="M"))]) 
      out.ymP_sevandolder[t+1]=mean(data_base_out$y[which((data_base_out$a>=7)&(data_base_out$o=="P")&(data_base_out$s=="M"))]) 
      out.yf_inf[t+1]=mean(data_base_out$y_inf[which(data_base_out$s=="F")]) 
      out.yfH_inf[t+1]=mean(data_base_out$y_inf[which((data_base_out$o=="H")&(data_base_out$s=="F"))]) 
      out.yfP_inf[t+1]=mean(data_base_out$y_inf[which((data_base_out$o=="P")&(data_base_out$s=="F"))]) 
      out.ym_inf[t+1]=mean(data_base_out$y_inf[which(data_base_out$s=="M")]) 
      out.ymH_inf[t+1]=mean(data_base_out$y_inf[which((data_base_out$o=="H")&(data_base_out$s=="M"))]) 
      out.ymP_inf[t+1]=mean(data_base_out$y_inf[which((data_base_out$o=="P")&(data_base_out$s=="M"))]) 
      out.bvf_f[t+1]=mean(data_base_out$bv_f[which((data_base_out$s=="F"))]) 
      out.bvfH_f[t+1]=mean(data_base_out$bv_f[which((data_base_out$o=="H")&(data_base_out$s=="F"))]) 
      out.bvfP_f[t+1]=mean(data_base_out$bv_f[which((data_base_out$o=="P")&(data_base_out$s=="F"))])
      out.bvf_m[t+1]=mean(data_base_out$bv_m[which((data_base_out$s=="F"))]) 
      out.bvfH_m[t+1]=mean(data_base_out$bv_m[which((data_base_out$o=="H")&(data_base_out$s=="F"))]) 
      out.bvfP_m[t+1]=mean(data_base_out$bv_m[which((data_base_out$o=="P")&(data_base_out$s=="F"))]) 
      out.bvm_f[t+1]=mean(data_base_out$bv_f[which((data_base_out$s=="M"))]) 
      out.bvmH_f[t+1]=mean(data_base_out$bv_f[which((data_base_out$o=="H")&(data_base_out$s=="M"))]) 
      out.bvmP_f[t+1]=mean(data_base_out$bv_f[which((data_base_out$o=="P")&(data_base_out$s=="M"))]) 
      out.bvm_m[t+1]=mean(data_base_out$bv_m[which((data_base_out$s=="M"))])
      out.bvmH_m[t+1]=mean(data_base_out$bv_m[which((data_base_out$o=="H")&(data_base_out$s=="M"))]) 
      out.bvmP_m[t+1]=mean(data_base_out$bv_m[which((data_base_out$o=="P")&(data_base_out$s=="M"))]) 
      out.bvf_f_yearl[t+1]=mean(data_base_out$bv_f[which((data_base_out$s=="F")&(data_base_out$a==1))]) 
      out.bvfH_f_yearl[t+1]=mean(data_base_out$bv_f[which((data_base_out$o=="H")&(data_base_out$s=="F")&(data_base_out$a==1))]) 
      out.bvfP_f_yearl[t+1]=mean(data_base_out$bv_f[which((data_base_out$o=="P")&(data_base_out$s=="F")&(data_base_out$a==1))])
      out.bvf_m_yearl[t+1]=mean(data_base_out$bv_m[which((data_base_out$s=="F")&(data_base_out$a==1))]) 
      out.bvfH_m_yearl[t+1]=mean(data_base_out$bv_m[which((data_base_out$o=="H")&(data_base_out$s=="F")&(data_base_out$a==1))]) 
      out.bvfP_m_yearl[t+1]=mean(data_base_out$bv_m[which((data_base_out$o=="P")&(data_base_out$s=="F")&(data_base_out$a==1))]) 
      out.bvm_f_yearl[t+1]=mean(data_base_out$bv_f[which((data_base_out$s=="M")&(data_base_out$a==1))]) 
      out.bvmH_f_yearl[t+1]=mean(data_base_out$bv_f[which((data_base_out$o=="H")&(data_base_out$s=="M")&(data_base_out$a==1))]) 
      out.bvmP_f_yearl[t+1]=mean(data_base_out$bv_f[which((data_base_out$o=="P")&(data_base_out$s=="M")&(data_base_out$a==1))]) 
      out.bvm_m_yearl[t+1]=mean(data_base_out$bv_m[which((data_base_out$s=="M")&(data_base_out$a==1))])
      out.bvmH_m_yearl[t+1]=mean(data_base_out$bv_m[which((data_base_out$o=="H")&(data_base_out$s=="M")&(data_base_out$a==1))]) 
      out.bvmP_m_yearl[t+1]=mean(data_base_out$bv_m[which((data_base_out$o=="P")&(data_base_out$s=="M")&(data_base_out$a==1))]) 
      
      ### Variance output values
      out.sd.a[t+1]=sd(data_base_out$a) 
      out.sd.aH[t+1]=sd(data_base_out$a[which(data_base_out$o=="H")]) 
      out.sd.aP[t+1]=sd(data_base_out$a[which(data_base_out$o=="P")]) 
      out.sd.yf[t+1]=sd(data_base_out$y[which(data_base_out$s=="F")]) 
      out.sd.yfH[t+1]=sd(data_base_out$y[which((data_base_out$o=="H")&(data_base_out$s=="F"))]) 
      out.sd.yfP[t+1]=sd(data_base_out$y[which((data_base_out$o=="P")&(data_base_out$s=="F"))]) 
      out.sd.ym[t+1]=sd(data_base_out$y[which(data_base_out$s=="M")]) 
      out.sd.ymH[t+1]=sd(data_base_out$y[which((data_base_out$o=="H")&(data_base_out$s=="M"))]) 
      out.sd.ymP[t+1]=sd(data_base_out$y[which((data_base_out$o=="P")&(data_base_out$s=="M"))]) 
      out.sd.yf_inf[t+1]=sd(data_base_out$y_inf[which(data_base_out$s=="F")]) 
      out.sd.yfH_inf[t+1]=sd(data_base_out$y_inf[which((data_base_out$o=="H")&(data_base_out$s=="F"))]) 
      out.sd.yfP_inf[t+1]=sd(data_base_out$y_inf[which((data_base_out$o=="P")&(data_base_out$s=="F"))]) 
      out.sd.ym_inf[t+1]=sd(data_base_out$y_inf[which(data_base_out$s=="M")]) 
      out.sd.ymH_inf[t+1]=sd(data_base_out$y_inf[which((data_base_out$o=="H")&(data_base_out$s=="M"))])
      out.sd.ymP_inf[t+1]=sd(data_base_out$y_inf[which((data_base_out$o=="P")&(data_base_out$s=="M"))]) 
      out.sd.bvf_f[t+1]=sd(data_base_out$bv_f[which(data_base_out$s=="F")]) 
      out.sd.bvfH_f[t+1]=sd(data_base_out$bv_f[which((data_base_out$o=="H")&(data_base_out$s=="F"))]) 
      out.sd.bvfP_f[t+1]=sd(data_base_out$bv_f[which((data_base_out$o=="P")&(data_base_out$s=="F"))]) 
      out.sd.bvf_m[t+1]=sd(data_base_out$bv_m[which(data_base_out$s=="F")])
      out.sd.bvfH_m[t+1]=sd(data_base_out$bv_m[which((data_base_out$o=="H")&(data_base_out$s=="F"))]) 
      out.sd.bvfP_m[t+1]=sd(data_base_out$bv_m[which((data_base_out$o=="P")&(data_base_out$s=="F"))])
      out.sd.bvm_f[t+1]=sd(data_base_out$bv_f[which(data_base_out$s=="M")]) 
      out.sd.bvmH_f[t+1]=sd(data_base_out$bv_f[which((data_base_out$o=="H")&(data_base_out$s=="M"))]) 
      out.sd.bvmP_f[t+1]=sd(data_base_out$bv_f[which((data_base_out$o=="P")&(data_base_out$s=="M"))]) 
      out.sd.bvm_m[t+1]=sd(data_base_out$bv_m[which(data_base_out$s=="M")]) 
      out.sd.bvmH_m[t+1]=sd(data_base_out$bv_m[which((data_base_out$o=="H")&(data_base_out$s=="M"))]) 
      out.sd.bvmP_m[t+1]=sd(data_base_out$bv_m[which((data_base_out$o=="P")&(data_base_out$s=="M"))]) 
      
      out.sd.bvf_f_yearl[t+1]=sd(data_base_out$bv_f[which((data_base_out$s=="F")&(data_base_out$a==1))]) 
      out.sd.bvfH_f_yearl[t+1]=sd(data_base_out$bv_f[which((data_base_out$o=="H")&(data_base_out$s=="F")&(data_base_out$a==1))]) 
      out.sd.bvfP_f_yearl[t+1]=sd(data_base_out$bv_f[which((data_base_out$o=="P")&(data_base_out$s=="F")&(data_base_out$a==1))]) 
      out.sd.bvf_m_yearl[t+1]=sd(data_base_out$bv_m[which((data_base_out$s=="F")&(data_base_out$a==1))])
      out.sd.bvfH_m_yearl[t+1]=sd(data_base_out$bv_m[which((data_base_out$o=="H")&(data_base_out$s=="F")&(data_base_out$a==1))]) 
      out.sd.bvfP_m_yearl[t+1]=sd(data_base_out$bv_m[which((data_base_out$o=="P")&(data_base_out$s=="F")&(data_base_out$a==1))])
      out.sd.bvm_f_yearl[t+1]=sd(data_base_out$bv_f[which((data_base_out$s=="M")&(data_base_out$a==1))]) 
      out.sd.bvmH_f_yearl[t+1]=sd(data_base_out$bv_f[which((data_base_out$o=="H")&(data_base_out$s=="M")&(data_base_out$a==1))]) 
      out.sd.bvmP_f_yearl[t+1]=sd(data_base_out$bv_f[which((data_base_out$o=="P")&(data_base_out$s=="M")&(data_base_out$a==1))]) 
      out.sd.bvm_m_yearl[t+1]=sd(data_base_out$bv_m[which((data_base_out$s=="M")&(data_base_out$a==1))]) 
      out.sd.bvmH_m_yearl[t+1]=sd(data_base_out$bv_m[which((data_base_out$o=="H")&(data_base_out$s=="M")&(data_base_out$a==1))]) 
      out.sd.bvmP_m_yearl[t+1]=sd(data_base_out$bv_m[which((data_base_out$o=="P")&(data_base_out$s=="M")&(data_base_out$a==1))]) 
      
      out.sd.ym_six[t+1]=sd(data_base_out$y[which((data_base_out$a==6)&(data_base_out$s=="M"))]) 
      out.sd.ymH_six[t+1]=sd(data_base_out$y[which((data_base_out$a==6)&(data_base_out$o=="H")&(data_base_out$s=="M"))]) 
      out.sd.ymP_six[t+1]=sd(data_base_out$y[which((data_base_out$a==6)&(data_base_out$o=="P")&(data_base_out$s=="M"))]) 
      out.sd.ym_sev[t+1]=sd(data_base_out$y[which((data_base_out$a==7)&(data_base_out$s=="M"))]) 
      out.sd.ymH_sev[t+1]=sd(data_base_out$y[which((data_base_out$a==7)&(data_base_out$o=="H")&(data_base_out$s=="M"))]) 
      out.sd.ymP_sev[t+1]=sd(data_base_out$y[which((data_base_out$a==7)&(data_base_out$o=="P")&(data_base_out$s=="M"))]) 
      
      ### Counting output values
      out.H[t+1]=H # 
      out.prop_m_45[t+1]=round(nrow(subset(data_base_out,(data_base_out$a%in%c(4,5))&(data_base_out$s=="M")))/nrow(subset(data_base_out,(data_base_out$s=="M"))),digits=2) 
      out.prop_mH_45[t+1]=round(nrow(subset(data_base_out,(data_base_out$o=="H")&(data_base_out$a%in%c(4,5))&(data_base_out$s=="M")))/nrow(subset(data_base_out,(data_base_out$o=="H")&(data_base_out$s=="M"))),digits=2) 
      out.prop_mP_45[t+1]=round(nrow(subset(data_base_out,(data_base_out$o=="P")&(data_base_out$a%in%c(4,5))&(data_base_out$s=="M")))/nrow(subset(data_base_out,(data_base_out$o=="P")&(data_base_out$s=="M"))),digits=2) 
      out.struct_a[t+1]=round(nrow(subset(data_base_out,data_base_out$a>=7))/nrow(data_base_out),digits=2) 
      out.struct_aH[t+1]=round(nrow(subset(data_base_out,(data_base_out$o=="H")&(data_base_out$a>=7)))/nrow(subset(data_base_out,(data_base_out$o=="H"))),digits=2) 
      out.struct_aP[t+1]=round(nrow(subset(data_base_out,(data_base_out$o=="P")&(data_base_out$a>=7)))/nrow(subset(data_base_out,(data_base_out$o=="P"))),digits=2) 
      out.q[t+1]=nrow(subset(data_base_out,(data_base_out$s=="M")&(data_base_out$y>y_min)))/nrow(subset(data_base_out,(data_base_out$s=="M"))) 
      out.qH[t+1]=nrow(subset(data_base_out,(data_base_out$o=="H")&(data_base_out$s=="M")&(data_base_out$y>y_min)))/nrow(subset(data_base_out,(data_base_out$o=="H")&(data_base_out$s=="M"))) 
      out.qP[t+1]=nrow(subset(data_base_out,(data_base_out$o=="P")&(data_base_out$s=="M")&(data_base_out$y>y_min)))/nrow(subset(data_base_out,(data_base_out$o=="P")&(data_base_out$s=="M"))) 
      out.sr[t+1]=nrow(subset(data_base_out,(data_base_out$s=="F")))/nrow(subset(data_base_out,(data_base_out$s=="M"))) 
      out.srH[t+1]=nrow(subset(data_base_out,(data_base_out$o=="H")&(data_base_out$s=="F")))/nrow(subset(data_base_out,(data_base_out$o=="H")&(data_base_out$s=="M"))) 
      out.srP[t+1]=nrow(subset(data_base_out,(data_base_out$o=="P")&(data_base_out$s=="F")))/nrow(subset(data_base_out,(data_base_out$o=="P")&(data_base_out$s=="M"))) 
      out.N[t+1]=nrow(data_base_out) 
      out.NH[t+1]=nrow(subset(data_base_out,(data_base_out$o=="H"))) 
      out.NP[t+1]=nrow(subset(data_base_out,(data_base_out$o=="P"))) 
      out.Nm[t+1]=nrow(subset(data_base_out,(data_base_out$s=="M"))) 
      out.NmH[t+1]=nrow(subset(data_base_out,(data_base_out$o=="H")&(data_base_out$s=="M"))) 
      out.NmP[t+1]=nrow(subset(data_base_out,(data_base_out$o=="P")&(data_base_out$s=="M")))
      out.Nf[t+1]=nrow(subset(data_base_out,(data_base_out$s=="F"))) 
      out.NfH[t+1]=nrow(subset(data_base_out,(data_base_out$o=="H")&(data_base_out$s=="F")))
      out.NfP[t+1]=nrow(subset(data_base_out,(data_base_out$o=="P")&(data_base_out$s=="F"))) 
      
      ##### Simulation breaks in case of empty population
      if(out.N[t+1]==0){
          stop=TRUE
          print("!!! FINAL ERROR: EMPTY METAPOP. !!!")
          break
        }else if(out.NP[t+1]==0){
            stop=TRUE
            print("!!! FINAL ERROR: EMPTY PROTECTED POP. !!!")
            break
        }else if(out.NH[t+1]==0){
            stop=TRUE
            print("!!! FINAL ERROR: EMPTY EXPLOITED POP. !!!")
          break}
      if (stop){break} # stop replicate of scenario because of simulation error
    
    } # end loop of time steps of each replicate of each scenario
  
    #######################################################################################################################################
    ########################################         Results savings                       #################################################    
    #######################################################################################################################################
    
    out.id=gsub("\\s", "_", Sys.time())
    output_directory=file.path(getwd(),paste("/Outputs/Outputs_General/",scen.label,sep=""))
    if(!dir.exists(output_directory)){dir.create(output_directory)}
    new_output_directory=file.path(paste(output_directory,out.id,sep="/"))
    dir.create(new_output_directory)
    
    ##### Saving of output variables in external files
    output_data=data.frame(t=out.t,
                           ### Average output values
                           a=out.a,aH=out.aH,aP=out.aP,
                           am=out.am,amH=out.amH,amP=out.amP,
                           yf=out.yf,yfH=out.yfH,yfP=out.yfP,
                           ym=out.ym,ymH=out.ymH,ymP=out.ymP,
                           ym_six=out.ym_six,ymH_six=out.ymH_six,ymP_six=out.ymP_six,
                           ym_sev=out.ym_sev,ymH_sev=out.ymH_sev,ymP_sev=out.ymP_sev,
                           ym_sevandolder=out.ym_sevandolder,ymH_sevandolder=out.ymH_sevandolder,ymP_sevandolder=out.ymP_sevandolder,
                           yf_inf=out.yf_inf,yfP_inf=out.yfP_inf,yfH_inf=out.yfH_inf,
                           ym_inf=out.ym_inf,ymP_inf=out.ymP_inf,ymH_inf=out.ymH_inf,
                           bvf_f=out.bvf_f,bvfH_f=out.bvfH_f,bvfP_f=out.bvfP_f,
                           bvf_m=out.bvf_m,bvfH_m=out.bvfH_m,bvfP_m=out.bvfP_m, 
                           bvm_f=out.bvm_f,bvmH_f=out.bvmH_f,bvmP_f=out.bvmP_f,
                           bvm_m=out.bvm_m,bvmH_m=out.bvmH_m,bvmP_m=out.bvmP_m,
                           bvf_f_yearl=out.bvf_f_yearl,bvfH_f_yearl=out.bvfH_f_yearl,bvfP_f_yearl=out.bvfP_f_yearl,
                           bvf_m_yearl=out.bvf_m_yearl,bvfH_m_yearl=out.bvfH_m,bvfP_m_yearl=out.bvfP_m_yearl, 
                           bvm_f_yearl=out.bvm_f_yearl,bvmH_f_yearl=out.bvmH_f_yearl,bvmP_f_yearl=out.bvmP_f_yearl,
                           bvm_m_yearl=out.bvm_m_yearl,bvmH_m_yearl=out.bvmH_m_yearl,bvmP_m_yearl=out.bvmP_m_yearl,
                           
                           ### Variance output values
                           a.sd=out.sd.a,aH.sd=out.sd.aH,aP.sd=out.sd.aP,
                           yf.sd=out.sd.yf,yfH.sd=out.sd.yfH,yfP.sd=out.sd.yfP,
                           ym.sd=out.sd.ym,ymH.sd=out.sd.ymH,ymP.sd=out.sd.ymP,
                           yf_inf.sd=out.sd.yf_inf,yfP_inf.sd=out.sd.yfP_inf,yfH_inf.sd=out.sd.yfH_inf,
                           ym_inf.sd=out.sd.ym_inf,ymP_inf.sd=out.sd.ymP_inf,ymH_inf.sd=out.sd.ymH_inf,
                           bvf_f.sd=out.sd.bvf_f,bvfH_f.sd=out.sd.bvfH_f,bvfP_f.sd=out.sd.bvfP_f,
                           bvf_m.sd=out.sd.bvf_m,bvfH_m.sd=out.sd.bvfH_m,bvfP_m.sd=out.sd.bvfP_m, 
                           bvm_f.sd=out.sd.bvm_f,bvmH_f.sd=out.sd.bvmH_f,bvmP_f.sd=out.sd.bvmP_f,
                           bvm_m.sd=out.sd.bvm_m,bvmH_m.sd=out.sd.bvmH_m,bvmP_m.sd=out.sd.bvmP_m,
                           bvf_f.sd_yearl=out.sd.bvf_f_yearl,bvfH_f.sd_yearl=out.sd.bvfH_f_yearl,bvfP_f.sd_yearl=out.sd.bvfP_f_yearl,
                           bvf_m.sd_yearl=out.sd.bvf_m_yearl,bvfH_m.sd_yearl=out.sd.bvfH_m_yearl,bvfP_m.sd_yearl=out.sd.bvfP_m_yearl, 
                           bvm_f.sd_yearl=out.sd.bvm_f_yearl,bvmH_f.sd_yearl=out.sd.bvmH_f_yearl,bvmP_f.sd_yearl=out.sd.bvmP_f_yearl,
                           bvm_m.sd_yearl=out.sd.bvm_m_yearl,bvmH_m.sd_yearl=out.sd.bvmH_m_yearl,bvmP_m.sd_yearl=out.sd.bvmP_m_yearl,
                           ym_six.sd=out.sd.ym_six,ymH_six.sd=out.sd.ymH_six,ymP_six.sd=out.sd.ymP_six,
                           ym_sev.sd=out.sd.ym_sev,ymH_sev.sd=out.sd.ymH_sev,ymP_sev.sd=out.sd.ymP_sev,
                           
                           ### Counting output values
                           H=out.H,
                           prop_m_45=out.prop_m_45,prop_mH_45=out.prop_mH_45,prop_mP_45=out.prop_mP_45,
                           struct_a=out.struct_a,struct_aH=out.struct_aH,struct_aP=out.struct_aP,
                           q=out.q,qH=out.qH,qP=out.qP,
                           sr=out.sr,srH=out.srH,srP=out.srP,
                           N=out.N,NH=out.NH,NP=out.NP,
                           Nm=out.Nm,NmH=out.NmH,NmP=out.NmP,Nf=out.Nf,NfH=out.NfH,NfP=out.NfP,
                           M_PtoH=out.M_PtoH,M_HtoP=out.M_HtoP,M_PtoH_Risk=out.M_PtoH_Risk,M_HtoP_NoRisk=out.M_HtoP_NoRisk,
                           r_PtoH=out.r_PtoH,r_HtoP=out.r_HtoP,r_PtoH_Risk=out.r_PtoH_Risk,r_HtoP_NoRisk=out.r_HtoP_NoRisk)
    
    output_file=paste("OutputFile",scen.label,out.id,sep="_")
    output_file=paste(output_file,".xlsx",sep="")
    write.xlsx(output_data,paste(new_output_directory,output_file,sep="/"),colNames=TRUE,rowNames=FALSE,startRow=1)
    
    ##### Saving of initial population in external files
    init_file=paste("InitFile",scen.label,out.id,sep="_")
    init_file=paste(init_file,".xlsx",sep="")
    write.xlsx(init_base,paste(new_output_directory,init_file,sep="/"),colNames=TRUE,rowNames=FALSE,startRow=1)
    
    ##### Saving of final population in external files
    data_file=paste("DataFile",scen.label,out.id,sep="_")
    data_file=paste(data_file,".xlsx",sep="")
    write.xlsx(data_base,paste(new_output_directory,data_file,sep="/"),colNames=TRUE,rowNames=FALSE,startRow=1)
    
    ##### Saving of migration in external files
    migration_file=paste("MigrationFile",scen.label,out.id,sep="_")
    migration_file=paste(migration_file,".xlsx",sep="")
    write.xlsx(migration_base,paste(new_output_directory,migration_file,sep="/"),colNames=TRUE,rowNames=FALSE,startRow=1)
    
    ##### Saving of harvest in external files
    harvest_file=paste("HarvestFile",scen.label,out.id,sep="_")
    harvest_file=paste(harvest_file,".xlsx",sep="")
    write.xlsx(harvest_base,paste(new_output_directory,harvest_file,sep="/"),colNames=TRUE,rowNames=FALSE,startRow=1)

    ##### Saving of recruitment in external files
    recruitment_file=paste("RecruitmentFile",scen.label,out.id,sep="_")
    recruitment_file=paste(recruitment_file,".xlsx",sep="")
    write.xlsx(recruitment_base,paste(new_output_directory,recruitment_file,sep="/"),colNames=TRUE,rowNames=FALSE,startRow=1)
    
  } # end loop of all replicates of each scenario
} # end loop of all possible scenarios

