#######################################################################################################################################
########################################                                                ################################################
########################################                 Simulation Model                #################################################     
########################################          Definition of Useful Functions         ################################################
########################################                                                 ################################################
#######################################################################################################################################

rm(list=ls()) 

#######################################################################################################################################
########################################          Definitions              #################################################     
#######################################################################################################################################

######## 1) Aging and growing ########
## >>> Length.Growth:: horn length growth
## Individual inputs:: s:sex; y: horn size; y_inf: asymptotic horn length
## Parameters inputs:: rho: growth rate
## Output:: y_new: new horn size
Length.Growth=function(s,y,y_inf,rho_females,rho_males){
  # Select sex-specific growth rate
  rho=ifelse(data_base$s[i]=="F",rho_females,rho_males)
  # horn growth increment from von Bertalanffy growth equation
  delta_y=ifelse((y_inf-y)>0,(y_inf-y)*(1-exp(-rho)),0) # null increment delta_y if y>=y_inf
  y_new=round(y+delta_y,digits=2)
  return(y_new)
}

######## 2) Natural mortality ########
## >>> Nat.Mort: natural mortality
## Individual inputs:: a: age; s: sex
## Parameters inputs:: mu_females/males: sex-specific death rate; max_a_females/males: sex-specific longevity
## Output:: probability of natural mortality
Nat.Mort=function(a,s,mu_females,mu_males,max_a_females,max_a_males){
  # Select sex-specific parameters
  if(s=="F"){
    mu=mu_females
    max_a=max_a_females
    }else{
    mu=mu_males
    max_a=max_a_males}
  # Age-dependent natural mortality
  p_dying=1-1/(1+exp(-(mu[1]+mu[2]*(a-8)+mu[3]*(a-8)^2))) # inverse logit transformation
  p_dying=ifelse(a>max_a,1,p_dying) # longevity threshold
  return(p_dying)
}

########  3) Harvesting ########
## >>> Harv.Mort: harvest mortality
## Individual inputs:: y: horn size; s: sex; o: origin
## Populational inputs:: H: nb. of males harvested during ongoing hunting season
## Regulation inputs:: H_max: maximum legal harvest, y_min: minimum legal trophy size
## Parameters inputs:: phi: harvest rate of legal rams
## Output:: probability of being killed
Harv.Mort=function(s,o,y,H,H_max,y_min,phi){
  # Null probability for non-legal males, protected males or females
  p_harvest=0
  # Harvest probability horn size dependent for legal males from harvested population
  if((s=="M")&(o=="H")){
    # Pre-compute the legality according to quota and minimum legal size
    if((H<H_max)&(y>y_min)){
      p_harvest=phi
      } 
  }
  return(p_harvest)
}

######## 4) Migration ########
# >>> Breed.Mig: rams breeding migration
## Individual inputs:: a: age
## Populational inputs:: Nm_comp_scale: scaled nb. of competitors; older_m: age of older male (to take into account age structure)
## Parameters inputs:: lambda: age and nb. of competitors influences; sr: female per male sex-ratio
## Output:: probability of migrating
Breed.Mig=function(a,lambda,sr,Nm_comp_scale,older_m){
    ages=seq(1,older_m,by=1)
    p_migrate=(1/(1+exp(-(lambda[1]+lambda[2]*bs(ages,degree=2)[,1]+lambda[3]*bs(ages,degree=2)[,2]+lambda[4]*sr+lambda[5]*Nm_comp_scale))))
    p_migrate=p_migrate[a]
  return(p_migrate)
}

######## 5) Reproduction ########
## >>> Female.Repro: probability of a female to produce a yearling
## Individual inputs:: a: age
## Populational inputs:: N_scale: scaled pop. size
## Parameters inputs:: theta: yearling production rate (depending on age and pop. size)
## Output:: probability of producing a yearling
Female.Repro=function(a,theta,N_scale){
  # Lambing for mature and non-mother females
  p_birth=1/(1+exp(-(theta[1]+theta[2]*a+theta[3]*a^2+theta[4]*N_scale))) # inverse logit transformation
  return(p_birth)
}

## >>> Male.Repro: probability of a male to sire lambs
## Individual inputs:: y: horn size
## Parameters inputs:: gamma: horn length influence
## Output:: probability of siring
Male.Repro=function(y,gamma){
  p_siring=1/(1+exp(-(gamma[1]+gamma[2]*y))) # inverse logit transformation
  return(p_siring)
}

######## 6) Offspring ########
## >>> Asymt.Length: asymptotic horn length of new yearling indiv.
## Individual inputs:: s: sex; bv_new_f,bv_new_m: offspring breeding value for sex-specific horn length
## Populational inputs:: y_inf_mean_females/males: populational mean asymptotic horn length by sex; E_ff, E_mm: intra-sexual environmental variance
## Output:: sex-specific asymptotic horn length
Asympt.Length=function(s,bv_new_f,bv_new_m,y_inf_mean_females,y_inf_mean_males,E_ff,E_mm){
  # Compute sex-specific asymptotic horn length
  if(s=="F"){y_inf=rnorm(1,y_inf_mean_females+bv_new_f,sqrt(E_ff))} 
  else if(s=="M"){y_inf=rnorm(1,y_inf_mean_males+bv_new_m,sqrt(E_mm))}
  else {print("!!! SEX ERROR !!!")}
  return(round(y_inf,digits=2))
}

## >>> Breeding.Value: breeding value of new yearling indiv.
## Individual inputs:: bv_father_f, bv_father_m, bv_mother_f, bv_mother_m: breeding value of mother and father for sex-specific horn length (in cm)
## Parameters inputs:: G_ff, G_fm, G_mm: intra and inter-sexual genetic variance
## Output:: breeding value for sex-specific horn length (bv_new=c(bv_new_f,bv_new_m))
Breeding.Value=function(bv_mother_f,bv_mother_m,bv_father_f,bv_father_m,G_mm,G_mf,G_ff){
  bv_mother=c(bv_mother_f,bv_mother_m)
  bv_father=c(bv_father_f,bv_father_m) 
  G=matrix(c(G_ff,G_mf,G_mf,G_mm),2,2)
  bv_new=rmvnorm(1,(bv_mother+bv_father)/2,G) # generates data from the multivariate normal distribution from mean vector and covariance matrix
  return(bv_new)
}

## >>> Sex.Assign: sex assignation
## Parameters inputs:: sigma: sex ratio
## Output:: sex of new born lambs
Sex.Assign=function(sigma){
  # 1 random value in [0,1] with proba. of sigma
  p_sex=rbinom(1,1,sigma) 
  if(p_sex==1){sex="F"}
  else{sex="M"}
  return(sex)
}
