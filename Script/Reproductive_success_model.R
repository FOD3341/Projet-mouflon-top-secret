# read packages

if (!require("brms")) install.packages("brms")

if (!require("ggplot2")) install.packages("ggplot2")

if (!require("tidyverse")) install.packages("tidyverse") 

if (!require("tidybayes")) install.packages("tidybayes")

if (!require("hrbrthemes")) install.packages("hrbrthemes")

# read data

data <- read.csv("Cleaned_data/cleaned_data.csv", header = T, sep = ",")
horn_data <- read.csv("Cleaned_data/adjhlg.csv", header = T, sep = ",")
ped <- read.csv("raw_data/Ram_Pedigree_W_Parents_Age.csv",header =T, sep =";")
LRS <- read.csv("raw_data/RamMtn_LifetimeReproSuccess.csv",header =T, sep =";")

# rename columns for merge
data <- data %>%
  rename(age = Age)
# merging the two data frame and keeping only one value for each value of interest

data <- data %>%
  left_join(horn_data,by = c("ID","age"))
  
  
# select only one line per year

data <- data %>%
  group_by(ID, Year) %>%  # Grouper by ID and year
  summarise( 
    age = first(age),
    density = first(density),
    RS = first(RS),
    hlg12 = first(adjhlg[JJ == 12]),  # Sélectionner adjhlg pour JJ = 12
    hlg114 = first(adjhlg[JJ == 114]),
    EBV = first(mu),
    age_class = first(age_class),
    age_scaled = first(age_scaled),
    density_scaled = first(density_scaled),
    mu_scaled = first(mu_scaled),
    cohort = first(cohort),
    wtd12 = first(wtd12)
  )

### data visualization

# RS by age
ggplot(data, aes(x = factor(age), y = RS)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) +
  labs(title = "RS by age",
       x = "age",
       y = "RS") +
  theme_ipsum()

# relationship between BV and reproductive success

# Plot with mu RS
ggplot(data, aes(x = EBV, y = RS)) +
  geom_point() +
  geom_smooth(method = "gam", color="red", fill="#69b3a2", se=TRUE) +  # Add regression lines
  labs(title = "Relationship Between EBV and RS",
       x = "EBV", 
       y = "RS") +
  theme_ipsum()

#Plot for horn_length 12 and reproductive succes
ggplot(data, aes(x = hlg12, y = RS)) +
  geom_point() +
  geom_smooth(method = "gam", color="red", fill="#69b3a2", se=TRUE) +  # Add regression lines
  labs(title = "Relationship Between hlg and RS",
       x = "horn_length", 
       y = "RS") +
  theme_ipsum()

#Plot for horn_length 114 and reproductive succes
ggplot(data, aes(x = hlg114, y = RS)) +
  geom_point() +
  geom_smooth(method = "gam", color="red", fill="#69b3a2", se=TRUE) +  # Add regression lines
  labs(title = "Relationship Between hlg and RS",
       x = "horn_length", 
       y = "RS") +
  theme_ipsum()


# relationship between density and reproductive success
# Plot with facets
ggplot(data, aes(x = density, y = RS)) +
  geom_point() +
  geom_smooth(method = "gam", color="red", fill="#69b3a2", se=TRUE) +
  labs(title = "Relationship Between density and RS",
       x = "Density", 
       y = "RS") +
  theme_ipsum()

# Reproductive success by Year

ggplot(data, aes(x = factor(Year), y = RS)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) +
  labs(title = "RS by Year",
       x = "Year",
       y = "RS") +
  theme(axis.text.x = element_text(size = 1)) +
  theme_ipsum()


# Reproductive success by cohort

ggplot(data, aes(x = factor(cohort), y = RS)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) +
  labs(title = "RS by cohort",
       x = "cohort",
       y = "RS") +
  theme(axis.text.x = element_text(size = 5)) +
  theme_ipsum()

#  Model horn length

prior <- c(
  prior(normal(-0.5, 0.5), class = "b", coef = "density"),
  prior(normal(0.5, 0.5), class = "b", coef = "age"),
#  prior(normal(-0.5, 0.5), class = "b", coef = "IageE2"),
  prior(normal(0.5, 0.5), class = "b", coef = "hlg12"),
  prior(normal(0.5, 1), class = "b", coef = "age:hlg12"),
#  prior(normal(-0.5, 1), class = "b", coef = "IageE2:hlg12"),
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "sd", group = "ID"),
  prior(normal(0, 1), class = "sd", group = "Year"),
  prior(normal(0, 1), class = "sd", group = "cohort")
)
+ I(age^2)*hlg12 + I(age^2)
model <- brm(
  RS ~ density + age  + hlg12 + hlg12*age  + (1 |  Year) + (1 |  ID) + (1 |  cohort) ,
  data = data,
  family = bernoulli(link = "logit"),
  prior = prior,
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 10),
  silent = 0,
  save_pars = save_pars(all = TRUE)
)

# Get posterior predictions
conditions <- expand.grid(
  age = 1:19, 
  hlg12 = c(18, 22, 26),  # Trois valeurs fixes de hlg12
  density = mean(data$density, na.rm = TRUE), # Une seule valeur moyenne de density
  cohort = "new",
  Year = "new",
  ID = "new"
)

preds <- posterior_epred(model, newdata = conditions,allow_new_levels = TRUE)

# Summarize predictions
pred_summary <- apply(preds, 2, function(x) c(mean = mean(x), lower = quantile(x, 0.025), upper = quantile(x, 0.975)))
pred_summary <- t(pred_summary)
colnames(pred_summary) <- c("mean", "lower", "upper")


#Convert predictions to a data frame and add row numbers
pred_df <- pred_summary %>%
  as.data.frame() %>%
  mutate(Rownumber = row_number()) 

# Add row numbers to conditions as well
conditions <- conditions %>%
  mutate(Rownumber = row_number()) 

# Perform the left join using the Rownumber
effect <- pred_df %>%
  left_join(conditions, by = "Rownumber")

plot(effect$age,effect$mean)


#Creating horn class
horn_class <- quantile(data$hlg12,na.rm = T)



# Create age classes based on fixed cutoffs
data <- data %>%
  filter(!is.na(hlg12))

data$hlg_class <- cut(data$hlg12,
                     breaks = c(0, 18, 24,40),  # Define the break points
                     labels = c("Petite cornes", "Moyenne cornes","Grosse cornes"), 
                     right = FALSE)  # right = FALSE means the interval is left-closed

effect <- as.data.frame(conditional_effects(model,"age:hlg12")[[1]])

effect$hlg_class <- factor(effect$hlg12, 
                           levels = c(17.5307079105918, 22.1390055802425, 26.7473032498933),
                           labels = c("Petite cornes", "Moyenne cornes", "Grosse cornes"))
# Assuming your data frame is called 'data' and contains columns 'RS', 'age', 'hlg12'

prop <- data %>%
  mutate(ntot = nrow(data)) %>%
  group_by(age, hlg_class) %>%
  mutate(proportion = n() / ntot) %>%
  distinct(age,hlg_class,proportion)
  


P <- ggplot() +
  geom_line(data = effect, aes(x = age, y = estimate__, colour = hlg_class)) +
  geom_ribbon(data = effect, aes(x = age, ymin = lower__, ymax = upper__, fill = factor(hlg_class)), alpha = 0.3) +
  labs(title = "Succès reproducteur en relation avec la taille des cornes et l'âge",
       x = "Âge",
       y = "Succès reproducteur") +
  theme_minimal() +
  scale_colour_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB"), 
                      name = "Classe de cornes",  # Customizing the legend title
                      labels = c("Petite cornes", "Moyenne cornes", "Grande cornes")) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB"), 
                    name = "Classe de cornes",  # Customizing the legend title
                    labels = c("Petite cornes", "Moyenne cornes", "Grande cornes"),
                    guide = guide_legend(override.aes = list(fill = c("#66C2A5", "#FC8D62", "#8DA0CB")))) # Ensures the same legend for fill

P + geom_bar(data = prop, aes(x = age, y = proportion, fill = factor(hlg_class)), 
             stat = "identity", position = "stack", alpha = 0.5, show.legend = FALSE) +
  scale_y_continuous(
    sec.axis = sec_axis(~., name = "proportion d'observation")  # Secondary axis for proportion
  ) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB"), 
name = "Classe de cornes",  # Customizing the legend title
labels = c("Petite cornes", "Moyenne cornes", "Grande cornes"),
guide = guide_legend(override.aes = list(fill = c("#66C2A5", "#FC8D62", "#8DA0CB"))))

  


  
  









conditional_effects(model)
pp_check(model)
summary(model)
plot(model)


# model 2 hlg age_class

prior2 <- c(
  prior(normal(-0.5, 0.5), class = "b", coef = "density"),
  prior(normal(0.5, 0.5), class = "b", coef = "age_classmadame"),
  prior(normal(0.5, 0.5), class = "b", coef = "age_classviellemadame"),
  prior(normal(0.5, 0.5), class = "b", coef = "hlg12"),
  prior(normal(0.5, 1), class = "b", coef = "hlg12:age_classmadame"),
  prior(normal(-0.5, 1), class = "b", coef = "hlg12:age_classviellemadame"),
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "sd", group = "ID"),
  prior(normal(0, 1), class = "sd", group = "Year"),
  prior(normal(0, 1), class = "sd", group = "cohort")
)

model2 <- brm(
  RS ~ density + hlg12 + age_class + age_class*hlg12 + (1 |  Year) + (1 |  ID) + (1 |  cohort),
  data = data,
  family = bernoulli(link = "logit"),
  prior = prior2,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 10),
  silent = 0,
  save_pars = save_pars(all = TRUE)
)

conditional_effects(model2)
pp_check(model2)
summary(model2)
plot(model2)


# model 3 EBV

prior3 <- c(
  prior(normal(-0.5, 0.5), class = "b", coef = "density"),
  prior(normal(0.5, 0.5), class = "b", coef = "age"),
  prior(normal(-0.5, 0.5), class = "b", coef = "IageE2"),
  prior(normal(0.5, 0.5), class = "b", coef = "EBV"),
  prior(normal(0.5, 1), class = "b", coef = "age:EBV"),
  prior(normal(-0.5, 1), class = "b", coef = "IageE2:EBV"),
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "sd", group = "ID"),
  prior(normal(0, 1), class = "sd", group = "Year"),
  prior(normal(0, 1), class = "sd", group = "cohort")
)

model3 <- brm(
  RS ~ density + age + I(age^2) + EBV + age*EBV + I(age^2)*EBV + (1 |  Year) + (1 |  ID) + (1 |  cohort),
  data = data,
  family = bernoulli(link = "logit"),
  prior = prior3,
  chains = 4,
  cores = 4,
  iter = 20000,
  warmup = 5000,
  control = list(adapt_delta = 0.95, max_treedepth = 10),
  silent = 0,
  save_pars = save_pars(all = TRUE)
)

conditional_effects(model3)
pp_check(model3)
summary(model3)
plot(model3)

loo_compare(loo(model3),loo(model(model(1))))
# model 4 EBV age class

prior4 <- c(
  prior(normal(-0.5, 0.5), class = "b", coef = "density"),
  prior(normal(0.5, 0.5), class = "b", coef = "age_classmadame"),
  prior(normal(0.5, 0.5), class = "b", coef = "age_classviellemadame"),
  prior(normal(0.5, 0.5), class = "b", coef = "EBV"),
  prior(normal(0.5, 1), class = "b", coef = "EBV:age_classmadame"),
  prior(normal(-0.5, 1), class = "b", coef = "EBV:age_classviellemadame"),
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "sd", group = "ID"),
  prior(normal(0, 1), class = "sd", group = "Year"),
  prior(normal(0, 1), class = "sd", group = "cohort")
)

model4 <- brm(
  RS ~ density  + EBV + EBV*age_class + (1 |  Year) + (1 |  ID) + (1 |  cohort),
  data = data,
  family = bernoulli(link = "logit"),
  prior = prior4,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 10),
  silent = 0,
  save_pars = save_pars(all = TRUE)
)

conditional_effects(model4)
pp_check(model4)
summary(model4)
plot(model4)


# model 5 bv et hlg

ggplot(data, aes(x = EBV, y = hlg12)) +
  geom_smooth(method = "lm") +
  geom_point() +
  labs(title = "cor btw EBV et hlg",
       x = "EBV",
       y = "hlg12") +
  theme_ipsum()

cor(data$EBV,data$hlg12,use = "complete.obs")


prior5 <- c(
  prior(normal(-0.5, 0.5), class = "b", coef = "density"),
  prior(normal(0.5, 0.5), class = "b", coef = "age"),
#  prior(normal(0.5, 0.5), class = "b", coef = "IageE2"),
  prior(normal(0.5, 0.5), class = "b", coef = "hlg12"),
  prior(normal(0.5, 0.5), class = "b", coef = "EBV"),
# prior(normal(0.5, 1), class = "b", coef = "age:EBV"),
#  prior(normal(-0.5, 1), class = "b", coef = "IageE2:EBV"),
#  prior(normal(0.5, 1), class = "b", coef = "age:hlg12"),
#  prior(normal(-0.5, 1), class = "b", coef = "IageE2:hlg12"),
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "sd", group = "ID"),
  prior(normal(0, 1), class = "sd", group = "Year"),
  prior(normal(0, 1), class = "sd", group = "cohort")
)

model5 <- brm(
  RS ~ density  + age + EBV + hlg12 + (1 |  Year) + (1 |  ID) + (1 |  cohort),
  data = data,
  family = bernoulli(link = "logit"),
  prior = prior5,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 2000,
  control = list(adapt_delta = 0.95, max_treedepth = 10),
  silent = 0,
  save_pars = save_pars(all = TRUE)
)

conditional_effects(model5)
pp_check(model5)
summary(model5)
plot(model5)


# looking for collinearity

loo_compare(loo(model3),loo(model(model1)))
# Extract posterior samples
posterior_samples <- as_draws_matrix(model5)

# Select only regression coefficients
posterior_betas <- posterior_samples[, c("b_age","b_density","b_EBV", "b_hlg12","b_age:EBV","b_age:hlg12","b_IageE2:hlg12","b_IageE2:EBV")]

# Compute correlation matrix
cor_matrix <- cor(posterior_betas)

# Print correlation matrix
print(cor_matrix)

library(bayesplot)

# Posterior scatter plot with correlation
mcmc_pairs(posterior_betas)



# animal model
ped <- ped %>% rename(
  ID = Id,
)

data <- data %>%
  left_join(ped, by = "ID")

data <- data %>%
  mutate(animal = ID)

library(nadiv)

ped <- ped[,c("ID","Dam","Sire")]

ped <- ped %>%
  rename(animal = ID)
  
  
Ped <- prepPed(ped)

Amat <- as.matrix(makeA(Ped))

data$age_class <-  factor(data$age_class,levels = c("tite jeune", "madame", "vielle madame"))

### adding LRS

data <- data %>%
  left_join(LRS %>% select(ID, LRS.wean),by = "ID")



bf_wtd <- bf(wt ~  1 + age_class*sq.date + (1 | a | gr(animal, cov = Amat)) + (1 | b |  Dam) + (1 | c | cohort) + (1 | d | ID) )
bf_hlg <- bf(hlg ~ 1 + age_class*sq.date + (1 | a | gr(animal, cov = Amat)) + (1 | b |  Dam) + (1 | c | cohort) + (1 | d | ID) )
bf_LRS <- bf(LRS.wean ~ 1 + (1 | a | gr(animal, cov = Amat)) + (1 | b |  Dam) + (1 | c | cohort) + (1 | d | ID) )



Model_test <- brm(
  bf_wtd,
  data = data,
  family = gaussian,
  prior = c(
    # Fixed Effects
    # prior(normal(0, 10), class = "b"),  # Fixed effect coefficients
    prior(student_t(3, 52.2, 8.4), class = "Intercept"),
  
    # coefficients for age class and julian days
    prior(student_t(3, 10, 8.4), class = "b", coef = "age_classviellemadame"),
    prior(student_t(3, 10, 8.4), class = "b", coef = "age_classmadame"),
    prior(student_t(3, 10, 8.4), class = "b", coef = "sq.date"),

    # interaction
    prior(student_t(3, 10, 8.4), class = "b", coef = "age_classmadame:sq.date"),
    prior(student_t(3, 10, 8.4), class = "b", coef = "age_classviellemadame:sq.date"),

    # Random Effects
    prior(student_t(3, 1, 8.4), class = "sd", group = "Dam"),  
    prior(student_t(3, 1, 8.4), class = "sd", group = "ID"),
    prior(student_t(3, 1, 8.4), class = "sd", group = "cohort"),
    prior(student_t(3, 1, 8.4), class = "sd", group = "animal"),
    
    # Residual Variance
    prior(student_t(3, 0, 8.4), class = "sigma")
),
data2 = list(Amat = Amat),
chains = 4, cores = 4, iter = 2000, warmup = 500,
control = list(adapt_delta = 0.95, max_treedepth = 15),
)

plot(Model_test)
summary(Model_test)
pp_check(Model_test)
conditional_effects(Model_test)

v_animal <- (VarCorr(Model_test, summary = FALSE)$animal$sd)^2
v_Dam <- (VarCorr(Model_test, summary = FALSE)$Dam$sd)^2
v_cohort <-  (VarCorr(Model_test, summary = FALSE)$cohort$sd)^2
v_ID <-  (VarCorr(Model_test, summary = FALSE)$ID$sd)^2
v_r <- (VarCorr(Model_test, summary = FALSE)$residual$sd)^2


h2_mass <- as.mcmc(v_animal / (v_animal + v_Dam + v_cohort + v_r))

plot(h2_mass)

summary(h2_mass)




hist(data$LRS.wean)





model1 <- brm(
  bf_wtd + bf_hlg + bf_LRS + set_rescor(TRUE),
  data = data,
   family = gaussian(),
  # priors_mv <- c(
  #   # Fixed Effects for wt
  #   prior(student_t(3, 61.3, 9.8), class = "Intercept", resp = "wt"),
  #   prior(student_t(3, 10, 8.4), class = "b", coef = "age_classmadame", resp = "wt"),
  #   prior(student_t(3, 10, 8.4), class = "b", coef = "age_classviellemadame", resp = "wt"),
  #   prior(student_t(3, 10, 8.4), class = "b", coef = "sq.date", resp = "wt"),
  #   prior(student_t(3, 10, 8.4), class = "b", coef = "age_classmadame:sq.date", resp = "wt"),
  #   prior(student_t(3, 10, 8.4), class = "b", coef = "age_classviellemadame:sq.date", resp = "wt"),
  #   
  #   # Fixed Effects for hlg
  #   prior(student_t(3, 23.2, 4.1), class = "Intercept", resp = "hlg"),
  #   prior(student_t(3, 10, 8.4), class = "b", coef = "age_classmadame", resp = "hlg"),
  #   prior(student_t(3, 10, 8.4), class = "b", coef = "age_classviellemadame", resp = "hlg"),
  #   prior(student_t(3, 10, 8.4), class = "b", coef = "sq.date", resp = "hlg"),
  #   prior(student_t(3, 10, 8.4), class = "b", coef = "age_classmadame:sq.date", resp = "hlg"),
  #   prior(student_t(3, 10, 8.4), class = "b", coef = "age_classviellemadame:sq.date", resp = "hlg"),
  #   
  #   # Random Effects for both responses
  #   prior(student_t(3, 0, 9.8), class = "sd", group = "animal", resp = "wt"),
  #   prior(student_t(3, 0, 4.1), class = "sd", group = "animal", resp = "hlg"),
  #   
  #   prior(student_t(3, 0, 9.8), class = "sd", group = "cohort", resp = "wt"),
  #   prior(student_t(3, 0, 4.1), class = "sd", group = "cohort", resp = "hlg"),
  #   
  #   prior(student_t(3, 0, 9.8), class = "sd", group = "Dam", resp = "wt"),
  #   prior(student_t(3, 0, 4.1), class = "sd", group = "Dam", resp = "hlg"),
  #   
  #   prior(student_t(3, 0, 9.8), class = "sd", group = "ID", resp = "wt"),
  #   prior(student_t(3, 0, 4.1), class = "sd", group = "ID", resp = "hlg"),
  #   
  #   # Residual Variance for Each Response
  #   prior(student_t(3, 0, 9.8), class = "sigma", resp = "wt"),
  #   prior(student_t(3, 0, 4.1), class = "sigma", resp = "hlg")
  #   
  # ),
  data2 = list(Amat = Amat),
  chains = 4, cores = 4, iter = 2000, warmup = 500,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
)


pp_check(model1,resp = "hlg")
pp_check(model1,resp = "wt")
pp_check(model1,resp = "LRS.wean")

summary(model1)
conditional_effects(model1)




print(ggplot(data = data1, aes(x = cohort, y = LRS)) + geom_smooth(method = "loess") + geom_point())

plot(data1$LRS ~ data1$cohort)

summary(model1)

v_animal <- (VarCorr(model1, summary = FALSE)$animal$sd)^2
v_ID <- (VarCorr(model1, summary = FALSE)$ID$sd)^2
v_Dam <- (VarCorr(model1, summary = FALSE)$Dam$sd)^2
v_cohort <-  (VarCorr(model1, summary = FALSE)$cohort$sd)^2
v_r <- (VarCorr(model1, summary = FALSE)$residual$sd)^2

h2_hlg <- as.mcmc(v_animal[, "hlg_Intercept"] / (v_animal[,"hlg_Intercept"] +  v_ID[,"hlg_Intercept"] + v_Dam[,"hlg_Intercept"] + v_cohort[,"hlg_Intercept"] +  v_r[,"hlg"]))
h2_mass <- as.mcmc(v_animal[, "wt_Intercept"] / (v_animal[, "wt_Intercept"] + v_ID[,"wt_Intercept"] + v_Dam[, "wt_Intercept"] + v_cohort[, "wt_Intercept"] + v_r[, "wt"]))

plot(h2_hlg)
plot(h2_mass)

summary(h2_hlg)
summary(h2_mass)

cor_g <- as.mcmc((VarCorr(model1, summary = FALSE)$animal$cor[, 1, 2]))

plot(cor_g)

summary(cor_g)
