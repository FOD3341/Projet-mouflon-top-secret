# read packages

if (!require("brms")) install.packages("brms")

if (!require("ggplot2")) install.packages("ggplot2")

if (!require("tidyverse")) install.packages("tidyverse") 

if (!require("tidybayes")) install.packages("tidybayes")

if (!require("hrbrthemes")) install.packages("hrbrthemes")

# read data

data <- read.csv("Cleaned_data/cleaned_data.csv", header = T, sep = ",")
horn_data <- read.csv("Cleaned_data/adjhlg.csv", header = T, sep = ",")

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
    cohort = first(cohort)
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
  prior(normal(-0.5, 0.5), class = "b", coef = "IageE2"),
  prior(normal(0.5, 0.5), class = "b", coef = "hlg12"),
  prior(normal(0.5, 1), class = "b", coef = "age:hlg12"),
  prior(normal(-0.5, 1), class = "b", coef = "IageE2:hlg12"),
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "sd", group = "ID"),
  prior(normal(0, 1), class = "sd", group = "Year"),
  prior(normal(0, 1), class = "sd", group = "cohort")
)

model <- brm(
  RS ~ density + age + I(age^2) + hlg12 + hlg12*age + I(age^2)*hlg12 + (1 |  Year) + (1 |  ID) + (1 |  cohort) ,
  data = data,
  family = bernoulli(link = "logit"),
  prior = prior,
  chains = 4,
  cores = 4,
  iter = 20000,
  warmup = 5000,
  control = list(adapt_delta = 0.95, max_treedepth = 10),
  silent = 0,
  save_pars = save_pars(all = TRUE)
)

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

ggplot(data, aes(x = EBV, y = hlg)) +
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
  prior(normal(0.5, 0.5), class = "b", coef = "IageE2"),
  prior(normal(0.5, 0.5), class = "b", coef = "hlg12"),
  prior(normal(0.5, 0.5), class = "b", coef = "EBV"),
  prior(normal(0.5, 1), class = "b", coef = "age:EBV"),
  prior(normal(-0.5, 1), class = "b", coef = "IageE2:EBV"),
  prior(normal(0.5, 1), class = "b", coef = "age:hlg12"),
  prior(normal(-0.5, 1), class = "b", coef = "IageE2:hlg12"),
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "sd", group = "ID"),
  prior(normal(0, 1), class = "sd", group = "Year"),
  prior(normal(0, 1), class = "sd", group = "cohort")
)

model5 <- brm(
  RS ~ density  + age + I(age^2) + EBV + hlg12 + age*EBV + I(age^2)*EBV + age*hlg12 + I(age^2)*hlg12 + (1 |  Year) + (1 |  ID) + (1 |  cohort),
  data = data,
  family = bernoulli(link = "logit"),
  prior = prior5,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 1000,
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