

### Loading packages

if (!require("tidyverse")) install.packages("tidyverse")
if (!require("here")) install.packages("here")
if (!require("ggplot2")) install.packages("ggplot2") 

### Read.data
df <- read.csv("Raw_data/All_ewes_23_corrected.csv",header = T,sep =";")
Repro <- read.csv("Raw_data/RamMtn_Repro_71-23.csv",header = T,sep = ";")
EBV <- read.csv("Raw_data/ebv_hlF.csv",header = T,sep = ",",dec =".")
surv <- read.csv("Raw_data/surv.csv",header = T,sep = ",",dec =".")
hornLength <- read.csv("Raw_data/hornLength.csv",header = T,sep = ",",dec =".")

### data cleaning

### merge data.frame
df_horn <- hornLength %>%
  left_join(df_ewes, by = c("ID", "Annee")) 

plot(df_horn$Lhorn, df_horn$Left.Horn.Lenght)

### merge data.frame
df1 <- df %>%
  left_join(Repro, by = c("ID", "Year"))

### calculate a density for each year

df2 <- df1 %>%
  group_by(Year) %>%  # Group data by year
  mutate(density = n_distinct(ID))

### Remove Ewes wich reproductive success is 8 or 9

df2 <- subset(df2,Reproduction < 8)

### Create a new column with the reproductive as a binary variable

df3 <- df2 %>%
  mutate(RS = if_else(Reproduction == 5, 1, 0))

### Group by year and ID
df4 <- df3 %>%
  group_by(ID, Year) %>%  # Grouper par ID et année
  summarise( 
    age = first(Age),
    density = first(density),
    RS = first(RS),
    sex = first(Lamb_Sex),
    hlg = max(hlg),
    wtd12 = first(wtd12),
    wtd114 = first(wtd114)
  )

# Left join EBV

EBV <- EBV %>%
  rename(ID = X)

df5 <- df4 %>%
  left_join(EBV,by = "ID")

# mean and sd of the standard distribution of the EBV

for(i in 1:nrow(df5)){
  df5[i,"mu"] <- mean(as.numeric(df5[i,9:(ncol(df5)-2)]))
  df5[i,"sd"] <- sd(as.numeric(df5[i,9:(ncol(df5)-2)]))
}

df5 <- df5 %>% select(-starts_with("V"))

# adding cohort
df5[,"cohort"] <- df5[,"Year"] - df5["age"]

# Create age classes based on fixed cutoffs
df5$age_class <- cut(df5$age,
                     breaks = c(2, 3, 8, max(df5$age)),  # Define the break points
                     labels = c("tite jeune", "madame","vielle madame"), 
                     right = FALSE)  # right = FALSE means the interval is left-closed
# Scaling the data

df5[,c("age_scaled","density_scaled","hlg_scaled","mu_scaled")] <- scale(df5[,c("age","density","hlg","mu")])


# looking for error in the data base

# horn length
ggplot(df5, aes(x = hlg)) +
  geom_histogram(binwidth = 0.5)

# age
ggplot(df5, aes(x = age)) +
  geom_histogram(binwidth = 0.5)

# age_class
ggplot(df5, aes(x = age_class)) +
  geom_bar()
  
# mu
ggplot(df5, aes(x = mu)) +
  geom_histogram(binwidth = 0.5)

# horn by age class
ggplot(df5, aes(x = age,y = hlg,group = ID)) +
  geom_line() +
  geom_point()

# removing Q3, because of horn breaking
df5 <- df5 %>%
  filter(!(ID == "Q3" & age > 3))

library(dplyr)

library(dplyr)

# Define the threshold for a "substantial" decline (e.g., 10% drop)
threshold <- 0.1

# Create a function to detect breaks
detect_breaks <- function(df) {
  df %>%
    arrange(ID, age) %>%  # Sort by individual ID and age
    group_by(ID) %>%
    mutate(
      # Calculate the difference between consecutive years
      horn_diff = hlg - lag(hlg),  # Direct difference between consecutive years
      # Calculate the percentage change
      horn_diff_pct = (hlg - lag(hlg)) / lag(hlg)
    ) %>%
    filter(!is.na(horn_diff_pct) & abs(horn_diff_pct) >= threshold)  # Filter by substantial decline
}

# Apply the function to your dataset
df_with_breaks <- detect_breaks(df5)

# View the individuals where a substantial decline in horn size was detected
print(df_with_breaks)


# Apply the function to your dataset
df_with_breaks <- detect_breaks(df)

# save the data frame

write.csv(df5, "Cleaned_data/cleaned_data.csv", row.names = FALSE)
