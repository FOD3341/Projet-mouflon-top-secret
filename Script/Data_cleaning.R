

### Loading packages

if (!require("tidyverse")) install.packages("tidyverse") 
if (!require("here")) install.packages("here")


### Read.data
df <- read.csv("Raw_data/All_ewes_23_corrected.csv",header = T,sep =";")
Repro <- read.csv("Raw_data/RamMtn_Repro_71-23.csv",header = T,sep = ";")
EBV <- read.csv("Raw_data/ebv_hlF.csv",header = T,sep = ",",dec =".")

### data cleaning

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

df5[,c("age","density","hlg","mu")] <- scale(df5[,c("age","density","hlg","mu")])

# save the data frame

write.csv(df5, "Cleaned_data/cleaned_data.csv", row.names = FALSE)
