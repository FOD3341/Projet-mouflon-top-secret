

### Loading packages

if (!require("tidyverse")) install.packages("tidyverse")
if (!require("here")) install.packages("here")
if (!require("ggplot2")) install.packages("ggplot2") 

### Read.data
df <- read.csv("Raw_data/All_ewes_23_corrected.csv",header = T,sep =";")
Repro <- read.csv("Raw_data/RamMtn_Repro_71-23.csv",header = T,sep = ";")
EBV <- read.csv("Raw_data/ebv_hlF.csv",header = T,sep = ",",dec =".")

### data cleaning

# mean and sd of the standard distribution of the EBV


# Define the threshold for a "substantial" decline (e.g., 10% drop)
threshold <- 0.1

# Create a function to detect breaks
detect_breaks <- function(df) {
  df %>%
    arrange(ID, Age,sq.date) %>%  # Sort by individual ID and age
    group_by(ID) %>%
    mutate(
      # Calculate the difference between consecutive years
      horn_diff = hlg - lag(hlg),  # Direct difference between consecutive years
      # Calculate the percentage change
      horn_diff_pct = (hlg - lag(hlg)) / lag(hlg)
    ) %>%
    filter(!is.na(horn_diff_pct) & abs(horn_diff_pct) >= threshold)   # Filter by substantial decline
}

# Apply the function to your dataset
df_with_breaks <- detect_breaks(df)

df_with_breaks_filter <- df_with_breaks %>%
  filter(horn_diff > 2.5 | horn_diff < -2.5) %>%
  distinct(ID)

### plot the subpect growth patern

# vector of the error ID growth
names_ID <- c("W17","U19","Q8","Q3","P13","O12","N4","M7","M26","28W","2N","2S","30V","37U","3M","3X","40U","58T","7I","9T","A50","E12","E12","J5","L6")

plots <- list()

for(i in names_ID){
  df6 <- df %>%
    filter(ID == i) %>%
    arrange(Year,sq.date) %>%
    mutate(rank = row_number())
  
  # Create plot and store it in list
  plots[[i]] <- ggplot(df6, aes(x = rank, y = hlg)) +
    geom_smooth(method = "loess") +
    geom_point() +
    ggtitle(paste("ID:", i))  # Label plot by ID
}

hlg_max_min <- df %>%
  group_by(ID) %>%
  summarise(max_hlg = max(hlg,na.rm = TRUE), min_hlg = min(hlg,na.rm = TRUE))

# removing data points because of horn breaking and measurment error
df_clean <- df %>%
  left_join(hlg_max_min, by = "ID") %>%
  filter(
    !(ID == "L6" & hlg > 27.5) &
      !(ID == "J5" & hlg < 20 & Age == 10) &
      !(ID == "E12" & hlg < 12) &
      !(ID == "A50" & Age >= 6) &  # Fixed column name
      !(ID == "9T" & Age > 6) &
      !(ID == "40U" & hlg < 17 & Age == 7) &
      !(ID == "3X" & Age > 5) &
      !(ID == "2N" & Age < 3) &
      !(ID == "M26" & hlg < 15 & Age == 5) &
      !(ID == "M7" & ((hlg == 19.1 & Age == 4) | (hlg == 17.1 & Age == 11))) &
      !(ID == "Q3" & Age > 3) &
      !(ID %in% c("Q8", "U19")) &  # Multiple IDs removed at once
      !(ID == "7I" & hlg == min_hlg) &
      !(ID == "58T" & hlg == max_hlg) &
      !(ID == "37U" & hlg == max_hlg) &
      !(ID == "30V" & hlg == max_hlg) &
      !(ID == "2S" & hlg == max_hlg) &
      !(ID == "28W" & hlg == max_hlg) &
      !(ID == "N4" & hlg == max_hlg) &
      !(ID == "O12" & hlg == min_hlg) &
      !(ID == "P13" & hlg == min_hlg) &
      !(ID == "W17" & hlg == max_hlg) &
      !(ID == "3X" & hlg == max_hlg) &
      !(ID == "18R")
  )


ggplot(df_clean, aes(x = Age, y = hlg, group = ID)) +
  geom_line() +
  geom_point()



### merge data.frame
df1 <- df_clean %>%
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

# save the data frame

write.csv(df5, "Cleaned_data/cleaned_data.csv", row.names = FALSE)
