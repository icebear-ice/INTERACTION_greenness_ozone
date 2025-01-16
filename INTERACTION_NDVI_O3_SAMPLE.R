# Required packages
library(dlnm) 
library(gnm) 
library(splines) 
library(dplyr)
library(data.table)
library(tsModel)

# Data import options: (choose one option)

# Option 1: Read your own data
data_import <- function(file_path) {
  # Read CSV file
  data <- read.csv(file_path, encoding = "UTF-8")
  
  #' Required columns in your data:
  #' - ALL: Daily all-cause mortality counts
  #' - O3: Daily ozone concentration
  #' - T: Daily mean temperature
  #' - RH: Daily relative humidity
  #' - NDVI: Greenness index
  #' - year: Year
  #' - month: Month
  #' - dow: Day of week
  #' - code: Location code
  
  # Convert to data.table
  return(data.table(data))
}

# Option 2: Generate random sample data
data_generate <- function(n = 1000) {
  set.seed(123)
  data <- data.frame(
    date = seq.Date(as.Date("2013-01-01"), by = "day", length.out = n),
    ALL = rpois(n, lambda = 10),  # Daily mortality
    O3 = pmax(rnorm(n, mean = 40, sd = 10), 0),  # O3 concentration
    T = rnorm(n, mean = 20, sd = 5),  # Temperature
    RH = pmin(pmax(rnorm(n, mean = 70, sd = 10), 0), 100),  # Relative humidity
    NDVI = pmin(pmax(rnorm(n, mean = 0.5, sd = 0.1), 0), 1),  # NDVI
    year = format(seq.Date(as.Date("2013-01-01"), by = "day", length.out = n), "%Y"),
    month = format(seq.Date(as.Date("2013-01-01"), by = "day", length.out = n), "%m"),
    dow = format(seq.Date(as.Date("2013-01-01"), by = "day", length.out = n), "%w"),
    code = rep("001", n)
  )
  return(data.table(data))
}

# Usage example:
if(FALSE) {  # Set to TRUE to run
  # Either import your data:
  # data <- data_import("path/to/your/data.csv")
  
  # Or use generated sample data:
  data <- data_generate(1000)
  
  # Check data structure
  str(data)
  head(data)
}


# Data preprocessing
data_processed <- data %>%
  filter(O3 > 0) %>%
  mutate(
    # Calculate moving averages
    lag0 = O3,
    lag01 = as.numeric(runmean(O3, 2, align = "right")),
    lag02 = as.numeric(runmean(O3, 3, align = "right")),
    lag03 = as.numeric(runmean(O3, 4, align = "right")),
    lag04 = as.numeric(runmean(O3, 5, align = "right")),
    lag05 = as.numeric(runmean(O3, 6, align = "right")),
    # Temperature moving average
    TEMP = runmean(T, 4, align = "right")
  )

# Create stratification
data_processed <- data.table(data_processed)
data_processed[, stratum := factor(paste(code, year, month, dow, sep = ":"))]
data_processed[, keep := sum(ALL) > 0, by = stratum]

# Define analysis parameters
PMs <- c("lag0", "lag01", "lag02", "lag03", "lag04", "lag05")
IQR <- diff(quantile(data_processed$O3, c(0.25, 0.75), na.rm = TRUE))

# Result storage matrix
result <- matrix(NA, nrow = length(PMs), ncol = 10)
colnames(result) <- c("lag", "RR", "RR_low", "RR_high", 
                      "RR_low_green", "RR_low_green_ci_low", "RR_low_green_ci_high",
                      "RR_high_green", "RR_high_green_ci_low", "RR_high_green_ci_high")

# Main analysis loop
for(i in seq_along(PMs)) {
  # Main model
  modfull <- gnm(ALL ~ data_processed[[PMs[i]]] + 
                   ns(TEMP, 6) + ns(RH, 3),
                 eliminate = stratum,
                 data = data_processed,
                 family = quasipoisson,
                 subset = keep)
  
  # Interaction analysis
  intval <- quantile(data_processed$NDVI, c(0.05, 0.95))
  cbint1 <- data_processed[[PMs[i]]] * (data_processed$NDVI - intval[1])
  cbint2 <- data_processed[[PMs[i]]] * (data_processed$NDVI - intval[2])
  
  # Update models with interactions
  modint1 <- update(modfull, . ~ . + cbint1)
  modint2 <- update(modfull, . ~ . + cbint2)
  
  # Calculate relative risks
  b0 <- summary(modfull)$coefficients[1, 1]
  se0 <- summary(modfull)$coefficients[1, 2]
  b1 <- summary(modint1)$coefficients[1, 1]
  se1 <- summary(modint1)$coefficients[1, 2]
  b2 <- summary(modint2)$coefficients[1, 1]
  se2 <- summary(modint2)$coefficients[1, 2]
  
  # Store results
  result[i,] <- c(
    i,
    exp(b0 * IQR),
    exp((b0 - 1.96 * se0) * IQR),
    exp((b0 + 1.96 * se0) * IQR),
    exp(b1 * IQR),
    exp((b1 - 1.96 * se1) * IQR),
    exp((b1 + 1.96 * se1) * IQR),
    exp(b2 * IQR),
    exp((b2 - 1.96 * se2) * IQR),
    exp((b2 + 1.96 * se2) * IQR)
  )
}

# Save results
write.csv(result, "analysis_results.csv", row.names = FALSE)