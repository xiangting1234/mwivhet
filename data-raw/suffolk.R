## code to prepare `suffolk` dataset goes here
# 1. Load haven to read Stata files
library(haven)

# 2. Read your .dta file from your laptop
# (Adjust the path to where your .dta file is currently located)
suffolk <- read_dta("suffolk_est_twoyears_demos.dta")

# 3. (Optional) Clean or format the data here
# e.g., my_demo_data <- as.data.frame(my_demo_data)

# 4. Save it into the package's data/ folder
usethis::use_data(suffolk, overwrite = TRUE)
