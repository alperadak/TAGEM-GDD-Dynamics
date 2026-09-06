#=========================================================
# 1. CLEAR WORKSPACE
#=========================================================
rm(list = ls())
gc()


#=========================================================
# 2. LOAD REQUIRED LIBRARIES
#=========================================================
library(dplyr)
library(envirotypeR)
library(readxl)


#=========================================================
# 3. DEFINE FILE PATHS
#=========================================================
coord_file <- "E:/ALPER_BURAK_ENV/Institues_with_coords.xlsx"
output_file <- "E:/ALPER_BURAK_ENV/env_data_1981_2025.rds"


#=========================================================
# 4. LOAD ENVIRONMENT LOCATIONS
#=========================================================
env_locations <- read_excel(coord_file)

# Check number of unique environments
length(unique(env_locations$Env))

# Inspect data
head(env_locations)


#=========================================================
# 5. PREPARE INPUT DATA FOR envirotypeR
#=========================================================
df <- env_locations %>%
  rename(
    env = Env,
    lat = Latitude,
    lon = Longitude
  ) %>%
  mutate(
    start = as.Date("1981-01-01"),
    end   = as.Date("2025-12-31")
  )

# Check prepared input
head(df)


#=========================================================
# 6. DOWNLOAD DAILY WEATHER DATA
#=========================================================
env.data <- envirotypeR::get_weather(
  env.id    = df$env,
  lat       = df$lat,
  lon       = df$lon,
  start.day = df$start,
  end.day   = df$end,
  parallel  = TRUE
)


#=========================================================
# 7. CHECK DOWNLOADED DATA
#=========================================================
head(env.data)

# Environments successfully downloaded
unique(env.data$env)

# Number of environments
length(unique(env.data$env))


#=========================================================
# 8. SAVE WEATHER DATA
#=========================================================
saveRDS(
  env.data,
  file = output_file
)


#=========================================================
# 9. CONFIRM SAVED FILE
#=========================================================
file.exists(output_file)
