# Single source of truth for the demo: which counties, which dates, which files.
# Sourced by every other script in this repo.

# 10 single-FIPS counties with dense pre-2023 NWSS sample coverage, mixing large
# metros and smaller counties, spread across regions.
FIPS_SUBSET <- c(
  "17031", # Cook County, IL
  "48201", # Harris County, TX
  "04013", # Maricopa County, AZ
  "36071", # Orange County, NY
  "26099", # Wayne County, MI
  "36029", # Erie County, NY
  "08069", # Larimer County, CO
  "49049", # Salt Lake County, UT
  "55025", # Dane County, WI
  "04027"  # Pima County, AZ
)

FIPS_COUNTY_NAMES <- c(
  "17031" = "Cook County, IL",
  "48201" = "Harris County, TX",
  "04013" = "Maricopa County, AZ",
  "36071" = "Orange County, NY",
  "26099" = "Wayne County, MI",
  "36029" = "Erie County, NY",
  "08069" = "Larimer County, CO",
  "49049" = "Salt Lake County, UT",
  "55025" = "Dane County, WI",
  "04027" = "Pima County, AZ"
)

# Demo analysis window: restricted to JHU's active reporting range so cases and
# wastewater data overlap. JHU stopped updating in March 2023.
START_DATE <- as.Date("2020-01-01")
END_DATE   <- as.Date("2023-03-31")

# Max consecutive missing weeks to interpolate across (longer gaps are left NA).
INTERP_THRESHOLD <- 3

# File paths
NWSS_RAW_FILE  <- "data/raw/nwss_county_subset.csv"
JHU_RAW_FILE   <- "data/raw/jhu_cases_subset.csv"
WW_CLEANED_FILE    <- "data/processed/ww_cleaned.csv"
CASES_CLEANED_FILE <- "data/processed/cases_cleaned.csv"

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
