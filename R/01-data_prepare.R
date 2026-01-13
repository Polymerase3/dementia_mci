## ----------------------------READING PACKAGES---------------------------------
# setting random generator seed
set.seed(16748991)
seed_before <- .Random.seed

# creating vector of necessary packages
packages <- c(
  "dplyr",
  "showtext",
  "stringr",
  "data.table",
  "DBI",
  "duckdb"
)

## load development version of phiper
library(phiper)

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())

if (any(installed_packages == FALSE)) {
  pak::pkg_install(packages[!installed_packages])
}

# packages loading
invisible(lapply(packages, library, character.only = TRUE))

# add font + phiper use font in all plots
font_add_google("Montserrat", "monte")
phip_use_montserrat()
showtext_auto()

# removing unnecessary variables
rm(list = c('installed_packages', 'packages'))

## ---------------------LOADING DATA--------------------------------------------

# listing filenames of all files
csv_files <- list.files(path = "data_original",
                        pattern = "\\.csv$",
                        full.names = TRUE)

# removing unnecessary infos (pre- + suffix)
dataframe_names <- csv_files %>%
  str_replace_all(c(".csv" = "", "data_original/" = ""))

## importing 800 files with the data
data <- list()
for(x in unique(csv_files)){
  data[[x]] <- fread(x)
}

## adding sampleID
data <- lapply(seq_along(data), function(i) {
  data[[i]]$sampleID <- dataframe_names[i]
  data[[i]]
})

## combining all files into one dataset
## CAN and probably will take some time --> completely normal
dementia <- Reduce(rbind, data)

## cleaning colnames
colnames(dementia)[1] <- "peptideID"

## reordering the columns
dementia <- dementia[, c("sampleID", "peptideID", "input", "count", "fold_change")]

dementia <- as.data.frame(dementia)

# read the metadata
metadata <- read.csv("data_original/other_data/metadata_complete.csv")

# add the group control_or_MCI
metadata <- metadata %>%
  mutate(
    control_or_MCI = if_else(
      coalesce(control, 0L) == 1L | coalesce(MCI, 0L) == 1L,
      1L, 0L
    )
  )

# join to the dementia on sampleID
dementia <- dementia %>%
  left_join(
    metadata,
    by = c("sampleID" = "sample_id")
  )

## ---------------------------------- DATA EXPORT ------------------------------
## exporting the data to a .parquet file --> for phiper use
dementia_export <- dementia

## safety-check
dementia_export <- dementia_export %>%
  distinct(sampleID, peptideID, .keep_all = TRUE)

# 1. open (or create) a DuckDB database ----------------------------------------
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:", read_only = FALSE)

# 2. expose the data frame to DuckDB -------------------------------------------
duckdb::duckdb_register(con, "df_tbl", dementia_export)

# -- df is now a virtual table called "df_tbl" inside the DB
# 3. Persist it to Parquet -----------------------------------------------------
dbExecute(
  con,
  "COPY df_tbl TO 'data_original/other_data/dementia_full.parquet' (FORMAT PARQUET);"
)

# 4. Clean up ------------------------------------------------------------------
dbDisconnect(con, shutdown = TRUE)
rm(list = c("dementia", "dementia_export", "con"))
gc()
