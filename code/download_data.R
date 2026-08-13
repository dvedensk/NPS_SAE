dir.create("data", recursive = TRUE, showWarnings = FALSE)

urls <- c(
  people = "https://www2.census.gov/programs-surveys/acs/data/pums/2023/5-Year/csv_pca.zip",
  households = "https://www2.census.gov/programs-surveys/acs/data/pums/2023/5-Year/csv_hca.zip"
)

zip_paths <- setNames(
  file.path("data", basename(unname(urls))),
  names(urls)
)

for (nm in names(urls)) {
  message("Downloading ", nm, " file...")

  download.file(
    url = urls[[nm]],
    destfile = zip_paths[[nm]],
    mode = "wb",
    quiet = FALSE
  )

  message("Unzipping ", basename(zip_paths[[nm]]), " into data/ ...")

  unzip(
    zipfile = zip_paths[[nm]],
    exdir = "data",
    overwrite = TRUE
  )
}
