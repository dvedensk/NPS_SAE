library(readr)
library(dplyr)

households <- read_csv(
    file = file.path("data", "psam_h06.csv"),
    col_select = c("SERIALNO", "ACCESSINET", "TEL")
)
persons <- read_csv(
    file = file.path("data", "psam_p06.csv"),
    col_select = c(
        "SERIALNO",
        "PWGTP",
        "PUMA",
        # "STATE", # not necessary while we consider just CA
        "AGEP",
        "CIT",
        "RAC1P", # Race (9 levels)
        "HISP", # Hispanic origin, "01" corresponds to non-hispanic
        "SEX",
        "WAGP", # Wages/salary past 12 months
        "SCHL", # Collapse into fewer categories
        "HICOV", # Health insurance (1 = covered, 2 = not covered)
        "PUBCOV", # Public health coverage (1 = no, 2 = yes)
        "POVPIP", # Income-to-poverty ratio (0-501)
        "SSP", # Social Security income
        "PERNP"
    )
) # Total earnings

acs_pop <- persons %>% left_join(households, by = c("SERIALNO"))
# exclude group quarters:
acs_pop <- acs_pop %>%
    filter(substr(SERIALNO, start = 5, stop = 6) == "HU") %>%
    filter(!is.na(WAGP))

# recode response so 0 is "No" and 1 is "Yes" and bin covariates
# Keep AGEP, POVPIP, SSP, PERNP as continuous (unbinned) for use in sampling weights
acs_pop_agep_continuous <- acs_pop$AGEP # Save continuous version before binning
acs_pop <- acs_pop %>%
    mutate(
        HICOV = ifelse(HICOV == 2, 0, 1),
        PUBCOV = ifelse(PUBCOV == 2, 1, 0)
    ) %>%
    mutate(
        SEX = ifelse(SEX == 1, "Male", "Female"),
        AGEP_binned = cut(AGEP, c(0, 18, 34, 65, Inf), right = F),
        RAC1P = case_when(RAC1P == 1 ~ "White",
            RAC1P == 2 ~ "Black",
            RAC1P == 6 ~ "Asian",
            .default = "Other"
        )
    ) %>%
    mutate(RAC1P = ifelse(HISP == "01", RAC1P, "Hispanic")) %>%
    mutate(
        SEX = as.factor(SEX),
        RAC1P = as.factor(RAC1P),
        PUMA = as.factor(PUMA),
        # Replace missing values with 0 for continuous variables used in sampling
        POVPIP = ifelse(is.na(POVPIP), 0, POVPIP),
        SSP = ifelse(is.na(SSP), 0, SSP),
        PERNP = ifelse(is.na(PERNP), 0, PERNP)
    ) %>%
    select(-HISP)

write_csv(acs_pop, file = file.path("data", "ACS_NPS_pop.csv"))
