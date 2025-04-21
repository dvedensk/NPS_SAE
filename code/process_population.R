require(readr)

households <- read_csv("psam_h06.csv", col_select=c("SERIALNO", "ACCESSINET", "TEL"))
persons <- read_csv("psam_p06.csv", col_select=c("SERIALNO",
                                                 "PWGTP",
                                                 "PUMA",
                                                 "STATE",
                                                 "AGEP",
                                                 "RAC1P", #Race (9 levels)
                                                 "SEX",
                                                 "WAGP", #Wages/salary past 12 months
                                                 "SCHL", #Collapse into fewer categories
                                                 "HICOV")) #Health insurance (1 = covered, 2 = not covered)

acs.pop <- persons %>% left_join(households, by=c("SERIALNO"))
acs.pop <- acs.pop %>% filter(substr(SERIALNO, start=5,stop=6)=="HU") %>% #exclude group quarters
                       filter(!is.na(WAGP))

acs.pop <- acs.pop <- acs.pop %>% mutate(HICOV = ifelse(HICOV==2, 0, 1)) #recode so 0 = "NO"

write_csv(acs.pop, file="ACS_NPS_pop.csv")
