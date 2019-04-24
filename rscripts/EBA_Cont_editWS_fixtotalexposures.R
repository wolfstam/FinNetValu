### EBA data set
# Info for the csv file: From the EBA stress test 2011 results website
# "https://eba.europa.eu/risk-analysis-and-data/eu-wide-stress-testing/2011/results"
# download the stress test tool for individual results, select all banks, and
# save the DPCache file. 

library(tidyverse)
library(readxl)
library(lubridate)

## Create portfolio matrices as in Cont paper
cont.data.raw <- read_csv("/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/2011 EU-wide stress test disclosure templates - revised DPCache.csv",
                          guess_max=2000)

# cont.data.parsed <- cont.data.raw %>%
#   mutate(date = ymd(DATE), exposure = parse_double(VALUE_NUM, locale=locale(decimal_mark=","))) %>%
#   mutate(bank_id = BANKING_GROUP_CODE, code = INFORMATION_CODE, country = C_COUNTRY_CODE)
cont.data.parsed <- cont.data.raw %>%
  mutate(date = ymd(DATE), exposure = VALUE_NUM) %>%
  mutate(bank_id = BANKING_GROUP_CODE, code = INFORMATION_CODE, country = C_COUNTRY_CODE)

ill.asset.codes <- c(residential.mortgage = 33013,
                     commercial.real.estate  = 33018,
                     retail.revoling = 33015,
                     retail.SME = 33016,
                     retail.other = 33017,
                     indirect.sovereign = 34017,
                     defaulted = 33020)
liq.asset.codes <- c(inst.client = 33010,
                     corporate = 33011,
                     sovereign = 34013,
                     sovereign = 34014,
                     sovereign = 34015,
                     direct.sovereign.derivatives = 34016)
liab.codes <- c(tier.1 = 30014)
total.exposure.code <- c(total.exposures = 33021)

eu.country.codes <- c(
    ## EU member states
    "AT", "BE", "BG", "CY", "CZ", "DK", "EE", "FI", "FR", "DE", "GR", "HU", "IE", "IT", "LV", "LT", "LU", "MT", "NL", "PL", "PT", "RO", "SK", "SI", "ES", "SE", "GB")
country.codes <- c(
    eu.country.codes, # EU member states
    "US", # United States
    "NO", # Norway
    "IS", # Iceland
    "LI", # Lichtenstein
    "JP", # Japan
    "A1", # Asia
    "E3", # Other non-EEA non-emerging countries
    "E5", # Eastern Europe non-EEA
    "M1", # Middle and South America
    "R5") # Restof the world

cont.data <- cont.data.parsed %>%
    filter(country %in% country.codes) %>%
    filter(code %in% c(ill.asset.codes, liq.asset.codes, liab.codes, total.exposure.code))

## Sanity check for maturities
cont.data %>%
    group_by(mat = MATURITY_CODE == 999, code, country, bank_id) %>%
    summarize(total = sum(exposure, na.rm=TRUE)) %>%
    ungroup() %>%
    arrange(bank_id, country)
    #spread(mat, total)

key <- function (nv, v) {
    ## Get name corresponding to value v in named vector nv
    sapply(v, function (x) names(nv)[nv == x])
}

cont.Pi <- cont.data %>%
    ## Get liquid assets of all maturities
    filter(is.na(MATURITY)) %>% ## Either total or no maturity given
    filter(code %in% liq.asset.codes) %>%
    ## combine sovereign exposures
    mutate(code_name = key(liq.asset.codes, code)) %>%
    group_by(bank_id, code_name, country) %>%
    summarize(exposure = 1e6*sum(exposure, na.rm=TRUE)) %>%
    ungroup() %>%
    ## spread asset-country combinations
    unite(asset, code_name, country) %>%
    spread(asset, exposure, fill=0)

Pi <- cont.Pi %>% select(-bank_id) %>% as.matrix()
rownames(Pi) <- cont.Pi$bank_id
colnames(Pi) <- colnames(cont.Pi %>% select(-bank_id))

aggregate.ill.assets <- function (asset_names) {
    sapply(asset_names, function (an) {
        if (grepl("residential.mortgage", an) |
            grepl("commercial.real.estate", an))
            an
        else
            "ill.asset"
    })
}

# compute difference between sum of individual exposures and total exposures in balance sheet
# add difference as additional category
test.totalexp <- cont.data %>%
  filter(code %in% c(33021))
test <- cont.data %>% 
  filter(code %in% c(33010, 33011, 33013, 33015, 33016, 33017, 33018, 33020)) %>% 
  group_by(bank_id, country) %>%
  mutate(cumsum = cumsum(exposure)) %>%
  mutate(sum.ind.exposures = max(cumsum))

# cont.Theta <- cont.data %>%
#     ## Get illiquid assets of all maturities
#     filter(is.na(MATURITY)) %>% ## Either total or no maturity given
#     filter(code %in% ill.asset.codes) %>%
#     ## Name asset classes-country combinations
#     mutate(code_name = key(ill.asset.codes, code)) %>%
#     unite(asset, code_name, country) %>%
#     ## Aggregate as in Cont paper
#     mutate(asset = aggregate.ill.assets(asset)) %>%
#     group_by(bank_id, asset) %>%
#     summarize(exposure = 1e6*sum(exposure, na.rm=TRUE)) %>%
#     ungroup() %>%
#     spread(asset, exposure, fill=0)


cont.Theta <- cont.data %>%
  ## Get illiquid assets of all maturities
  filter(is.na(MATURITY)) %>% ## Either total or no maturity given
  filter(code %in% c(ill.asset.codes, total.exposure.code)) %>%
  ## Name asset classes-country combinations
  mutate(code_name = key(c(ill.asset.codes, total.exposure.code), code)) %>%
  unite(asset, code_name, country) %>%
  ## Aggregate as in Cont paper
  mutate(asset = aggregate.ill.assets(asset)) %>%
  group_by(bank_id, asset) %>%
  summarize(exposure = 1e6*sum(exposure, na.rm=TRUE)) %>%
  ungroup() %>%
  spread(asset, exposure, fill=0)


cont.Capital <- cont.data.parsed %>%
    filter(code %in% liab.codes) %>%
    filter(is.na(SCENARIO)) %>%
    select(bank_id, capital = exposure) %>%
    mutate(capital = 1e6*capital)

cap <- cont.Capital$capital
names(cap) <- cont.Capital$bank_id

write_csv(cont.Pi, "/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/Cont_Pi.csv")
write_csv(cont.Theta, "/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/Cont_Theta.csv")
write_csv(cont.Capital, "/Users/wolfgang/Desktop/Uni/HIWI/Bertschinger/SysRisk/Data/EBA_2011/Cont_C.csv")