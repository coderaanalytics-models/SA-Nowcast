#######################################
# Nowcasting Data EconData Clean
#######################################

options(repos = c(CRAN = "https://cran.mirror.ac.za"))
if(!requireNamespace("fastverse", quietly = TRUE)) install.packages("fastverse")
if(!requireNamespace("econdatar", quietly = TRUE)) {
  if(!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  library("remotes")
  install_github("coderaanalytics/econdatar", ref = "3.1.0")
}
library(fastverse)
fastverse_extend(econdatar, seasonal, writexl, install = TRUE) # seastests -> many dependencies, also tsbox, xts, tseries


# Seasonally Adjusting Relevant Indicators
seasadj <- function(x) {
  cc <- whichNA(x, invert = TRUE)
  if(length(cc) < 10) return(rep(NA_real_, length(x)))
  x[cc] <- tryCatch(final(seas(x)), error = function(e) tryCatch(final(seas(x, outlier = NULL)),
                                                                 error = function(e2) forecast::seasadj(stl(na.omit(x), "periodic"))))
  x
}

# Imputing internal NA's for seasonal adjustment
spline_impute <- function(X) {
  for(i in seq_col(X)) {
    x <- X[, i]
    nnai <- whichNA(x, invert = TRUE)
    ln <- length(nnai)
    t1 <- nnai[1L]
    t2 <- nnai[ln]
    if (ln != t2 - t1 + 1L)
      x[t1:t2] <- spline(nnai, x[nnai], xout = t1:t2)$y
    X[, i] = x
  }
  X
}

# Helper for log-differencing
adjust_negative <- function(x) if(any((min <- fmin(x)) <= 0)) x %r+% iif(min <= 0, -min+1, 0) else x


#
### Final Selection of Series from EconData --------------------------------------------------------------------------------------------
#

econdata_monthly <- list(
  Real = list(
    Production = list( # Multiple R-squared:  0.3613
      KBP7085N = c("MANUFACTURING:MAN001.I.S",
                   "MANUFACTURING:MAN001.S.S"), # Total Manufacturing (Business cycles + Rand)
      KBP7062N = c("MINING:MIN001.I.S",
                   "MINING:MIN001.S.S"), # Total Mining Production, Seas. Adj.
      KBP7068N = c("ELECTRICITY:ELE002.I.S",
                   "ELECTRICITY:ELE001.S.S") # Electricity Generation and availability for distr. in SA # ELE003_S_N
    ),
    Sales = list( # Multiple R-squared:  0.1832
      KBP7067N = c("MOTOR_TRADE:MTS003.S",
                   "MOTOR_TRADE:MTS005.N"), # Total motor trade and new vehicle sales (replaces Number of vehicles sold, Seas Adj. Index)
      KBP7086T = c("RETAIL:RET008.I.S",
                   "RETAIL:RET008.S.S"), # Retail Sales
      KBP7087T = c("WHOLESALE:WHO001.I.S",
                   "WHOLESALE:WHO001.S.S") # Wholesale Sales
    ),
    Prices = list( # Multiple R-squared:  0.1085
      KBP7155N = c("CPI_ANL_SERIES:CPI60001",
                   "BUSINESS_CYCLES:CPI1000.M.N"), # CPI Headline
      KBP7198M = c("PPI:PPI001",
                   "PPI:PPI027",
                   "PPI:PPI028",
                   "PPI:PPI041") # Producer prices. Replaced with total, final manufactures, petrol and motor vehicles
    ),
    Tourism = list( # Multiple R-squared:  0.8654
      MIGRATION = c("MIGRATION:MIG001.A.N0.TA",
                    "MIGRATION:MIG001.A.A0.TA",
                    "MIGRATION:MIG011.N.A0.TX",
                    "MIGRATION:MIG011.N.N0.TX"), # Total + Total Air + Total overnight tourists + Air
      TOURIST_ACCOMMODATION = c("TOURIST_ACCOMMODATION:TOU036.S",
                                "TOURIST_ACCOMMODATION:TOU006.S",
                                "TOURIST_ACCOMMODATION:TOU011.S") # Total income, Stay Units Nights Sold and Occupancy rate
    ),
    `Other Real` = list( # Multiple R-squared:  0.4357
      KBP7090N = "BUSINESS_CYCLES:DIFN003.M.S", # Leading Indicator
      LAND_TRANSPORT = c("LAND_TRANSPORT:LAN001.S",
                         "LAND_TRANSPORT:LAN002.S",
                         "LAND_TRANSPORT:LAN018.S",
                         "LAND_TRANSPORT:LAN019.S") # Not really very current information...
    )
  ),
  Financial = list( # Multiple R-squared:  0.005111
    `Money and Credit` = list(
      FINANCIAL_SECTOR = "FINANCIAL_SECTOR:MON0088.M", # M0
      KBP1374M = "FINANCIAL_SECTOR:MON0300.M", # M3
      KBP1347M = "FINANCIAL_SECTOR:MON0023.M", # PSC
      KBP1367M = "FINANCIAL_SECTOR:MON0191.M" # Credit to the Government
    ),
    `Other Fiancial` = list( # Multiple R-squared:  0.01558
      FINANCIAL_SECTOR = "FINANCIAL_SECTOR:MON0263.M", # NFA
      LIQUIDATIONS = "LIQUIDATIONS:LIQ002.A.L.A.N" # Total Liquidations and Insolvencies (could ad more, any signal value??)
    )
  ),
  External = list(
    Trade = list( # Multiple R-squared:  0.2809
      EXTERNAL_SECTOR = c("EXTERNAL_SECTOR:CURX600.M",
                          "EXTERNAL_SECTOR:CURM600.M") # Exports and Imports
    ),
    `Exchange Rates` = list( # Multiple R-squared:  0.003367
      KBP5339M = "EXTERNAL_SECTOR:BOP5329.M", # Rand <-> USD Exchange Rate
      KBP5393M = "EXTERNAL_SECTOR:BOP5393.M" # NEER
    ),
    Reserves = list( # Multiple R-squared:  0.02715
      KBP1021M = "EXTERNAL_SECTOR:BOP5806.M", # Total reserves
      EXTERNAL_SECTOR = "EXTERNAL_SECTOR:BOP5272.M" # Foreign currency reserves
    )
  ),
  Fiscal = list(
    `Cash Flow` = list( # Multiple R-squared:  0.8279
        KBP4597M = "FISCAL_SECTOR:NGFC020.M", # Total Revenue (Replaced with Cash Flow Revenue)
        KBP4601M = "FISCAL_SECTOR:NGFC040.M", # Total Expenditure (Replaced with Cash Flow Expenditure)
        KBP4050M = "FISCAL_SECTOR:NGFC050.M" # Cash Flow Balance
    ),
    Financing = list( # Multiple R-squared:  0.01947
        KBP4022M = "FISCAL_SECTOR:NGFC102.M", # Financing: Domestic Government Bonds
        KBP4026M = "FISCAL_SECTOR:NGFC103.M", # Financing: Foreign Bonds and Loans
        KBP4023M = "FISCAL_SECTOR:NGFC101.M", # Financing: Treasury Bills and Short-Term Loans
        KBP4003M = "FISCAL_SECTOR:NGFC006.M", # Financing: Change in Cash Balances
        KBP4030M = "FISCAL_SECTOR:NGFC100.M" # Total financing of national government
    ),
    Debt = list( # Multiple R-squared:  0.005562
        KBP4114M = "FISCAL_SECTOR:NGD1213.M", # Total loan debt of national government: Total gross loan debt
        KBP4105M = "FISCAL_SECTOR:NGD1209.M", # Total loan debt of national government: Total domestic debt (replaced with marketable debt)
        KBP4108M = "FISCAL_SECTOR:NGD7900.M", # Total loan debt of national government: Total foreign debt
        FISCAL_SECTOR = "FISCAL_SECTOR:NGD4500.M" # Domestic non-marketable debt
    )
  )
)

# Note: these series are renamed, so the that DFM models reading the excel file continue to work, even if we change the data source
econdata_quarterly <- list(Real = list(`Other Real`= list(BUSINESS_CYCLES = c(UNEMP = "BUSINESS_CYCLES:LABT079.Q.S")),
                                        Production = list(NATL_ACC = c(GDP = "NATL_ACC:KBP6006.N.S",
                                                                       RGDP = "NATL_ACC:KBP6006.R.S")))) # KBP6006_R_N,

#
### Creating Nowcasting Dataset ----------------------------------------------------------------------------------------------
#

nc_ind <- c(econdata_monthly, econdata_quarterly) %>%
  rapply(qDF, how = "list") %>%
  unlist2d(c("broad_sector", "topic", "QB"), "series_alt", DT = TRUE) %>%
  ftransform(series = X, QB = NULL, X = NULL)

nc_data <- nc_ind %>%
  fmutate(data = apply(do.call(rbind, strsplit(series, ":")), 1,
                       function(x) read_dataset(x[1],
                                                series_key = x[2],
                                                tidy = TRUE,
                                                wide = FALSE,
                                                prettymeta = FALSE)),
          updated = apply(do.call(rbind, strsplit(series, ":")), 1,
                          function(x) read_release(x[1], tidy = TRUE)$release[1]),
          metadata = lapply(data, function(x) x$metadata),
          data = lapply(data, function(x) x$data))
  
nc_series <- cbind(nc_data, rbindlist(nc_data$metadata, fill = TRUE)) %>%
  frename(data_set_name = dataset,
          LABEL = label,
          UNIT_MULT = unit,
          FREQ = freq,
          SEASONAL_ADJUST = seas_adj,
          SOURCE_IDENTIFIER = src_code,
          COMMENT = comment) %>%
  fmutate(updated = as.Date(as.POSIXct(updated)),
          from_date = as.Date(sapply(data, function(x) head(x, 1)$time_period)),
          to_date = as.Date(sapply(data, function(x) tail(x, 1)$time_period)),
          n_obs = sapply(data, nrow),
          dsid = do.call(rbind, strsplit(series, ":"))[,1],
          series = gsub("\\.", "_", do.call(rbind, strsplit(series, ":"))[,2]),
          series_orig = series,
          seas_adj = seas_adj == "S") %>%
  fselect(dsid, dataset, series, label, freq, unit, seas_adj,
          n_obs, from_date, to_date, src_code, comment,
          updated, broad_sector, topic, series_alt)

settransform(nc_series,
   minimal = topic %in% c("Production", "Sales", "Prices", "Tourism", "Other Real", "Trade", "Cash Flow"),
   series_orig = series,
   series = iif(nchar(series_alt) > 2L, series_alt, series),
   series_alt = NULL
)

names_alt <- nc_series %$% set_names(series, series_orig)[series_orig != series]

nc_data_m <- rbindlist(nc_data$data[sapply(nc_data$metadata, function(x) x$FREQ) == "M"]) %>%
  fmutate(series_key = gsub("\\.", "_", series_key)) %>%
  frename(time_period = date) %>%
  pivot(how = "wider", values = "obs_value", names = "series_key") %>%
  roworderv("date", decreasing = FALSE)

nc_data_q <- rbindlist(nc_data$data[sapply(nc_data$metadata, function(x) x$FREQ) == "Q"]) %>%
  fmutate(series_key = gsub("\\.", "_", series_key)) %>%
  frename(time_period = date) %>%
  pivot(how = "wider", values = "obs_value", names = "series_key") %>%
  roworderv("date", decreasing = FALSE) %>%
  frename(names_alt, .nse = FALSE)

# Adjusting Monthly Indicators
# nc_data_m %>% num_vars() %>% sapply(isSeasonal, freq = 12) %>% which() %>% names()
# -> not run because seastest has many dependencies
nc_seas_m <- c("MTS005_N", "CPI60001", "CPI1000_M_N", "MIG001_A_N0_TA", "MIG001_A_A0_TA",
  "MIG011_N_A0_TX", "MIG011_N_N0_TX", "MON0088_M", "MON0300_M",
  "MON0023_M", "MON0191_M", "LIQ002_A_L_A_N", "CURX600_M", "CURM600_M",
  "NGFC020_M", "NGFC040_M", "NGFC050_M", "NGFC102_M", "NGFC101_M",
  "NGFC006_M", "NGFC100_M", "NGD1213_M", "NGD1209_M", "NGD4500_M")
nc_data_sa_ts <- nc_data_m %>% get_vars(nc_seas_m) %>% qM() %>%
  ts(start = c(year(nc_data_m$date[1L]), month(nc_data_m$date[1L])), frequency = 12) %>%
  spline_impute()

for (i in seq_col(nc_data_sa_ts)) {
  cat(colnames(nc_data_sa_ts)[i], "\n")
  nc_data_sa_ts[, i] <- seasadj(nc_data_sa_ts[, i])
}

get_vars(nc_data_m, colnames(nc_data_sa_ts)) <- mctl(nc_data_sa_ts)
nc_series[series %in% colnames(nc_data_sa_ts), seas_adj := TRUE]

# Transformed datasets
nc_rates_m <- nc_series[freq == "M" & unit %ilike% "Percentage", series]
nc_data_m_logdiff <- nc_data_m %>%
  ftransform(fselect(., -date) %>%
               adjust_negative() %>%
               tfmv(nc_rates_m, function(x) x/100 + 1) %>%
               fgrowth(logdiff = TRUE))

nc_rates_q <- nc_series[freq == "Q" & unit %ilike% "Percentage", series]
nc_data_q_logdiff <- nc_data_q %>%
  ftransform(fselect(., -date) %>%
               adjust_negative() %>%
               tfmv(nc_rates_q, function(x) x/100 + 1) %>%
               fgrowth(logdiff = TRUE))

# Percent of distinct values
print(sort(fndistinct(nc_data_m_logdiff)/fnobs(nc_data_m_logdiff)) * 100)
print(sort(fndistinct(nc_data_q_logdiff)/fnobs(nc_data_q_logdiff)) * 100)

# Correlations with quarterly indicators
nc_corrs <- nc_data_m_logdiff %>% merge(nc_data_q_logdiff, by = "date") %>%
  fselect(-date) %>% pwcor(., gv(., names_alt)) %>% qDT("series") %>%
  add_stub("corr_", cols = -1) %>% frename(tolower) %>%
  fmutate(avg_abs_corr = pmean(abs(corr_unemp), abs(corr_gdp), abs(corr_rgdp)))

if (!all(nc_corrs$series %in% nc_series$series)) stop("missing series")
nc_series = nc_series[nc_corrs, on = "series"]

# Saving Data
list(series = nc_series,
     data_m = nc_data_m,
     data_q = nc_data_q,
     data_logdiff_m = nc_data_m_logdiff,
     data_logdiff_q = nc_data_q_logdiff) %>%
  write_xlsx(sprintf("vintages/econdata_nowcast_data_%s.xlsx", format(Sys.Date(), "%d_%m_%Y")))
