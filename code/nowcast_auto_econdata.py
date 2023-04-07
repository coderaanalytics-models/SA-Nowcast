# %%
import os
import numpy as np
import pandas as pd
from statsmodels.tsa.api import DynamicFactorMQ
from datetime import datetime

# %%
def load_vintage(path):

    # Reading 
    series = pd.read_excel(path, sheet_name="series", engine='openpyxl')
    data = pd.read_excel(path, sheet_name="data_m", 
                        index_col="date", parse_dates=True, engine='openpyxl')
    data_logdiff = pd.read_excel(path, sheet_name="data_logdiff_m",
                                index_col="date", parse_dates=True, engine='openpyxl')
    gdp = pd.read_excel(path, sheet_name="data_q", 
                       index_col="date", parse_dates=True, engine='openpyxl')
    gdp_logdiff = pd.read_excel(path, sheet_name="data_logdiff_q", 
                        index_col="date", parse_dates=True, engine='openpyxl')
    
    # Subsetting 
    series_q = series.loc[series.freq == "Q"]
    series_m = series.loc[series.freq == "M"]
    data = data.loc[data.index >= "2000-01-01"]
    data_logdiff = data_logdiff.loc[data_logdiff.index >= "2000-01-01"]
    gdp = gdp.loc[gdp.index >= "2000-01-01"]
    gdp_logdiff = gdp_logdiff.loc[gdp_logdiff.index >= "2000-01-01"]

    # Period Index
    data.index = data.index.to_period()
    data_logdiff.index = data_logdiff.index.to_period()
    gdp.index = gdp.index.to_period()
    gdp_logdiff.index = gdp_logdiff.index.to_period()

    return dict(series = series, series_m = series_m, series_q = series_q,
                data = data, data_logdiff = data_logdiff, gdp = gdp, gdp_logdiff = gdp_logdiff)

# %%
vintages = os.listdir("vintages")
# %%
vintage_dates = {x: datetime.strptime(x[22:32], "%d_%m_%Y") for x in vintages}
# %%
latest_vintage = max(vintage_dates, key = vintage_dates.get)
print(latest_vintage)
# %%
def removekey(d, key):
    r = dict(d)
    del r[key]
    return r

previous_vintage = max(removekey(vintage_dates, latest_vintage), key = vintage_dates.get)
print(previous_vintage)
# %%
latest_data = load_vintage("vintages/" + latest_vintage)
previous_data = load_vintage("vintages/" + previous_vintage)

# %%
latest_data["series"].groupby('broad_sector', sort=False)["series"].count()
# %%
latest_data["series"][["series", "series_orig", "label", "freq", "unit", "seas_adj", "broad_sector", "topic", 
    "code_src", "comment"]].to_csv("nowcast/series.csv", index=False)
latest_data["gdp"].to_csv("nowcast/gdp.csv", index_label="quarter")
latest_data["gdp_logdiff"].to_csv("nowcast/gdp_logdiff.csv", index_label="quarter")
# %%
series = latest_data["series"]
labels = {k: v for k, v in zip(series.series, series.label)}
factors = {l: ['Global', v] for l, v in zip(series.series, series.broad_sector)}
# %%
dfm_latest = DynamicFactorMQ(endog=latest_data["data_logdiff"], 
                             endog_quarterly=latest_data["gdp_logdiff"][["UNEMP", "GDP", "RGDP"]], 
                             factors=factors, 
                             factor_multiplicities=dict(Global = 2, Real = 2, Financial = 1, Fiscal = 2, External = 2), 
                             factor_orders=2)
dfm_previous = DynamicFactorMQ(endog=previous_data["data_logdiff"], 
                               endog_quarterly=previous_data["gdp_logdiff"][["UNEMP", "GDP", "RGDP"]], 
                               factors=factors, 
                               factor_multiplicities=dict(Global = 2, Real = 2, Financial = 1, Fiscal = 2, External = 2), 
                               factor_orders=2)
# %%
dfm_latest_results = dfm_latest.fit(disp=10)
dfm_previous_results = dfm_previous.fit(disp=10)

# %% 
dfm_latest_results.summary()
# %%
today = datetime.date(max(vintage_dates.values())) # datetime.today()
today_q = pd.PeriodIndex([today], freq = "Q")[0]
print(f"\nToday is {str(today)}. We are nowcasting for Quarter {today_q}.")
# %% 
gdp_now = dfm_latest_results.get_prediction(start = today_q)
print(gdp_now.predicted_mean[["UNEMP", "GDP", "RGDP"]])
# %%
# gdp_now.summary_frame(endog=-1)
nowcast = gdp_now.predicted_mean[["UNEMP", "GDP", "RGDP"]].resample("Q").last()
nowcast["date"] = today
nowcast["quarter"] = str(today_q)
nowcast = nowcast[["date", "quarter", "UNEMP", "GDP", "RGDP"]]
print(nowcast)
# %% 
nowcast_old = pd.read_csv("nowcast/nowcast.csv")
nowcast_all = pd.concat([nowcast_old.loc[pd.to_datetime(nowcast_old.date) < pd.to_datetime(today)], 
                         nowcast]).reset_index(drop = True)
nowcast_all.tail()
# %% 
nowcast_all.to_csv("nowcast/nowcast.csv", index = False)
# %%
news = dfm_latest_results.news(dfm_previous_results, 
                               impact_date = str(today_q), 
                               impacted_variable = ["UNEMP", "GDP", "RGDP"], 
                               comparison_type = "previous")
# %%
news.summary()
# %%
news_df = news.details_by_impact.reset_index() 
news_df["impact date"] = pd.PeriodIndex(news_df["impact date"]).to_timestamp(freq="Q")
news_df["quarter"] = str(today_q)
news_df["date"] = today
news_df.head()
# %%
news_old = pd.read_csv("nowcast/news.csv")
# %%
news_all = pd.concat([news_old.loc[pd.to_datetime(news_old.date) < pd.to_datetime(today)], 
                      news_df[news_old.columns]])
news_all.tail()
# %%
news_all.to_csv("nowcast/news.csv", index = False)
