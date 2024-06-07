# %% 
# Import packages
from flask import Flask
from dash import Dash, html, dash_table, dcc, callback, Output, Input
import plotly.express as px
import plotly.graph_objects as go
# %%
import dash_bootstrap_components as dbc
import numpy as np

# %% 
import pandas as pd
nowcast = pd.read_csv("./nowcast/nowcast.csv")    # index_col="date", parse_dates=True
gdp_ld = pd.read_csv("./nowcast/gdp_logdiff.csv") # index_col="quarter", parse_dates=True
gdp_ld.index = pd.PeriodIndex(gdp_ld.quarter, freq="Q")
gdp_ld = gdp_ld.loc[gdp_ld.index >= pd.PeriodIndex(nowcast.quarter, freq="Q").min()]
news = pd.read_csv("./nowcast/news.csv")
series = pd.read_csv("./nowcast/series.csv")
news = news.merge(series, 
                  left_on = "updated variable", 
                  right_on = "series", how = "left")

news_labels = {"Series" : "series", 
               "Dataset" : "dataset",
               "Label" : "label",
               "Release" : "observed",
               "Forecast" : "forecast (prev)",  # "date", "update date",
               "News" : "news",
               "Weight" : "weight",
               "Impact" : "impact",
               "Sector" : "broad_sector",
               "Topic" : "topic"}

news_labels_rev = {v:k for k, v in news_labels.items()}
news_labels_df = pd.DataFrame(news_labels, index = ["1"])

# format this to month-day
q = nowcast.quarter.max()
nowcast_latest_quarter = nowcast.loc[nowcast.quarter == q]
nowcast_dates = dict(zip(nowcast_latest_quarter.date, 
                         pd.to_datetime(nowcast_latest_quarter.date).dt.strftime("%b-%d")))
all_nowcast_dates = list(nowcast.date)



# %%
nowcast_final = nowcast.copy()
nowcast_final.date = pd.to_datetime(nowcast_final.date)
final_ids = nowcast_final.groupby('quarter').date.idxmax()
nowcast_final = nowcast.iloc[final_ids] # nowcast.groupby("quarter").last().reset_index()
nowcast_other = nowcast.drop(index = final_ids)


# %% 
news["sector_topic"] = news.broad_sector + ": " + news.topic

# %% 
# external_stylesheets = ['https://raw.githubusercontent.com/plotly/dash-app-stylesheets/master/dash-docs-base.css']
# external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
# external_stylesheets = ['https://raw.githubusercontent.com/tcbegley/dash-bootstrap-css/main/dist/cyborg/bootstrap.css']
server = Flask(__name__)
app = Dash(server = server, external_stylesheets=[dbc.themes.CYBORG]) #[dbc.themes.CYBORG]) #  # QUARTZ

app.title = "SA Nowcast" # Set website title
# # https://dash.plotly.com/dash-daq/darkthemeprovider
# theme = {
#     'dark': True,
#     'detail': '#007439',
#     'primary': '#00EA64',
#     'secondary': '#6E6E6E',
# }

# app.layout = daq.DarkThemeProvider(theme = theme, children = 
app.layout = dbc.Container([
    # html.Div(children='My First App with Data'),
    html.Br(),
    html.H3("South Africa Nowcast"),
    html.Hr(),
    dbc.Row([
        html.Div([
            html.H6("Select a Variable to Nowcast : "),
            dcc.RadioItems(options=[{'label': 'Real GDP', 'value': 'RGDP'}, 
                                {'label': 'Nominal GDP', 'value': 'GDP'}, 
                                {'label': 'Unemployment', 'value': 'UNEMP'}], 
                                value='RGDP', id='nc-variable', 
                    inline=True, inputStyle={"margin-right": "15px", "margin-left": "30px"})
        ], className = "select-var-block")
    ]),
    html.Hr(),
    # Add a navbar with 2 tabs
    dcc.Tabs(
        id="tabs-with-classes",
        # value='tab-2',
        parent_className='custom-tabs',
        className='custom-tabs-container',
        # className="dash-bootstrap",
        children=[
        dcc.Tab(label='Latest Nowcast Quarter', children=[
            # dash_table.DataTable(data=nowcast.to_dict('records'), page_size=10),
            html.Br(),
            html.H5("Nowcast Evolution for Latest Quarter"),
            html.Hr(),
            dbc.Row([
                dcc.Graph(figure = {}, id='nowcast-qx')
            ]),
            html.Br(),
            html.Div([
                html.H5("News Releases ", style={"margin" : "0px", "padding" : "0px"}), 
                html.Span(" ", style={"width" : "50px", "margin" : "0px", "padding" : "0px"}),
                html.P("[ Impact  =  News  x  Weight  =  (Release - Forecast)  x  Weight ]", style={"margin" : "0px", "padding" : "0px"})
            ], className = "select-var-block"),
            # html.H5("News Releases [Impact = News x Weight = (Release - Forecast) x Weight]"),
            html.Hr(),
            dbc.Row([
                # # make a select input for the vinage of the nowcast
                # dcc.Dropdown(options=nowcast_dates, value = list(nowcast.date)[-1], id='nc-date', 
                #             style={"margin-bottom": "10px"}),
                dcc.RadioItems(options=nowcast_dates, value = all_nowcast_dates[-1], id='nc-date', 
                               inline=True, inputStyle={"margin-right": "5px", "margin-left": "20px"}, 
                               style={"margin-bottom": "10px"}),
                dash_table.DataTable(data = news_labels_df.to_dict("records"), 
                                     columns = [{"name": i, "id": i} for i in news_labels_df.columns],
                                     page_size=50, id='nowcast-qx-news', 
                                     style_table={'overflowX': 'scroll'},
                                     style_header={
                                         'backgroundColor': 'rgb(30, 30, 30)',
                                         'border': "0px",
                                         'fontWeight': 'bold'
                                         },
                                     style_cell={
                                         'backgroundColor': '#111111',
                                         'border': "0px",
                                         'color': 'white',
                                         'padding-right': '15px'
                                         }
                                     )
            ]) #,
            # html.Br()
        ]),
        dcc.Tab(label='All Nowcasts', children=[
                html.Br(),
                html.H5("All Nowcasts and News Releases"),
                html.Hr(),
                html.Div([
                    html.H6("Restrict Nowcasting / News Aggregation Range : "),
                    dcc.DatePickerRange(
                        id='nowcasts-date-picker-range',
                        month_format = 'D MMM YYYY',
                        display_format='DD/MM/YYYY',
                        # min_date_allowed = all_nowcast_dates[0],
                        start_date = all_nowcast_dates[0],
                        # max_date_allowed = all_nowcast_dates[-1],
                        end_date = all_nowcast_dates[-1],
                        style={"margin-bottom": "10px", "margin-left": "30px"}
                    )], className="select-var-block"),
                # dcc.RangeSlider(
                #     min=0,
                #     max=len(all_nowcast_dates)-1,
                #     step=None,
                #     marks= dict(zip(range(len(all_nowcast_dates)), 
                #                 pd.to_datetime(nowcast.date).dt.strftime("%b-%d %Y"))) #,
                #     # value=[all_nowcast_dates[0], all_nowcast_dates[-1]]
                # ),
                dbc.Row([
                    dbc.Col([
                        dcc.Graph(figure = {}, id='all-nowcasts-ts')
                    ], width=6),
                    dbc.Col([
                        dcc.Graph(figure = {}, id='all-nowcasts-news')
                    ], width=6),
                ])
        ]),
        dcc.Tab(label='About the Nowcast', children=[
            # html.Hr(),
            html.Br(),
            html.H5("About the SA Nowcast"),
            html.P(["The South Africa Nowcast is a project that aims to provide a timely and accurate estimate of the current state of the South African economy. It is updated on a weekly basis and released every Friday."]),
            html.P(["The nowcast draws its data from ", html.A("EconData", href='https://www.econdata.co.za'), ". EconData is a platform that enables automation of analytical workflows that depend on public domain or third-party data. It is also a leading-edge forecast management system – enabling data and model automation, within a best practice data and model governance framework. EconData supports data-sharing across databases and within institutions, codifies modelling process flows and provides user-level access control. EconData makes it easy to securely manage and share model scenarios and forecast vintages."]),
            html.P(["The nowcast is based on a mixed-frequency dynamic factor model following Banbura & Modugno (2014) and Bok et al. (2018) (see ", 
                    html.A("New York Fed Nowcasting Model", href='https://www.newyorkfed.org/research/policy/nowcast'), ") which is estimated using 54 monthly and 3 quarterly series. The model uses a Kalman filter, which allows for the inclusion of new data as it becomes available. All monthly/quarterly series are seasonally adjusted (using X13) and transformed to monthly/quarterly growth rates (in percentage terms) via log-differencing. Model-based news is computed as the difference between a new data release and its model-based forecast from the previous period (i.e. news = actual minus predicted growth rate of updated series). The impact of this news on the nowcast is given by a model-based weight following Banbura & Modugno (2014), such that nowcast revision = weight x news."]),        
             html.P(["The source code and data is publically available on ",
                     html.A("GitHub", href = "https://github.com/Stellenbosch-Econometrics/SA-Nowcast"), " which includes weekly vintages of the nowcasting dataset available as Excel files and all nowcasts and news releases as CSV files."]), 
            html.P(["More information about the database and nowcasting methodology is provided in the accompanying ", 
                    html.A("presentation slides", href = "https://raw.githubusercontent.com/Stellenbosch-Econometrics/SA-Nowcast/main/presentation/SAMADB_Nowcasting.pdf"), "."]),
            html.Br(),
            html.H5("Authors"), 
            html.P(["The database and nowcasting model was built by ", 
                    html.A("Sebastian Krantz", href = "https://github.com/SebKrantz"), 
                    ". It is generously hosted by ", 
                    html.A("Codera Analytics", href = "https://codera.co.za/"), 
                    ", which also maintains ", 
                    html.A("EconData", href = "https://www.econdata.co.za/"), ". Credit is also due to ", 
                    html.A("Chad Fulton", href = "https://github.com/ChadFulton"), " who implemented the routines to estimate dynamic factor models for nowcasting in the ",
                    html.A("statsmodels", href = "https://www.statsmodels.org/dev/statespace.html#dynamic-factor-models"), " Python library."]),
            html.Br(),
            html.H5("References"),
            html.P(["Banbura, M., & Modugno, M. (2014). Maximum likelihood estimation of factor models on datasets with arbitrary pattern of missing data. ",
                    html.I("Journal of Applied Econometrics, 29"), "(1), 133-160."]),
            html.P(["Bok, B., Caratelli, D., Giannone, D., Sbordone, A. M., & Tambalotti, A. (2018). Macroeconomic nowcasting and forecasting with big data. ",
                    html.I("Annual Review of Economics, 10"), ", 615-643."])
        ]),
    ]),
    html.Footer(["© 2023 Sebastian Krantz and Codera Analytics"], 
                style={'text-align': 'center', 'margin-top': '30px', 'margin-bottom': '20px', 'color': '#737373'})
], style={'max-width': '90%', 'margin': 'auto'})

@callback(
    Output('nowcast-qx', 'figure'),
    Output('nowcast-qx-news', 'data'),
    Input('nc-variable', 'value'),
    Input('nc-date', 'value')
)
def update_nccq_graphs(var, date):
    ## Nowcast for latest quarter
    news_qx = news.loc[(news.quarter == q) & (news["impacted variable"] == var)]
    news_latest_quarter = news_qx.groupby(["date", "sector_topic"]).agg({"impact": "sum"}) \
               .reset_index().sort_values("sector_topic", ascending=False)
    # news_latest_quarter.impact = news_latest_quarter.impact / 10
    fig_qx = px.bar(news_latest_quarter, x="date", y="impact", color="sector_topic", barmode="relative")
    fig_qx.add_trace(go.Scatter(x=nowcast_latest_quarter.date, y=nowcast_latest_quarter[var], 
                                line=dict(color="white"), mode='lines+markers', name='Nowcast'))
    # Edit the layout
    fig_qx.update_layout(title='Nowcast for ' + q,   
                        xaxis_title='Date', yaxis_title='Quarterly Log-Difference Growth Rate (%)',
                        legend_title = "Sector: Topic",
                        hovermode="x", hoverlabel = dict(namelength = -1), # https://github.com/plotly/plotly.js/issues/460
                        autosize=False, width=1200, height=500,
                        margin=dict(l=20, r=20, t=40, b=20), template="plotly_dark")
    # delete the hover template
    fig_qx.update_yaxes(hoverformat=".2f")
    fig_qx.update_traces(hovertemplate=None)
    
    # %%
    news_vars = list(news_labels.values())
    news_latest_quarter_dict = news_qx.loc[news_qx.date == date, news_vars]
    news_latest_quarter_dict[news_vars[3:8]] = news_latest_quarter_dict[news_vars[3:8]].transform(lambda x: x.round(3))
    news_latest_quarter_dict = news_latest_quarter_dict.rename(columns=news_labels_rev).to_dict('records')
    if len(news_latest_quarter_dict) == 0:
        news_latest_quarter_dict = [{v : None for k, v in news_labels_rev.items()}]
    return fig_qx, news_latest_quarter_dict

@callback(
    Output('all-nowcasts-ts', 'figure'),
    Output('all-nowcasts-news', 'figure'),
    Input('nc-variable', 'value'),
    Input('nowcasts-date-picker-range', 'start_date'),
    Input('nowcasts-date-picker-range', 'end_date')
)
def update_allnc_graphs(var, start_date, end_date):
    ## Adjusting the date range
    nowcast_final_range = nowcast_final.loc[(nowcast_final.date >= start_date) & (nowcast_final.date <= end_date)]
    nowcast_other_range = nowcast_other.loc[(nowcast_other.date >= start_date) & (nowcast_other.date <= end_date)]
    gdp_dates = gdp_ld.index.to_timestamp()
    gdp_ld_range = gdp_ld.loc[(gdp_dates >= start_date) & (gdp_dates <= end_date)]

    ## Time series plot of all nowcasts 
    nowcast_other_range["nc_date"] = nowcast_other_range.date.astype(str) + " : " + ((np.exp(nowcast_other_range[var]/100)**4-1)*100).round(2).astype(str) 
    nowcast_other_range["nc_date_all"] = nowcast_other_range.groupby("quarter").nc_date.transform(lambda x: "<br>   " + "   <br>   ".join(reversed(x.tolist())))

    fig = go.Figure()   

    fig.add_trace(go.Scatter(x=nowcast_final_range.quarter,
                             y=(np.exp(nowcast_final_range[var]/100).rolling(4).apply(np.prod, raw = True)-1)*100, 
                             customdata=nowcast_final_range.date, 
                             hovertemplate="<br>   %{customdata} : %{y:.2f}",
                             mode='lines+markers', marker=dict(color='cyan'), name='Nowcast (YoY%)'))

    fig.add_trace(go.Scatter(x=nowcast_other_range.quarter, y=(np.exp(nowcast_other_range[var]/100)**4-1)*100, 
                            customdata=nowcast_other_range.nc_date_all, 
                            # stackgroup = nowcast_other_range.quarter,
                            hovertemplate= "%{customdata}", # "%{y:.2f} (%{customdata})",
                            mode='markers', marker=dict(color='red'), name='Older Nowcasts'))
    
    fig.add_trace(go.Scatter(x=nowcast_final_range.quarter, y=(np.exp(nowcast_final_range[var]/100)**4-1)*100, 
                    customdata=nowcast_final_range.date, 
                    hovertemplate="<br>   %{customdata} : %{y:.2f}",
                    mode='lines+markers', marker=dict(color='orange'), name='Nowcast'))

    fig.add_trace(go.Scatter(x=gdp_ld_range.quarter, y=(np.exp(gdp_ld_range[var]/100)**4-1)*100, 
                             hovertemplate="%{y:.2f}", 
                             mode='lines+markers', marker=dict(color='green'), name='Outcome'))
    
    # Edit the layout
    fig.update_layout(title='All Nowcasts (+ Backtesting 2019Q2-2023Q1)', 
                      xaxis_title='Quarter', 
                      yaxis_title='Quarterly Annualised Growth Rate (%)', 
                      hovermode="x unified", 
                      # hoverlabel = dict(namelength = 0),
                      legend_traceorder="reversed",
                      autosize=False, width=750, height=500, margin=dict(l=20, r=20, t=40, b=20), 
                      template="plotly_dark")
    
    fig.update_xaxes(categoryorder='array', 
                     categoryarray=nowcast_final_range.quarter)
    
    # News digest
    news_tot = news.loc[(news.date >= start_date) & (news.date <= end_date) & (news["impacted variable"] == var)] \
                .groupby(["sector_topic", "broad_sector", "topic"]) \
                .agg({"weight": "mean", "impact": "mean"}).reset_index() 
    news_tot["abs_impact"] = news_tot.impact.abs() # * 100

    ## Barcharts of average news impacts
    tot_news_fig = go.Figure()
    tot_news_fig.add_bar(x=news_tot.sector_topic, y=news_tot.abs_impact * 4)
    tot_news_fig.update_layout(title='Average Absolute News Impact',      
                    xaxis_title='', yaxis_title='Average Absolute Impact',
                    barmode='stack', hovermode="x", autosize=False, width=750, height=500,
                    margin=dict(l=20, r=20, t=40, b=20), template="plotly_dark")
    return fig, tot_news_fig

if __name__ == '__main__':
    app.run_server(debug=True)
# 

# %%
