FROM rocker/r-ver:4

RUN apt-get update
RUN apt-get -y install libcurl4-openssl-dev libssl-dev zlib1g-dev

RUN echo 'options(repos = c(CRAN = "https://cran.mirror.ac.za"))' >> "${R_HOME}/etc/Rprofile.site"


RUN Rscript -e "install.packages(c('fastverse', 'remotes', 'seasonal', 'writexl'));"

RUN Rscript -e "library(remotes); install_github('coderaanalytics/econdatar', ref = '3.1.0');"

RUN mkdir -p /app/vintages
RUN mkdir -p /app/nowcast

COPY SA-Nowcast/code/R /app
COPY SA-Nowcast/vintages /app/vintages
COPY SA-Nowcast/nowcast /app/nowcast

CMD [ "Rscript", "/app/nowcast_data_econdata_clean.R" ]
