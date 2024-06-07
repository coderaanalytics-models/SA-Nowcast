FROM rocker/r-ver:4

RUN echo 'options(repos = c(CRAN = "https://cran.mirror.ac.za"))' >> "${R_HOME}/etc/Rprofile.site"


Rscript -e "install.packages(c('fastverse', 'DBI', 'RMySQL', 'writexl', 'seasonal', 'remotes')); \
            library(fastverse); \
            library(remotes); \
            install_github('coderaanalytics/econdatar', ref = '3.0.2'); \
            fastverse_extend(econdatar, seasonal, writexl, install = TRUE)"

COPY ./code .
COPY ./vintages .

CMD ["Rscript" "./code/R/nowcast_data_econdata_clean.R"]
