FROM python:3.10

WORKDIR /usr/src/app

RUN python -m pip install --upgrade pip
RUN pip install pandas statsmodels openpyxl

RUN mkdir -p /app/vintages
RUN mkdir -p /app/nowcast

COPY SA-Nowcast/code/python /app
COPY SA-Nowcast/vintages /app/vintages
COPY SA-Nowcast/nowcast /app/nowcast

CMD [ "python", "/app/nowcast_auto_econdata.py" ]
