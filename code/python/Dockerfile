FROM python:3.10

WORKDIR /usr/src/app

RUN python -m pip install --upgrade pip
RUN pip install pandas statsmodels openpyxl

COPY ./code .
COPY ./vintages .
COPY ./nowcast .

CMD [ "python", "./code/python/nowcast_auto_econdata.py" ]
