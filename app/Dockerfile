FROM python:3.7

WORKDIR /usr/src/app

COPY SA-Nowcast/app/requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

COPY SA-Nowcast/app .

RUN mkdir vintages
RUN mkdir nowcast

EXPOSE 8050

CMD [ "gunicorn", "--workers=5", "--threads=1", "-b 0.0.0.0:8050", "app:server" ]
