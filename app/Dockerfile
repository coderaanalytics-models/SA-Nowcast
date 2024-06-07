FROM python:3.7

WORKDIR /usr/src/app

COPY ./app .
COPY ./nowcast ./nowcast

RUN pip install --no-cache-dir -r requirements.txt

EXPOSE 8050

CMD [ "gunicorn", "--workers=5", "--threads=1", "-b 0.0.0.0:8050", "app:server" ]
