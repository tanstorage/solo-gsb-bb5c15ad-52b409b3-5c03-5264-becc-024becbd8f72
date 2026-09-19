FROM python:3.12-slim

ENV PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /app

# Application code and dependency manifest; dependencies themselves are
# installed when the container starts (see start.sh).
COPY requirements.txt ./
COPY app ./app
COPY start.sh ./

RUN chmod +x start.sh

EXPOSE 8000

CMD ["./start.sh"]
