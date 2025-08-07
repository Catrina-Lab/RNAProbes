FROM python:3.9 AS builder

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /app

RUN pip install poetry
RUN poetry config virtualenvs.in-project true
COPY pyproject.toml poetry.lock ./
RUN poetry install
FROM python:3.9-slim
WORKDIR /app
COPY --from=builder /app/.venv .venv/
COPY . .

ENV IS_WEB_APP=TRUE
RUN find /app/src/RNAStructure_Binaries/Linux64 -type f -exec chmod +x {} \;
CMD ["/app/.venv/bin/python3", "-m", "flask", "run", "--host=0.0.0.0", "--port=8080"]
