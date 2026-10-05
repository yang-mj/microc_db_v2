FROM python:3.12-slim-bookworm AS builder
ARG WITH_CONTACT=0
RUN apt-get update && apt-get install -y --no-install-recommends build-essential libcurl4-openssl-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY pyproject.toml README.md constraints.txt ./
COPY microc_explorer ./microc_explorer
RUN if [ "$WITH_CONTACT" = "1" ]; then pip wheel --wheel-dir /wheels -c constraints.txt '.[contact]'; \
    else pip wheel --wheel-dir /wheels -c constraints.txt .; fi

FROM python:3.12-slim-bookworm
ENV PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/microc-matplotlib
RUN apt-get update && apt-get install -y --no-install-recommends libcurl4 ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && useradd --create-home --uid 1000 appuser
COPY --from=builder /wheels /wheels
RUN pip install --no-cache-dir --no-index /wheels/*.whl && rm -rf /wheels
WORKDIR /app
COPY app.py config.json README.md ./
COPY .streamlit ./.streamlit
RUN chmod -R a+rX /app
USER appuser
EXPOSE 8501
HEALTHCHECK --interval=30s --timeout=5s CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8501/_stcore/health',timeout=3)"
CMD ["microc-web", "--host", "0.0.0.0", "--port", "8501"]
