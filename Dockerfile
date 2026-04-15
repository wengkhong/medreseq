# ── Stage 1: build ────────────────────────────────────────────────────────────
# Compile pysam (needs htslib headers) and install all Python deps into a venv.
FROM python:3.12-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc g++ \
        zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Isolated venv so it can be copied cleanly to the runtime stage
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install Python dependencies first (layer-cached unless pyproject.toml changes)
COPY pyproject.toml /build/pyproject.toml
WORKDIR /build
RUN pip install --no-cache-dir "pysam>=0.22" "primer3-py>=2.0" "click>=8.0"

# Install the medreseq package itself
COPY medreseq/ /build/medreseq/
RUN pip install --no-cache-dir --no-deps .


# ── Stage 2: runtime ──────────────────────────────────────────────────────────
# Minimal image — only the shared libraries pysam needs at runtime.
FROM python:3.12-slim AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
        zlib1g libbz2-1.0 liblzma5 libcurl4 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Non-root user for safety
RUN useradd -m -u 1000 medreseq
USER medreseq

# /data is the working directory; mount reference and I/O volumes here.
WORKDIR /data

ENTRYPOINT ["medreseq"]
CMD ["--help"]
