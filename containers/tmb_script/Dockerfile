# Multi-stage build using official pixi image
FROM python:3.12-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends build-essential gcc

RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Copy pixi configuration and modify for linux platform only
COPY requirements.txt .
RUN pip install -r requirements.txt

# Final stage - minimal runtime image
FROM python:3.12-slim

# Copy pixi environment from builder
COPY --from=builder /opt/venv /opt/venv

# Set working directory
WORKDIR /opt/venv/bin/

# Copy the TMB calculator script
COPY calculate_tmb.py gtf_to_coding_bed.py .
COPY gencode_46_coding.bed.gz /opt/

# Set environment variables for pixi
ENV PATH="/opt/venv/bin:$PATH"
ENV PYTHONPATH="/opt/venv/lib/python3.11/site-packages"

# Set entrypoint to the TMB calculator
ENTRYPOINT ["python", "calculate_tmb.py"]
CMD ["--help"]
