# Multi-stage Docker build for Owl MSI analysis tool
# Stage 1: Build stage
FROM rust:1.89-slim AS builder

# Install system dependencies required for rust-htslib
RUN apt-get update && apt-get install -y \
    pkg-config \
    libclang-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    build-essential \
    g++ \
    cmake \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# Copy and extract source code
COPY owl-0.2.0.tar.gz .
RUN tar -xzf owl-0.2.0.tar.gz && mv owl-0.2.0/* . && rm -rf owl-0.2.0 owl-0.2.0.tar.gz

# Build the application in release mode with optimizations
ENV RUSTFLAGS="-C target-cpu=x86-64 -C strip=symbols"
RUN cargo build --release --target x86_64-unknown-linux-gnu

# Stage 2: Runtime stage
FROM debian:bookworm-slim AS runtime

# Install only runtime dependencies
RUN apt-get update && apt-get install -y \
    libcurl4 \
    libbz2-1.0 \
    liblzma5 \
    zlib1g \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Create non-root user for security
RUN useradd -m -u 1001 owl

# Copy the binary from builder stage
COPY --from=builder /build/target/x86_64-unknown-linux-gnu/release/owl /usr/local/bin/owl

# Copy the data file
COPY --from=builder /build/data/Simple-repeats-50k.filt.bed.gz /opt/owl/data/

# Set ownership and permissions
RUN chown -R owl:owl /opt/owl && chmod +x /usr/local/bin/owl

# Switch to non-root user
USER owl

# Set working directory
WORKDIR /opt/owl

# Set default command
ENTRYPOINT ["/usr/local/bin/owl"]
CMD ["--help"]
