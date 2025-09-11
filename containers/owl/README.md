# Owl Docker Image

This directory contains a Docker build for [Owl](https://github.com/PacificBiosciences/owl), a tool for microsatellite instability (MSI) analysis for HiFi sequencing data.

## Overview

Owl analyzes microsatellite instability from PacBio HiFi sequencing data by profiling repeats in BAM files and scoring samples for MSI characteristics.

## Docker Image Details

- **Base Image:** Debian Bookworm Slim
- **Final Image Size:** ~99MB
- **Architecture:** linux/amd64
- **User:** Non-root user `owl` (UID: 1001)
- **Working Directory:** `/opt/owl`

## Prerequisites

- Docker with buildx support
- `owl-0.1.2.tar.gz` source archive in this directory

## Building the Image

### Quick Build
```bash
docker buildx build --platform linux/amd64 -t owl:latest .
```

### Build with Custom Tag
```bash
docker buildx build --platform linux/amd64 -t your-registry/owl:0.1.2 .
```

## Multi-Stage Build Process

The Dockerfile uses a multi-stage build for optimal size and security:

1. **Builder Stage (`rust:1.82-slim`)**
   - Installs build dependencies (cmake, g++, pkg-config, etc.)
   - Extracts and compiles Rust source code
   - Builds optimized release binary

2. **Runtime Stage (`debian:bookworm-slim`)**
   - Installs only runtime dependencies
   - Copies compiled binary and data files
   - Runs as non-root user

## Usage

### Basic Usage
```bash
# Show help
docker run --rm owl:latest
```

## Included Data

The image includes the default simple repeat catalog at:
- `/opt/owl/data/Simple-repeats-50k.filt.bed.gz`

## Image Contents

### Binary Location
- **Owl executable:** `/usr/local/bin/owl`
