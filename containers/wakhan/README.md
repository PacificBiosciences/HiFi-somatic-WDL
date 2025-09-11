# Wakhan Docker Image

This directory contains the Docker configuration for building a containerized version of [Wakhan](https://github.com/KolmogorovLab/Wakhan), a tool for analyzing haplotype-specific chromosome-scale somatic copy number aberrations and aneuploidy using long reads (Oxford Nanopore, PacBio).

## Image Details

- **Base Image**: `mambaorg/micromamba:1.5.3`
- **Wakhan Version**: 0.1.2-3 (commit: `52d59c56be547bd93245a855c7aa5b6da35da691`)
- **Python Version**: 3.8
- **Platform**: linux/amd64 (x86_64)
- **Final Image Size**: ~1.39GB

## Building the Image

### Prerequisites

- Docker with buildx support
- Internet connection for downloading dependencies

### Build Commands

Build the image for x86_64 architecture:

```bash
docker buildx build --platform linux/amd64 -t wakhan:latest .
```

Build with custom tag:

```bash
docker buildx build --platform linux/amd64 -t your-registry/wakhan:v0.1.2 .
```
