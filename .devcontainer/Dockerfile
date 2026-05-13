FROM mcr.microsoft.com/devcontainers/base:ubuntu-24.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    clangd \
    clang-format \
    cmake \
    gdb \
    git \
    libboost-all-dev \
    libopenbabel-dev \
    libtbb-dev \
    ninja-build \
    openbabel \
    pkg-config \
    python3 \
    ripgrep \
    && rm -rf /var/lib/apt/lists/*
