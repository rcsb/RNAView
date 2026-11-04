# --- Build Stage ---
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    && rm -rf /var/lib/apt/lists/*

# Clone and compile RNAView
WORKDIR /opt
RUN git clone https://github.com/rcsb/RNAView rnaview-src
WORKDIR /opt/rnaview-src
RUN make

# --- Runtime Stage ---
FROM ubuntu:22.04

# Copy compiled binaries and required reference folders
WORKDIR /opt/RNAView
COPY --from=builder /opt/rnaview-src/bin ./bin
COPY --from=builder /opt/rnaview-src/BASEPARS ./BASEPARS

# Install additional dependencies
RUN apt-get update && apt-get install -y \
    inkscape \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Set the required RNAView environment variables
ENV RNAVIEW=/opt/RNAView
ENV PATH="/opt/RNAView/bin:${PATH}"

# Set a default directory for user input/output data
WORKDIR /data

# Default to showing the help menu
CMD ["rnaview", "-h"]
