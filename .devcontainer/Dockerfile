# Start with the rocker/tidyverse base
FROM rocker/tidyverse:4.6

# Define the build argument for vscode users (defaults to false)
ARG INSTALL_RADIAN="false"

# Enable Ubuntu 'universe' repository for ease of package installation
RUN if [ -f /etc/apt/sources.list.d/ubuntu.sources ]; then \
    sed -i 's/Components: main restricted/Components: main restricted universe multiverse/g' /etc/apt/sources.list.d/ubuntu.sources; \
elif [ -f /etc/apt/sources.list ]; then \
    sed -i 's/main/main universe/g' /etc/apt/sources.list; \
fi

# Install system dependencies required for compiling common R/Bioinformatics packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Libraries for r packages
    libfontconfig1-dev \
    libxml2-dev \
    libcairo2-dev \
    libtiff-dev \
    libuv1-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libxt-dev \
    zlib1g-dev \
    # GLPK library required for igraph (dependency for clusterProfiler plots)
    libglpk-dev \
    # .NET Runtime requirements
    libicu-dev \
    libgssapi-krb5-2 \
    # Python
    python3 \
    python3-pip \
    # Downloading tools needed for setup scripts
    curl \
    ca-certificates \
    # Deb packages
    tmux \
    bcftools \
    && rm -rf /var/lib/apt/lists/*

# Install .NET Runtime (v6.0). Required to run Illumina Nirvana.
# Set .NET Environment Variables so the system can find the runtime
ENV DOTNET_ROOT=/usr/share/dotnet
ENV PATH="${PATH}:${DOTNET_ROOT}"

# Use MS install script.
RUN curl -fsSL https://dot.net/v1/dotnet-install.sh -o /tmp/dotnet-install.sh \
    && chmod +x /tmp/dotnet-install.sh \
    && /tmp/dotnet-install.sh --version 6.0.36 --runtime dotnet --install-dir /usr/share/dotnet

# Download and install Illumina Nirvana
ENV NIRVANA_VERSION=3.18.1
RUN mkdir -p /opt/nirvana \
    && curl -L https://github.com/Illumina/Nirvana/releases/download/v${NIRVANA_VERSION}/Nirvana-${NIRVANA_VERSION}-net6.0.zip -o /tmp/nirvana.zip \
    && unzip /tmp/nirvana.zip -d /opt/nirvana \
    && rm /tmp/nirvana.zip \
    # Create handy shell wrappers so you can call `nirvana` and `nirvana-downloader` globally
    && printf '#!/bin/bash\nexec dotnet /opt/nirvana/Nirvana.dll "$@"\n' > /usr/local/bin/nirvana \
    && printf '#!/bin/bash\nexec dotnet /opt/nirvana/Downloader.dll "$@"\n' > /usr/local/bin/nirvana-downloader \
    && chmod +x /usr/local/bin/nirvana /usr/local/bin/nirvana-downloader

# Conditionally install pipx and radian if INSTALL_RADIAN is "true"
# We configure pipx to install globally so all container users can access it
RUN if [ "$INSTALL_RADIAN" = "true" ]; then \
        apt-get update && apt-get install -y --no-install-recommends \
            pipx \
        && rm -rf /var/lib/apt/lists/* \
        && PIPX_HOME=/opt/pipx PIPX_BIN_DIR=/usr/local/bin pipx install radian \
        && R -e "install.packages(c('languageserver', 'httpgd'), repos='https://cloud.r-project.org')"; \
    fi

# Pre-install Bioconductor and key genomic packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" \
    && R -e "BiocManager::install(c('GenomicRanges', 'limma', 'enrichplot', 'clusterProfiler'))"

WORKDIR /workspace