# UQFF Star-Magic — Docker image
#
# Build:  docker build -t uqff:5.27.1 .
# Run:    docker run --rm uqff:5.27.1 version
#         docker run --rm uqff:5.27.1 predict hubble_tension --json
#         docker run --rm uqff:5.27.1 list --filter neutrino
#         docker run --rm uqff:5.27.1 gate          # run the fidelity gate
#
# Interactive Python shell:
#         docker run --rm -it --entrypoint python uqff:5.27.1
#
# Mount a local notebook:
#         docker run --rm -it -p 8888:8888 -v $(pwd):/work uqff:5.27.1 \
#             --entrypoint jupyter notebook --ip=0.0.0.0 --no-browser

# Use the official slim Python 3.12 base image for a small footprint.
FROM python:3.12-slim AS base

LABEL org.opencontainers.image.title="UQFF Star-Magic"
LABEL org.opencontainers.image.description="Vacuum-first unified physics framework"
LABEL org.opencontainers.image.version="5.27.1"
LABEL org.opencontainers.image.authors="Daniel T. Murphy <daniel.murphy00@enrgyone.com>"
LABEL org.opencontainers.image.source="https://github.com/Daniel8Murphy0007/Star-Magic"
LABEL org.opencontainers.image.licenses="AGPL-3.0-or-later"

# Install only the bare essentials. The pure-Python uqff package has no
# runtime dependencies beyond the standard library.
WORKDIR /uqff

# Install UQFF from PyPI (the published package).
RUN pip install --no-cache-dir --upgrade pip \
 && pip install --no-cache-dir uqff==5.27.1

# Create a non-root user (security best practice).
RUN useradd --uid 1000 --create-home --shell /bin/bash uqff
USER uqff
WORKDIR /home/uqff

# Default entrypoint is the `uqff` CLI. Pass subcommands as arguments.
ENTRYPOINT ["uqff"]
CMD ["version"]
