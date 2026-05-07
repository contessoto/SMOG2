FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC

RUN apt-get update && apt-get install -y --no-install-recommends \
    perl \
    libxml-simple-perl \
    libxml-validator-schema-perl \
    pdl \
    default-jre-headless \
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    build-essential \
    git \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace/SMOG2

# editable install of smog3 package
COPY . /workspace/SMOG2
RUN python3 -m pip install --break-system-packages -U pip setuptools wheel && \
    python3 -m pip install --break-system-packages -e .[dev] || \
    python3 -m pip install --break-system-packages -e .

# baseline perl dependency check required by parity harness
RUN perl -MXML::Simple -MXML::Validator::Schema -MPDL -e 'print "SMOG2 Perl deps ok\\n"'

CMD ["bash"]
