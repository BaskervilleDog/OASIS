# syntax=docker/dockerfile:1

FROM alpine:3.20 AS builder

ARG BLAST_VERSION=2.13.0
ARG BLAST_URL=https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz
ARG DATASETS_URL=https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets

RUN apk add --no-cache ca-certificates curl tar gzip

WORKDIR /opt/oasis

RUN curl -fsSL "${DATASETS_URL}" -o datasets \
    && chmod +x datasets \
    && curl -fsSL "${BLAST_URL}" -o /tmp/blast.tar.gz \
    && tar -xzf /tmp/blast.tar.gz -C /opt/oasis \
    && rm /tmp/blast.tar.gz

FROM debian:bookworm-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        curl \
        grep \
        libbz2-1.0 \
        libgomp1 \
        libstdc++6 \
        tini \
        unzip \
        zlib1g \
    && rm -rf /var/lib/apt/lists/*

ENV HOME=/opt/oasis
ENV PATH=/opt/oasis/ncbi-blast-2.13.0+/bin:/opt/oasis:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

WORKDIR /data

COPY --from=builder /opt/oasis /opt/oasis
COPY OASIS.sh /app/OASIS.sh

RUN chmod +x /app/OASIS.sh

ENTRYPOINT ["tini", "--", "/app/OASIS.sh"]
