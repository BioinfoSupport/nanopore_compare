#!/bin/sh

podman build -t sniffles:2.2 -t sniffles:latest sniffles/ && \
    podman push sniffles:2.2 docker://docker.io/fedxa/sniffles:2.2 && \
    podman push sniffles:latest docker://docker.io/fedxa/sniffles:latest


podman build -t bioconductor:3.18 -t bioconductor:latest bioconductor/ && \
    podman push bioconductor:3.18 docker://docker.io/fedxa/bioconductor:3.18 && \
    podman push bioconductor:latest docker://docker.io/fedxa/bioconductor:latest


podman build -t nanopore_tools:0.5 -t nanopore_tools:latest nanopore_tools/ && \
    podman push nanopore_tools:0.5 docker://docker.io/fedxa/nanopore_tools:0.5 && \
    podman push nanopore_tools:latest docker://docker.io/fedxa/nanopore_tools:latest
