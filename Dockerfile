# Use an appropriate base image for aarch64
FROM arm64v8/ubuntu:20.04

# Set a working directory
WORKDIR /app

# Install system dependencies and ncbi-blast+
RUN apt-get update && \
    apt-get install  -y curl gcc build-essential ncbi-blast+ && \
    apt-get clean

# Download and install Miniforge (a minimal conda installer)
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh" && \
    bash Miniforge3-Linux-aarch64.sh -b -p /opt/conda && \
    rm Miniforge3-Linux-aarch64.sh

# Add Conda to the PATH
ENV PATH=/opt/conda/bin:$PATH

# Copy your environment.yml
COPY environment.yml /app/

# Install Conda dependencies using mamba
RUN /opt/conda/bin/conda install -c conda-forge mamba -n base --yes && \
    mamba env create -f environment.yml

# Copy your application code
COPY . /app

# Command to run your application
CMD ["conda", "run", "--no-capture-output", "-n", "bio", "python", "wsgi.py"]
