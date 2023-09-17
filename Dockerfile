FROM ubuntu:18.04
SHELL ["/bin/bash", "-c"]

# Set some timezone or otherwise tzdata hangs the build.
ENV TZ=Asia/Dubai
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Add repository for g++-10
RUN apt-get update && apt-get install -y software-properties-common
RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test 

# Add repository for latest Cmake
RUN apt-get update && apt-get install -y curl
RUN curl https://apt.kitware.com/keys/kitware-archive-latest.asc | apt-key add -
RUN apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'

# Install dependencies
RUN apt-get update && apt-get install -y g++ gcc cmake git python3-dev g++-10 curl

# Install Rust Nightly
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:$PATH"
RUN echo $PATH
RUN rustup default nightly

# Compile Themisto
RUN git clone https://github.com/algbio/themisto --recursive
WORKDIR /themisto
WORKDIR /themisto/build
RUN cmake .. -DCMAKE_CXX_COMPILER=g++-10 -DMAX_KMER_LENGTH=31 -DCMAKE_BUILD_TYPE=Release -DCMAKE_BUILD_ZLIB=1
RUN make -j8

# Check that the compiled executable runs
RUN /themisto/build/bin/themisto
