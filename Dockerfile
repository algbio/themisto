FROM ubuntu:20.04
SHELL ["/bin/bash", "-c"]

# Set some timezone or otherwise tzdata hangs the build.
ENV TZ=Asia/Dubai
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y g++ gcc cmake git python3-dev g++-10 curl

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:$PATH"
RUN echo $PATH
RUN rustup default nightly

RUN git clone https://github.com/iosfwd/themisto --recursive
WORKDIR /themisto
#RUN git checkout dev

WORKDIR /themisto/build
RUN cmake .. -DCMAKE_CXX_COMPILER=g++-10 -DMAX_KMER_LENGTH=31 -DCMAKE_BUILD_TYPE=Release -DCMAKE_BUILD_ZLIB=1 -DCMAKE_BUILD_BZIP2=1
RUN make -j8
RUN /themisto/build/bin/themisto
