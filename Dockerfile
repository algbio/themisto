FROM ubuntu:20.04

# Set some timezone or otherwise tzdata hangs the build.
ENV TZ=Asia/Dubai
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y g++ gcc cmake git python3-dev g++-10 libz-dev

RUN git clone https://github.com/algbio/Themisto --recursive
WORKDIR /Themisto
RUN git checkout dev

WORKDIR /Themisto/build
#RUN cmake .. -DCMAKE_CXX_COMPILER=g++-8 -DMAX_KMER_LENGTH=32 -DCMAKE_BUILD_TYPE=Release -DCMAKE_BUILD_ZLIB=1 -DCMAKE_BUILD_BZIP2=1
RUN cmake .. -DCMAKE_CXX_COMPILER=g++-10 -DMAX_KMER_LENGTH=255 -DCMAKE_BUILD_TYPE=Debug -DCMAKE_BUILD_ZLIB=1 # Debug build
RUN make -j8
run /Themisto/build/bin/themisto
