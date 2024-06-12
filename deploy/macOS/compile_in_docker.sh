#!/bin/sh
## Wrapper for calling the build script `build.sh` inside the
## macoss-cross-compiler docker container. 
#
# Arguments
## 1: version number to build (checked out from the git source tree)
## 2: architecture (one of x86-64,arm64)

set -eo pipefail

VER=$1
if [[ -z $VER ]]; then
  echo "Error: specify version as argument 1"
  exit;
fi

ARCH=$2
if [[ -z $ARCH ]]; then
  echo "Error: specify architecture (one of x86-64,arm64) as argument 2"
  exit;
fi

set -ux

docker run \
  -v `pwd`:/io \
  --rm \
  -it \
  ghcr.io/shepherdjerred/macos-cross-compiler:latest \
  /bin/bash /io/build.sh $1 $2

