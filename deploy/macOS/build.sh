#!/bin/bash
## Build script for cross-compiling Themisto for macOS x86-64 or arm64.
## Call this from `compile_in_docker.sh` unless you know what you're doing.

set -exo pipefail

VER=$1
if [[ -z $VER ]]; then
  echo "Error: specify version"
  exit;
fi

ARCH=$2
if [[ -z $ARCH ]]; then
  echo "Error: specify architecture (one of x86-64,arm64)"
  exit;
fi

apt update
apt install -y cmake git libomp5 libomp-dev curl

rustup default stable

mkdir /io/tmp
cd /io/tmp

# Extract and enter source
git clone https://github.com/algbio/Themisto.git
cd Themisto
git checkout ${VER}
git submodule update --init --recursive

mkdir -p ggcat/.cargo

sed 's/cargo build/cargo build/g' ggcat/crates/capi/ggcat-cpp-api/Makefile > Makefile.tmp
mv Makefile.tmp ggcat/crates/capi/ggcat-cpp-api/Makefile

cd build
target_arch=""
if [ "$ARCH" = "x86-64" ]; then
    # Rust toolchain
    rustup target add x86_64-apple-darwin

    echo "[build]" >> ../ggcat/.cargo/config.toml
    echo "target = \"x86_64-apple-darwin\"" >> ../ggcat/.cargo/config.toml
    echo "[target.x86_64-apple-darwin]" >> ../ggcat/.cargo/config.toml
    echo "linker = \"x86_64-apple-darwin22-gcc\"" >> ../ggcat/.cargo/config.toml

    export CC="x86_64-apple-darwin22-gcc"
    export CXX="x86_64-apple-darwin22-g++"

    ## Setup ggcat-cpp-api cargo config files for cross compilation
    sed "s/cargo build/RUSTFLAGS='-L \/osxcross\/SDK\/MacOSX13.0.sdk\/usr\/lib' cargo build --target x86_64-apple-darwin/g" ../ggcat/crates/capi/ggcat-cpp-api/Makefile | sed 's/target\/release/target\/x86_64-apple-darwin\/release/g' | sed 's/fPIE/fPIE -march=x86-64 -mtune=generic -m64 -fPIC -static-libstdc++ -static-libgcc/g' | sed 's/ar cr/\/gcc\/bin\/x86_64-apple-darwin22-gcc-ar cr/g' > Makefile.tmp
    mv Makefile.tmp ../ggcat/crates/capi/ggcat-cpp-api/Makefile

    ## Prevent KMC from linking statically
    sed 's/-static-libgcc//g' ../SBWT/KMC/CMakeLists.txt | sed 's/-static-libstdc++//g' | sed 's/-static//g' > CMakeLists.txt.tmp
    mv CMakeLists.txt.tmp ../SBWT/KMC/CMakeLists.txt

    ## Prevent sdsl-lite from building with native instructions
    sed 's/-march=native/-march=x86-64/g' ../SBWT/sdsl-lite/CMakeLists.txt > CMakeLists.txt.tmp
    mv CMakeLists.txt.tmp ../SBWT/sdsl-lite/CMakeLists.txt

    ## Setup ggcat-cpp-api cargo config files for cross compilation
    cat ../ggcat/crates/capi/ggcat-cpp-api/Makefile | sed 's/fPIE/fPIE -march=x86-64 -mtune=generic -m64 -fPIC/g' > Makefile.tmp
    mv Makefile.tmp ../ggcat/crates/capi/ggcat-cpp-api/Makefile

    # compile x86_64
    cmake -DCMAKE_TOOLCHAIN_FILE="/io/$ARCH-toolchain_GNU.cmake" \
          -DCMAKE_C_FLAGS="-march=$ARCH -mtune=generic -m64 -fPIC -fPIE -static-libstdc++ -static-libgcc" \
          -DCMAKE_CXX_FLAGS="-march=$ARCH -mtune=generic -m64 -fPIC -fPIE -static-libstdc++ -static-libgcc" \
	  -DROARING_DISABLE_NATIVE=ON \
          -DBZIP2_LIBRARIES="/osxcross/SDK/MacOSX13.0.sdk/usr/lib/libbz2.tbd" -DBZIP2_INCLUDE_DIR="/osxcross/SDK/MacOSX13.0.sdk/usr/include" \
          -DZLIB="/osxcross/SDK/MacOSX13.0.sdk/usr/lib/libz.tbd" -DZLIB_INCLUDE_DIR="/osxcross/SDK/MacOSX13.0.sdk/usr/include" \
	  -DZLIB_LIBRARY=="/osxcross/SDK/MacOSX13.0.sdk/usr/lib/libz.tbd" \
	  -DMAX_KMER_LENGTH=31 \
	  ..
    target_arch="x86_64-apple-darwin22"
elif [ "$ARCH" = "arm64" ]; then
    # Rust toolchain
    rustup target add aarch64-apple-darwin

    echo "[build]" >> ../ggcat/.cargo/config.toml
    echo "target = \"aarch64-apple-darwin\"" >> ../ggcat/.cargo/config.toml
    echo "[target.aarch64-apple-darwin]" >> ../ggcat/.cargo/config.toml
    echo "linker = \"aarch64-apple-darwin22-gcc\"" >> ../ggcat/.cargo/config.toml

    export CC="aarch64-apple-darwin22-gcc"
    export CXX="aarch64-apple-darwin22-g++"

    ## Prevent KMC from compiling with x86 instructions
    sed 's/^UNAME_S.*$/UNAME_S=Darwin/g' ../SBWT/KMC/Makefile | sed 's/^UNAME_M.*$/UNAME_M=arm64/g' | sed 's/^[[:space:]]*CC[[:space:]]*=.*$//g' > Makefile.tmp
    sed 's/if(APPLE)/if(TRUE)/g' ../SBWT/KMC/CMakeLists.txt | sed 's/CMAKE_SYSTEM_PROCESSOR MATCHES "arm64"/TRUE/g' > CMakeLists.txt.tmp
    mv Makefile.tmp ../SBWT/KMC/Makefile
    mv CMakeLists.txt.tmp ../SBWT/KMC/CMakeLists.txt

    ## Prevent KMC from linking statically
    sed 's/-static-libgcc//g' ../SBWT/KMC/CMakeLists.txt | sed 's/-static-libstdc++//g' | sed 's/-static//g' > CMakeLists.txt.tmp
    mv CMakeLists.txt.tmp ../SBWT/KMC/CMakeLists.txt

    ## Prevent sdsl-lite from building with native instructions
    sed 's/-msse4.2[[:space:]]*-march=native/-march=armv8-a/g' ../SBWT/sdsl-lite/CMakeLists.txt > CMakeLists.txt.tmp
    mv CMakeLists.txt.tmp ../SBWT/sdsl-lite/CMakeLists.txt

    ## Setup ggcat-cpp-api cargo config files for cross compilation
    sed "s/cargo build/RUSTFLAGS='-L \/osxcross\/SDK\/MacOSX13.0.sdk\/usr\/lib' cargo build --target aarch64-apple-darwin/g" ../ggcat/crates/capi/ggcat-cpp-api/Makefile | sed 's/target\/release/target\/aarch64-apple-darwin\/release/g' | sed 's/fPIE/fPIE -march=armv8-a -mtune=generic -m64 -fPIC -static-libstdc++ -static-libgcc/g' | sed 's/ar cr/\/gcc\/bin\/aarch64-apple-darwin22-gcc-ar cr/g' > Makefile.tmp
    mv Makefile.tmp ../ggcat/crates/capi/ggcat-cpp-api/Makefile

    # compile aarch64
    cmake -DCMAKE_TOOLCHAIN_FILE="/io/$ARCH-toolchain_GNU.cmake" \
          -DCMAKE_C_FLAGS="-march=armv8-a -mtune=generic -m64 -fPIC -fPIE -static-libstdc++ -static-libgcc" \
          -DCMAKE_CXX_FLAGS="-march=armv8-a -mtune=generic -m64 -fPIC -fPIE -static-libstdc++ -static-libgcc" \
	  -DAPPLE=1 \
	  -DCMAKE_SYSTEM_PROCESSOR="arm64" \
	  -DROARING_DISABLE_NATIVE=ON \
          -DBZIP2_LIBRARIES="/osxcross/SDK/MacOSX13.0.sdk/usr/lib/libbz2.tbd" -DBZIP2_INCLUDE_DIR="/osxcross/SDK/MacOSX13.0.sdk/usr/include" \
          -DZLIB="/osxcross/SDK/MacOSX13.0.sdk/usr/lib/libz.tbd" -DZLIB_INCLUDE_DIR="/osxcross/SDK/MacOSX13.0.sdk/usr/include" \
	  -DZLIB_LIBRARY=="/osxcross/SDK/MacOSX13.0.sdk/usr/lib/libz.tbd" \
	  -DMAX_KMER_LENGTH=31 \
	  ..
    target_arch="aarch64-apple-darwin22"
fi

make VERBOSE=1 -j

## gather the stuff to distribute
target=themisto-${VER}-$target_arch
path=/io/tmp/$target
mkdir $path
cp ../build/bin/themisto $path/
cp ../README.md $path/
cp ../LICENSE.txt $path/
cd /io/tmp
tar -zcvf $target.tar.gz $target
mv $target.tar.gz /io/
cd /io/
rm -rf tmp cache

