docker run -t -i --rm \
  -v `pwd`:/io \
  -e CC='ccache gcc' \
  -e CXX='ccache g++' \
  -e CCACHE_DIR='/io/cache' \
  phusion/holy-build-box-64:3.0.2 \
  bash /io/build.sh $1
