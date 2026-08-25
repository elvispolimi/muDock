rm -rf build &&
  cmake -S . -B ./build -DCMAKE_INSTALL_PREFIX=./install -DCMAKE_BUILD_TYPE=Debug -DMUDOCK_ENABLE_TEST=ON &&
  cmake --build ./build
