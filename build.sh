#./clean.sh

cmake -S . -B build
cmake --build build -- -j 4 