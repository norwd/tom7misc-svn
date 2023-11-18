#!/bin/bash

make clean
make -j quad.exe || exit -1

set -x
set -e

rm -f coverage.out
./quad.exe 13526043882269 >> coverage.out
./quad.exe 100 >> coverage.out
./quad.exe 100000000000 >> coverage.out

# a factor of 3^1 does trigger a different code path
./quad.exe 120 >> coverage.out
./quad.exe 62768369664000 >> coverage.out

./quad.exe 199506591167822449 >> coverage.out

# 2^62
./quad.exe 4611686018427387904 >> coverage.out
# 2^61
./quad.exe 2305843009213693952 >> coverage.out
# 16!
./quad.exe 20922789888000 >> coverage.out

# 2^31-1 * 65537
./quad.exe 140739635773439 >> coverage.out

# (2^31-1)^2
./quad.exe 4611686014132420609 >> coverage.out

# Test prime factors > 2^32; overflow becomes a risk.

# (2^33 - 303)
# x=49367, y=78440
./quad.exe 8589934289 >> coverage.out

# (2^33 - 303) * 31337
# two solutions
./quad.exe 269182770814393 >> coverage.out

dos2unix -q coverage.out
diff coverage.golden coverage.out
exit_status=$?
if [ $exit_status -eq 0 ]; then
    echo OK
else
    echo DIFFERENCES
fi