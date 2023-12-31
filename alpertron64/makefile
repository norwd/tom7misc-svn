default: quad.exe test.exe montgomery64_test.exe

UNAME := $(shell uname -o)
ifeq ($(UNAME),Cygwin)
  # for 64 bits on windows (mingw).
  # For locally compiled libgmp; pick the right directory:
  LIBGMP=../../gmp-6.3.0

  GMP_LFLAGS=$(LIBGMP)/.libs/libgmp.a
  GMP_CFLAGS=-I$(LIBGMP)/ -DBIG_USE_GMP

  CXX=x86_64-w64-mingw32-g++
  CC=x86_64-w64-mingw32-g++
  PLATFORM_CFLAGS=
  PLATFORM_LFLAGS=-static
else ifeq ($(UNAME),Darwin)
  # Mac OS X
  CXX=g++
  CC=gcc
  GMP_LFLAGS=-lgmp
  GMP_CFLAGS=-DBIG_USE_GMP
else
  # Linux
  CXX=g++
  CC=gcc
  # Note! asan is enabled on linux, since I'm usually just using linux
  # for improved debugging tools.
  #
  # I turned off format string warnings on linux, since it complains about
  # %lld and I can't figure out any way to disable just that warning. :(
  # PLATFORM_CFLAGS=-fsanitize=undefined -fsanitize=address -Wno-format
  # For profiling
  PLATFORM_CFLAGS=-pg -Wno-format -DRUN_LOOP=1
  GMP_LFLAGS=-lgmp
  GMP_CFLAGS=-DBIG_USE_GMP
endif

CC_LIB=../cc-lib
SOS=../sos

# Add this to see dead code: -ffunction-sections -fdata-sections -Wl,--gc-sections -Wl,--print-gc-sections
CXXFLAGS=-std=c++20 -I $(SOS) -I $(CC_LIB) $(GMP_CFLAGS) -g -O2 -march=native -Wall $(PLATFORM_CFLAGS)
LFLAGS= $(GMP_LFLAGS) $(PLATFORM_LFLAGS)

HEADERS=bignbr.h modmult.h quadmodll.h bigconv.h quad.h montgomery64.h

# XXX util temporary
CC_LIB_OBJECTS=$(CC_LIB)/base/logging.o $(CC_LIB)/base/stringprintf.o $(CC_LIB)/base/int128.o $(CC_LIB)/bignum/big.o $(CC_LIB)/factorization.o $(CC_LIB)/util.o

TEST_CC_LIB_OBJECTS=$(CC_LIB_OBJECTS) $(CC_LIB)/ansi.o $(CC_LIB)/util.o $(CC_LIB)/arcfour.o $(CC_LIB)/subprocess.o $(CC_LIB)/atomic-util.o

SOS_OBJECTS=$(SOS)/sos-util.o

ALPERTRON_OBJECTS=modmult.o quad.o quadmodll.o bignbr.o bigconv.o

quad.exe:  main.o $(ALPERTRON_OBJECTS) $(CC_LIB_OBJECTS)
	@echo -n "@"
	@$(CXX) $^ $(CXXFLAGS) $(LFLAGS) -lm -o $@

%.o : %.cc $(HEADERS)
	@echo -n "."
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

grid.exe : grid.o $(ALPERTRON_OBJECTS) $(SOS_OBJECTS) $(TEST_CC_LIB_OBJECTS)
	@echo -n "@"
	@$(CXX) $^ $(CXXFLAGS) $(LFLAGS) -lm -o $@

maxvalues.exe : maxvalues.o $(ALPERTRON_OBJECTS) $(SOS_OBJECTS) $(TEST_CC_LIB_OBJECTS)
	@echo -n "@"
	@$(CXX) $^ $(CXXFLAGS) $(LFLAGS) -lm -o $@

montgomery64_test.exe : montgomery64_test.o $(TEST_CC_LIB_OBJECTS)
	@echo -n "@"
	@$(CXX) $^ $(CXXFLAGS) $(LFLAGS) -lm -o $@

test.exe : test.o $(ALPERTRON_OBJECTS) $(TEST_CC_LIB_OBJECTS)
	@echo -n "@"
	@$(CXX) $^ $(CXXFLAGS) $(LFLAGS) -lm -o $@

clean:
	rm -f *.o $(TEST_CC_LIB_OBJECTS) coverage.out *.exe


# This file is part of Alpertron Calculators.
#
# Copyright 2018 Dario Alejandro Alpern
#
# Alpertron Calculators is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# Alpertron Calculators is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Alpertron Calculators. If not, see
# <http://www.gnu.org/licenses/>.
