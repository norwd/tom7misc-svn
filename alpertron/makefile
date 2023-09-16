default: quad.exe afactor.exe test.exe

UNAME := $(shell uname -o)
ifeq ($(UNAME),Cygwin)
  # for 64 bits on windows (mingw).
  # For locally compiled libgmp; pick the right directory:
  LIBGMP=../../gmp-6.3.0

  GMP_LFLAGS=$(LIBGMP)/.libs/libgmp.a
  GMP_CFLAGS=-I$(LIBGMP)/ -DBIG_USE_GMP

  CXX=x86_64-w64-mingw32-g++
  CC=x86_64-w64-mingw32-g++
  PLATFORM_CFLAGS=-Wl,--stack,167772160
  PLATFORM_LFLAGS=-static
else
  # Linux
  CXX=g++
  CC=gcc
  PLATFORM_CFLAGS=-Wl,-z,stack-size=167772160
  GMP_LFLAGS=-lgmp
  GMP_CFLAGS=-DBIG_USE_GMP
endif

CC_LIB=../cc-lib

# PERF: gotta reduce this stack usage!
# XXX: Opt flags, -g
# Add this to see dead code: -ffunction-sections -fdata-sections -Wl,--gc-sections -Wl,--print-gc-sections
CXXFLAGS=-I $(CC_LIB) $(GMP_CFLAGS) -g -O -march=native -Wall $(PLATFORM_CFLAGS)
LFLAGS= $(GMP_LFLAGS) $(PLATFORM_LFLAGS)

HEADERS=bignbr.h enums.h factor.h globals.h modmult.h quadmodll.h baseconv.h bigconv.h quad.h

CC_LIB_OBJECTS=$(CC_LIB)/base/logging.o $(CC_LIB)/base/stringprintf.o $(CC_LIB)/bignum/big.o

ALPERTRON_OBJECTS=copystr.o bigint.o baseconv.o modmult.o quad.o quadmodll.o bignbr.o globals.o bigconv.o bigfactor.o

quad.exe:  main.o $(ALPERTRON_OBJECTS) $(CC_LIB_OBJECTS)
	@echo -n "@"
	@$(CXX) $^ $(CXXFLAGS) $(LFLAGS) -lm -o $@

afactor.exe:  afactor.o $(ALPERTRON_OBJECTS) $(CC_LIB_OBJECTS)
	@echo -n "@"
	@$(CXX) $^ $(CXXFLAGS) $(LFLAGS) -lm -o $@

test.exe : test.o $(ALPERTRON_OBJECTS) $(CC_LIB_OBJECTS)
	@echo -n "@"
	@$(CXX) $^ $(CXXFLAGS) $(LFLAGS) -lm -o $@

%.o : %.cc $(HEADERS)
	@echo -n "."
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -f *.o coverage.out *.exe


#    This file is part of Alpertron Calculators.
#
#    Copyright 2018 Dario Alejandro Alpern
#
#    Alpertron Calculators is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Alpertron Calculators is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Alpertron Calculators.  If not, see <http://www.gnu.org/licenses/>.