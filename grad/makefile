
default: iter.exe train-mnist.exe optimize.exe chopreduce.exe makesubst.exe colordiff.exe train-chess.exe evaluate-chess.exe train-cifar10.exe hfluint8_test.exe makechop.exe manualchop2.exe histo.exe race.exe eval-table-mnist.exe eval-table-cifar10.exe eval-table-chess.exe zt.exe improvezt.exe

all: grad.exe iter.exe make-data.exe mnist_test.exe train-mnist.exe evaluate-mnist.exe train-chess.exe evaluate-chess.exe smallcrush.exe smallcrush-gen.exe permbench.exe substz3.exe xor.exe goodrandom.exe optimize.exe manualchop.exe manualchop2.exe makesubst.exe chopreduce.exe findchop.exe op_test.exe expression_test.exe makefn.exe hash.exe network-gpu_test.exe network_test.exe hcomplex_test.exe fractal.exe render.exe hash-util_test.exe cifar10_test.exe evaluate-cifar10.exe hfluint8_bench.exe

CC_LIB=../cc-lib
CHESS=../chess

CXX=x86_64-w64-mingw32-g++
CC=x86_64-w64-mingw32-gcc
LFLAGS=-flto -static

# Note, this builds against AMD's distribution of OpenCL, since that's
# the card I had when I started. It also works fine with an NVIDIA GPU,
# and nvidia libs are a bit hard to find since they love CUDA.
AMD_ARCH=x86_64
# AMDSDK="/c/Program Files (x86)/AMD APP SDK"
# mingw wants /c/, but TDM likes c:/
# AMDSDK=c:/Program Files (x86)/AMD APP SDK/3.0
AMDSDK=c:/code/AMD_SDK/

CLINCLUDES="-I$(AMDSDK)/include"
CLLIBS='-L${AMDSDK}/lib/${AMD_ARCH}'

# no fast math!
# this worked for clang
# CXXFLAGS=-mfpmath=sse -msse2 -frounding-math -fno-fast-math -march=native -m64 -Wall -Wno-format -Wno-unused-function -Wno-deprecated -Wno-sign-compare -Wno-mismatched-tags -I. -I$(CC_LIB) -I$(CC_LIB)/re2 -std=c++20
# TODO: This may help if not already enabled?
# -DHALF_ENABLE_F16C_INTRINSICS=1
OPT=-O2
CXXFLAGS=-march=native -m64 -Wall -Wno-format -Wno-unused-function -Wno-deprecated -Wno-sign-compare -Wno-mismatched-tags -Wunsafe-loop-optimizations -I. -I$(CC_LIB) -I$(CC_LIB)/re2 -I$(CHESS) $(CLINCLUDES) -std=c++20 $(OPT)

CCLIB_OBJECTS=$(CC_LIB)/util.o $(CC_LIB)/arcfour.o $(CC_LIB)/base/stringprintf.o $(CC_LIB)/base/logging.o $(CC_LIB)/stb_image.o $(CC_LIB)/stb_image_write.o $(CC_LIB)/color-util.o $(CC_LIB)/image.o $(CC_LIB)/opt/opt.o $(CC_LIB)/bounds.o $(CC_LIB)/bitbuffer.o $(CC_LIB)/bignum/bigz.o $(CC_LIB)/bignum/bign.o $(CC_LIB)/bignum/bigq.o $(CC_LIB)/bignum/big.o $(CC_LIB)/geom/marching.o $(CC_LIB)/edit-distance.o $(CC_LIB)/ansi.o

RE2_OBJECTS=$(CC_LIB)/re2/bitstate.o $(CC_LIB)/re2/compile.o $(CC_LIB)/re2/dfa.o $(CC_LIB)/re2/filtered_re2.o $(CC_LIB)/re2/mimics_pcre.o $(CC_LIB)/re2/nfa.o $(CC_LIB)/re2/onepass.o $(CC_LIB)/re2/parse.o $(CC_LIB)/re2/perl_groups.o $(CC_LIB)/re2/prefilter.o $(CC_LIB)/re2/prefilter_tree.o $(CC_LIB)/re2/prog.o $(CC_LIB)/re2/re2.o $(CC_LIB)/re2/regexp.o $(CC_LIB)/re2/set.o $(CC_LIB)/re2/simplify.o $(CC_LIB)/re2/stringpiece.o $(CC_LIB)/re2/tostring.o $(CC_LIB)/re2/unicode_casefold.o $(CC_LIB)/re2/unicode_groups.o $(CC_LIB)/re2/util/rune.o $(CC_LIB)/re2/util/strutil.o

CHESS_OBJECTS=$(CHESS)/chess.o $(CHESS)/pgn.o

CL_OBJECTS=network.o network-gpu.o error-history.o clutil.o

%.o : %.cc *.h makefile
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@
	@echo -n "."

grad.s : grad.cc makefile
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) -S -fverbose-asm -c $< -o $@

OBJECTS=$(CCLIB_OBJECTS) $(RE2_OBJECTS) $(CHESS_OBJECTS) hash-util.o diff.o hfluint8.o hfluint16.o

# also increase stack size
LFLAGS= -L. -m64 -Wl,--subsystem,console -Wl,--stack,335544320 $(CLLIBS) -lz -lOpenCL $(OPT) $(FLTO) -lpsapi -static

grad.exe : grad.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

iter.exe : iter.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

make-data.exe : make-data.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

network_test.exe : network.o network-test-util.o network_test.o $(CCLIB_OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

network-gpu_test.exe : network.o network-gpu.o network-test-util.o network-gpu_test.o clutil.o $(CCLIB_OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

train-chess.exe : train-chess.o $(CL_OBJECTS) $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

train-mnist.exe : train-mnist.o $(CL_OBJECTS) $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

train-cifar10.exe : train-cifar10.o $(CL_OBJECTS) $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

mnist_test.exe : mnist_test.o $(CCLIB_OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

cifar10_test.exe : cifar10_test.o $(CCLIB_OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

nnchess_test.exe : nnchess_test.o $(CCLIB_OBJECTS) $(CHESS_OBJECTS) $(RE2_OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

evaluate-chess.exe : evaluate-chess.o network.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

evaluate-mnist.exe : evaluate-mnist.o $(CL_OBJECTS) $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

evaluate-cifar10.exe : evaluate-cifar10.o $(CL_OBJECTS) $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

hash.exe : hash.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

makefn.exe : makefn.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

expression_test.exe : expression_test.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

op_test.exe : op_test.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

findchop.exe : findchop.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

chopreduce.exe : chopreduce.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

makechop.exe : makechop.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

makesubst.exe : makesubst.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

manualchop.exe : manualchop.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

manualchop2.exe : manualchop2.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

optimize.exe : optimize.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

improvezt.exe : improvezt.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

goodrandom.exe : goodrandom.o $(OBJECTS) $(CC_LIB)/crypt/cryptrand.o
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

xor.exe : xor.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

substz3.exe : substz3.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

permbench.exe : permbench.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

hcomplex_test.exe : hcomplex_test.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

convert-old.exe : convert-old.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

fractal.exe : fractal.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

render.exe : render.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

hash-util_test.exe : hash-util_test.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

frac3d.exe : frac3d.o meshopt-simplifier.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"


# this does not work. delete
simplify-mesh.exe : simplify-mesh.o
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

colordiff.exe : colordiff.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

hfluint8_test.exe : hfluint8_test.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

hfluint8_bench.exe : hfluint8_bench.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

hfluint16_test.exe : hfluint16_test.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

makeshift.exe : makeshift.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

histo.exe : histo.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

race.exe : race.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

eval-table-mnist.exe : eval-table-mnist.o $(CL_OBJECTS) $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

eval-table-cifar10.exe : eval-table-cifar10.o $(CL_OBJECTS) $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

eval-table-chess.exe : eval-table-chess.o $(CL_OBJECTS) $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

old.exe : old.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

zt.exe : zt.o $(OBJECTS)
	@$(CXX) $^ -o $@ $(LFLAGS)
	@echo -n "!"

# external dependency on TestU01-2009
TESTU01=../../TestU01-2009

testu01.o : testu01.cc
	@$(CXX) $(CXXFLAGS) $^ -o $@ -c -I $(TESTU01)/include

smallcrush-main.o : smallcrush-main.cc
	@$(CXX) $(CXXFLAGS) $^ -o $@ -c -I $(TESTU01)/include

smallcrush-gen.o : smallcrush-gen.cc
	@$(CXX) $(CXXFLAGS) $^ -o $@ -c -I $(TESTU01)/include

smallcrush.exe : smallcrush-main.o $(OBJECTS)
	@$(CXX) $^ -o $@ -L $(TESTU01)/mylib/.libs -L $(TESTU01)/probdist/.libs -L $(TESTU01)/testu01/.libs $(LFLAGS) -ltestu01 -lprobdist -lmylib -lm -lwsock32
	@echo -n "!"

smallcrush-gen.exe : smallcrush-gen.o testu01.o $(OBJECTS)
	@$(CXX) $^ -o $@ -L $(TESTU01)/mylib/.libs -L $(TESTU01)/probdist/.libs -L $(TESTU01)/testu01/.libs $(LFLAGS) -ltestu01 -lprobdist -lmylib -lm -lwsock32
	@echo -n "!"

clean :
	rm -f *.o $(OBJECTS) *.exe

modelclean :
	rm -f grad.val grad.*.val error-*.png error-history.tsv train-*.png modelinfo.png grad-layer*.png
