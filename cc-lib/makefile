
default : heap_test.exe minmax-heap_test.exe rle_test.exe interval-tree_test.exe threadutil_test.exe color-util_test.exe textsvg_test.exe lines_test.exe image_test.exe util_test.exe randutil_test.exe json_test.exe arcfour_test.exe lastn-buffer_test.exe list-util_test.exe edit-distance_test.exe re2_test.exe webserver_test.exe stb_image_bench.exe process-util_test.exe stb_truetype_test.exe packrect_test.exe xml_test.exe mp3_test.exe top_test.exe bounds_test.exe optional-iterator_test.exe tuple-util_test.exe bitbuffer_test.exe ansi_test.exe factorization_test.exe atomic-util_test.exe interval-cover_test.exe set-util_test.exe hashing_test.exe subprocess_test.exe city_test.exe autoparallel_test.exe int128_test.exe do-not-optimize_test.exe $(TESTCOMPILE)

TESTCOMPILE=stb_image_write.o stb_image.o dr_wav.o bounds.o simpledxf.o bitbuffer.o

BASE=base/logging.o base/stringprintf.o

UNAME := $(shell uname -o)
ifeq ($(UNAME),Cygwin)
  # for 64 bits on windows (mingw).
  CXX=x86_64-w64-mingw32-g++
  CC=x86_64-w64-mingw32-g++
  CXXFLAGS=-march=native -pipe -Wall -I. --std=c++20 -O3 -static -flto=auto
  SOCK_LFLAGS= -lws2_32
else
  # Linux
  CXX=g++
  CC=gcc
  # no-psabi suppresses warnings all over STL in GCC 7.1 for ARM (e.g. pi)
  # it can probably be removed in a few years?
  CXXFLAGS=-march=native -Wall -I. --std=c++20 -O3 -Wno-psabi
  SOCK_LFLAGS=
endif

RE2_OBJECTS=re2/bitstate.o re2/compile.o re2/dfa.o re2/filtered_re2.o re2/mimics_pcre.o re2/nfa.o re2/onepass.o re2/parse.o re2/perl_groups.o re2/prefilter.o re2/prefilter_tree.o re2/prog.o re2/re2.o re2/regexp.o re2/set.o re2/simplify.o re2/stringpiece.o re2/tostring.o re2/unicode_casefold.o re2/unicode_groups.o re2/util/rune.o re2/util/strutil.o

%.o : %.cc *.h makefile
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@
	@echo -n "."

ansi_test.exe : ansi_test.o ansi.o util.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

stb_truetype_test.exe : stb_truetype.o stb_truetype_test.o
	$(CXX) $(CXXFLAGS) $^ -o $@

stb_image_bench.o : stb_image_bench.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

stb_image_bench.exe : stb_image.o stb_image_write.o stb_image_bench.o image.o arcfour.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

heap_test.o : heap_test.cc heap.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

heap_test.exe : heap_test.o
	$(CXX) $(CXXFLAGS) $^ -o $@

minmax-heap_test.o : minmax-heap_test.cc minmax-heap.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

minmax-heap_test.exe : minmax-heap_test.o arcfour.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

list-util_test.o : list-util_test.cc list-util.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

list-util_test.exe : list-util_test.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

tuple-util_test.o : tuple-util_test.cc tuple-util.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

tuple-util_test.exe : tuple-util_test.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

rle_test.o : rle_test.cc rle.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

rle_test.exe : rle_test.o rle.o $(BASE) arcfour.o
	$(CXX) $(CXXFLAGS) $^ -o $@

arcfour_test.exe : arcfour_test.o arcfour.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

arcfour_bench.exe : arcfour_bench.o arcfour.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

lastn-buffer_test.o : lastn-buffer_test.cc lastn-buffer.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

lastn-buffer_test.exe : lastn-buffer_test.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

interval-tree_test.exe : interval-tree_test.o $(BASE) arcfour.o
	$(CXX) $(CXXFLAGS) $^ -o $@

simpledxf.o : simpledxf.cc simpledxf.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

threadutil_test.exe : threadutil.h threadutil_test.o ansi.o util.o $(BASE)
	$(CXX) $(CXXFLAGS) threadutil_test.o ansi.o util.o $(BASE) -o $@ -lpthread

color-util_test.exe : color-util.o color-util_test.o stb_image.o image.o stb_image_write.o arcfour.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

textsvg_test.exe : textsvg.o textsvg_test.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

lines_test.exe : lines_test.o arcfour.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

image_test.exe : image_test.o arcfour.o image.o stb_image.o stb_image_write.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

util_test.exe : util_test.o util.o arcfour.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

randutil_test.exe : randutil.h randutil_test.o arcfour.o $(BASE)
	$(CXX) $(CXXFLAGS) randutil_test.o arcfour.o $(BASE) -o $@

json_test.exe : json_test.o $(BASE)
	$(CXX) $(CXXFLAGS) json_test.o $(BASE) -o $@

edit-distance_test.exe : edit-distance_test.o edit-distance.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

re2_test.exe : re2_test.o $(RE2_OBJECTS) $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@ -lpthread

webserver_test.exe : webserver_test.o webserver.o util.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(SOCK_LFLAGS) -lpthread

process-util_test.exe : process-util_test.o process-util.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

packrect_test.exe : packrect_test.o packrect.o opt/opt.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

xml_test.exe : xml_test.o xml.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

mp3_test.exe : mp3_test.o mp3.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

top_test.exe : top.o top_test.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

bounds_test.exe : bounds.o bounds_test.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

optional-iterator_test.exe : optional-iterator_test.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

periodically_test.exe : periodically_test.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

image_bench.exe : image_bench.o image.o stb_image.o stb_image_write.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

bitbuffer_test.exe : bitbuffer.o bitbuffer_test.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

image-frgba_test.exe : image-frgba.o image-frgba_test.o tinyexr.o image.o stb_image.o stb_image_write.o color-util.o util.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

vector-util_test.exe : vector-util_test.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

atomic-util_test.exe : atomic-util_test.o atomic-util.o ansi.o util.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

interval-cover_test.exe : interval-cover.o interval-cover_test.o ansi.o util.o arcfour.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

factorization_test.exe : factorization.o factorization_test.o ansi.o util.o arcfour.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

set-util_test.exe : set-util_test.o ansi.o util.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

hashing_test.exe : hashing_test.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

city_test.exe : city_test.o city/city.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

autoparallel_test.exe : autoparallel_test.o arcfour.o util.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

divbench.exe : divbench.o ansi.o util.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

int128_test.exe : int128_test.o base/int128.o ansi.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Note: Only implemented on windows.
subprocess_test.exe : subprocess_test.o subprocess.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

do-not-optimize_test.exe : do-not-optimize_test.o util.o $(BASE)
	$(CXX) $(CXXFLAGS) $^ -o $@

clean :
	rm -f *.o *.exe *.deleteme test-*.png $(BASE) $(RE2_OBJECTS)
