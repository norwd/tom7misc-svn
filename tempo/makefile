all : tempo.exe am2315_test.exe

# Ideally we would rebuild .o files if SVN_REVISION changes; is that possible?
# (not last-changed-revision, since we want to include a change to
# e.g. cc-lib)

CXX=g++
CXXFLAGS=-march=native -Wall -Wno-sign-compare -Wno-deprecated-declarations -O --std=c++17 -I ../cc-lib `mysql_config --cflags` -DSVN_REVISION=`svn info --show-item revision` -I/usr/include/mysql++ -Wno-psabi -ffunction-sections
LFLAGS=`mysql_config --libs` -lmysqlpp -lpthread

%.o : %.cc *.h
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
	@bash -c "echo -n '.'"

CCLIB_OBJECTS=../cc-lib/util.o ../cc-lib/webserver.o ../cc-lib/base/logging.o ../cc-lib/base/stringprintf.o ../cc-lib/stb_image_write.o ../cc-lib/stb_image.o ../cc-lib/image.o ../cc-lib/color-util.o ../cc-lib/pi/netutil.o ../cc-lib/pi/bcm2835.o ../cc-lib/arcfour.o ../cc-lib/process-util.o ../cc-lib/pi/pi-util.o

OBJECTS=tempo-util.o preserver.o server.o database.o onewire.o am2315.o logic.o

tempo.exe : tempo.o $(OBJECTS) $(CCLIB_OBJECTS)
	$(CXX) $^ -o $@ $(LFLAGS)

am2315_test.exe : am2315_test.o am2315.o $(CCLIB_OBJECTS)
	$(CXX) $^ -o $@ $(LFLAGS)

backfill.exe : backfill.cc
	x86_64-w64-mingw32-g++.exe -O backfill.cc -o backfill.exe

clean :
	rm -f *.o *.exe core $(CCLIB_OBJECTS)

veryclean : clean
	rm -f *~
