default : otolith.exe tfcompiler.exe

# Note the DLLs in the current directory are 32 bit, so if you change
# the architecture (I couldn't get it to work anyway), then you need to
# change the DLLs too. If the DLLs are the wrong architecture, it will
# immediately exit with no feedback.
# SDLARCH=x64
SDLARCH=x86

# Showing console output: The most important thing is setting the
# subsystem to "console", otherwise, no amount of redirection or
# defines or environment variables will work (I don't think this
# has anything to do with SDL in fact). 

# MLTON=/c/mlton/bin/mlton.bat -target x86_64-w64-mingw32 -codegen amd64
MLTON=/c/mlton/bin/mlton.bat
# -target x86_64-w64-mingw32 -codegen amd64

MINGWGCC=/c/mlton/bin/gcc -std=c99 -O1 -fomit-frame-pointer -fno-strict-aliasing -w -fschedule-insns -fschedule-insns2 -malign-functions=5 -malign-jumps=2 -I SDL/include -I/c/mlton/lib/mlton/include
# x86_64-w64-mingw32-gcc -m64

%.o : %.c makefile *.h
	${MINGWGCC} -c $< -o $@

%.o : ../sdlml/%.c makefile
	${MINGWGCC} -c -D WIN32 $< -o $@

sdl_win32_main.o : sdl_win32_main.c makefile
	${MINGWGCC} -D_GNU_SOURCE=1 -c $< -o $@

OBJECTS=blit.o sdlml.o

world-tf.sml : world.tfdesc tfcompiler.exe
	./tfcompiler.exe world.tfdesc

SMLS=otolith.sml prelude.sml screen-sig.sml screen.sml world-tf.sml draw.sml draw-sig.sml keyed-tesselation.sml keyed-tesselation-sig.sml render-sig.sml render.sml images.sml intmaths-sig.sml intmaths.sml world-sig.sml world.sml physics-sig.sml physics.sml

SDLML=../sdlml/*.sml ../sdlml/*.c

#EXNHISTORY=-const 'Exn.keepHistory true' 
#  -mwindows works, disables console?
otolith.exe : sdl_win32_main.o ${OBJECTS} otolith.mlb ${SMLS} ${SDLML} makefile
	${MLTON} -cc-opt " -I SDL/include -Dmain=SDL_main " -link-opt "-Wl,--subsystem,console -LSDL/lib/${SDLARCH} -lSDL -lkernel32 -lmingw32" -output $@ -default-ann 'allowFFI true' ${EXNHISTORY} otolith.mlb sdl_win32_main.o ${OBJECTS}

tfcompiler.exe : ../sml-lib/textformat/*
	${MLTON} -output $@ ../sml-lib/textformat/tfcompiler.mlb

intmaths-test.exe : intmaths-sig.sml intmaths.sml intmaths-test.sml intmaths-test.mlb
	${MLTON} -output $@ intmaths-test.mlb

clean :
	rm -f *.obj *.exe *.manifest *.ilk *.pdb *.o

veryclean : clean
	rm -f *-tf.sml
