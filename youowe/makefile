all : music.js

upload :
	pscp *.js *.png *.wav index.html root@spacebar.org:/var/www/youowe/

music.js : overworld.mid vampires.mid escape.mid boss.mid store.mid makefile
	echo "/* Generated file, do not edit */" > music.js
	../renderhero/dumpmidi overworld.mid song_overworld >> music.js
	../renderhero/dumpmidi vampires.mid song_vampires >> music.js
	../renderhero/dumpmidi escape.mid song_escape >> music.js
	../renderhero/dumpmidi boss.mid song_boss >> music.js
	../renderhero/dumpmidi store.mid song_store >> music.js

