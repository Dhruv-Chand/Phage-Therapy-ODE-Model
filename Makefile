phageTherapy: phageTherapy.o
	c++ -o phageTherapy phageTherapy.o -ltrapfpe -lpgplot -lcpgplot -lX11 -lm 

phageTherapy.o: phageTherapy.cpp
	c++ -c phageTherapy.cpp

junk: junk.o
	c++ -o junk junk.o -ltrapfpe -lpgplot -lcpgplot -lX11 -lm 

junk.o: junk.cpp
	c++ -c junk.cpp

junk2: junk2.o
	c++ -o junk2 junk2.o -ltrapfpe -lpgplot -lcpgplot -lX11 -lm 

junk2.o: junk2.cpp
	c++ -c junk2.cpp