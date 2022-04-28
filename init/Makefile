CC = g++

run: main.o empfuncs.o Makefile
	$(CC) -o run main.o empfuncs.o

main.o: main.cpp Makefile
	$(CC) -c main.cpp

empfuncs.o: empfuncs.cpp Makefile
	$(CC) -c empfuncs.cpp

clean:
	rm -f *~

distclean:
	make clean
	rm -f *.o run
