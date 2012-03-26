all: qini

INI_INC = ../simpleini

INC = -I${INI_INC} 

CFLAGS = -Wall ${INC}
LFLAGS = -lm

objs = qini.o qini_utils.o

qini: ${objs}
	g++  ${objs} ${LFLAGS} -o qini 
	cp -v qini ../bin/
clean:
	rm -f *.o

.cpp.o:
	indent $<
	g++ ${CFLAGS}  $< -c




