all: qini

APP3_CPP_INC = /lab/software/apparatus3/cpp/

INC = -I${APP3_CPP_INC} 

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




