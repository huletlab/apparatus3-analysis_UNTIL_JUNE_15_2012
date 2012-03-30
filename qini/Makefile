all: qini

APP3_CPP_INC = /home/pd3/apparatus3-analysis/
CPP_TOOLS_INC = /home/pd3/cpptools/

INC = -I${APP3_CPP_INC} -I${CPP_TOOLS_INC}

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




