all: qini

INI_INC = /home/pedro/cpp/simpleini
GSL_INC = /home/pedro/cpp/gsl/include
GSL_LIB = /home/pedro/cpp/gsl/lib

UTILS_INC = /home/pedro/cpp/utils

INC = -I${GSL_INC} -I${INI_INC} -I${UTILS_INC}

CFLAGS = -Wall ${INC}
LFLAGS = -lgsl -lgslcblas -lm

objs = qini.o qini_utils.o

qini: ${objs}
	g++  -L${GSL_LIB} ${objs} ${LFLAGS} -o qini 
	cp -v qini /lab/software/apparatus3/bin/qini
clean:
	rm -f *.o

.cpp.o:
	indent $<
	g++ ${CFLAGS}  $< -c




