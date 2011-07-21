all: testutils  

INI_INC = /lab/software/apparatus3/cpp/qini
CCFITS_INC = /usr/local/include/CCfits
GSL_INC = /lab/software/apparatus3/cpp/gsl/include

INC = -I${INI_INC} -I${CCFITS_INC} -I${GSL_INC} 

GSL_LIB = /lab/software/apparatus3/cpp/gsl/lib

CFLAGS = -Wall ${INC} 
LFLAGS = -L${GSL_LIB} -lgsl -lgslcblas -lm -lCCfits -ltiff

objs = testutils.o /lab/software/apparatus3/cpp/qini/qini_utils.o 

testutils: /lab/software/apparatus3/cpp/utils/testutils.o /lab/software/apparatus3/cpp/utils/utils.o /lab/software/apparatus3/cpp/qini/qini_utils.o
	g++  $^ ${LFLAGS} -o testutils

com: com.o utils.o /lab/software/apparatus3/cpp/qini/qini_utils.o
	g++ $^ ${LFLAGS} -o com
	cp  -v com /lab/software/apparatus3/bin/com

tiff: tiff.o utils.o /lab/software/apparatus3/cpp/qini/qini_utils.o
	g++ $^ ${LFLAGS} -o tiff
	cp  -v tiff /lab/software/apparatus3/bin/tiff

clean:
	rm -f *.o

.cpp.o:
	indent $<
	g++ $< ${CFLAGS} -c




