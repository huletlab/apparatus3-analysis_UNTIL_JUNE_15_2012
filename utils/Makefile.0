all: testutils  

GSL_INC = APP3-CPP/gsl-1.15/
APP3_CPP_INC = APP3-CPP/
INC = -I${APP3_CPP_INC} -I${GSL_INC}

GSL_LIB = -LAPP3-CPP/gsl-1.15/.libs/ -LAPP3-CPP/gsl-1.15/cblas/.libs/
CCFITS_LIB = APP3-CPP/CCfits/.libs/libCCfits.so
TIFF_LIB = APP3-CPP/tiff-3.9.4/libtiff/.libs/libtiff.so

RUN_TIME_PATHS = -RAPP3-CPP/gsl-1.15/.libs/:APP3-CPP/gsl-1.15/cblas/.libs/:APP3-CPP/CCfits/.libs/:APP3-CPP/tiff-3.9.4/libtiff/.libs/

CFLAGS = -Wall ${INC} 
LFLAGS = ${GSL_LIB} -lgsl -lgslcblas -lm ${CCFITS_LIB} ${TIFF_LIB} -Xlinker ${RUN_TIME_PATHS}

objs = testutils.o APP3-CPP/qini/qini_utils.o 

testutils: testutils.o utils.o APP3-CPP/qini/qini_utils.o
	g++ $^ ${LFLAGS} -o testutils

com: com.o utils.o APP3-CPP/qini/qini_utils.o
	g++ $^ ${LFLAGS} -o com
	cp  -v com /lab/software/apparatus3/bin/com

tiff: tiff.o utils.o APP3-CPP/qini/qini_utils.o
	g++ $^ ${LFLAGS} -o tiff
	cp  -v tiff /lab/software/apparatus3/bin/tiff

clean:
	rm -f *.o

.cpp.o:
	indent $<
	g++ $< ${CFLAGS} -c




