all: testutils  

GSL_INC = ../gsl-1.15/
CCFITS_INC = ../CCfits/
TIFF_INC = ../tiff-3.9.4/libtiff/

INC = -I../simpleini -I../ -I${CCFITS_INC} -I${GSL_INC} -I${TIFF_INC} 

GSL_LIB = -L../gsl-1.15/.libs/ -L../gsl-1.15/cblas/.libs/
CCFITS_LIB = ../CCfits/.libs/libCCfits.so
TIFF_LIB = ../tiff-3.9.4/libtiff/.libs/libtiff.so

RUN_TIME_PATHS = -R/home/pd3/apparatus3-analysis/gsl-1.15/.libs/:/home/pd3/apparatus3-analysis/gsl-1.15/cblas/.libs/:/home/pd3/apparatus3-analysis/CCfits/.libs/:/home/pd3/apparatus3-analysis/tiff-3.9.4/libtiff/.libs/

CFLAGS = -Wall ${INC} 
LFLAGS = ${GSL_LIB} -lgsl -lgslcblas -lm ${CCFITS_LIB} ${TIFF_LIB} -Xlinker ${RUN_TIME_PATHS}

objs = testutils.o ../qini/qini_utils.o 

testutils: testutils.o utils.o ../qini/qini_utils.o
	g++ $^ ${LFLAGS} -o testutils

com: com.o utils.o ../qini/qini_utils.o
	g++ $^ ${LFLAGS} -o com
	cp  -v com /lab/software/apparatus3/bin/com

tiff: tiff.o utils.o ../qini/qini_utils.o
	g++ $^ ${LFLAGS} -o tiff
	cp  -v tiff /lab/software/apparatus3/bin/tiff

clean:
	rm -f *.o

.cpp.o:
	indent $<
	g++ $< ${CFLAGS} -c




