all: testutils  

GSL_INC = /home/pd3/cpptools/gsl-1.15/
APP3_CPP_INC = /home/pd3/apparatus3-analysis/
CPP_TOOLS_INC = /home/pd3/cpptools/
CFITSIO_INC = /home/pd3/cpptools/cfitsio/
INC = -I${APP3_CPP_INC} -I${CPP_TOOLS_INC} -I${CFITSIO_INC} -I${GSL_INC}

GSL_LIB = -L/home/pd3/cpptools/gsl-1.15/.libs/ -L/home/pd3/cpptools/gsl-1.15/cblas/.libs/
CCFITS_LIB = /home/pd3/cpptools/CCfits/.libs/libCCfits.so
TIFF_LIB = /home/pd3/cpptools/tiff-4.0.0/libtiff/.libs/libtiff.so

RUN_TIME_PATHS = -R/home/pd3/cpptools/gsl-1.15/.libs/:/home/pd3/cpptools/gsl-1.15/cblas/.libs/:/home/pd3/cpptools/CCfits/.libs/:/home/pd3/cpptools/tiff-4.0.0/libtiff/.libs/

CFLAGS = -Wall ${INC} 
LFLAGS = ${GSL_LIB} -lgsl -lgslcblas -lm ${CCFITS_LIB} ${TIFF_LIB} -Xlinker ${RUN_TIME_PATHS}

objs = testutils.o /home/pd3/apparatus3-analysis/qini/qini_utils.o 

testutils: testutils.o utils.o /home/pd3/apparatus3-analysis/qini/qini_utils.o
	g++ $^ ${LFLAGS} -o testutils

com: com.o utils.o /home/pd3/apparatus3-analysis/qini/qini_utils.o
	g++ $^ ${LFLAGS} -o com
	cp  -v com /lab/software/apparatus3/bin/com

tiff: tiff.o utils.o /home/pd3/apparatus3-analysis/qini/qini_utils.o
	g++ $^ ${LFLAGS} -o tiff
	cp  -v tiff /lab/software/apparatus3/bin/tiff

clean:
	rm -f *.o

.cpp.o:
	indent $<
	g++ $< ${CFLAGS} -c




