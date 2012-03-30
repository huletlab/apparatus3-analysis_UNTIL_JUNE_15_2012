all: testfits

GSL_INC = /home/pd3/cpptools/gsl-1.15/
APP3_CPP_INC = /home/pd3/apparatus3-analysis/
CPP_TOOLS_INC = /home/pd3/cpptools/
INC = -I${APP3_CPP_INC} -I${CPP_TOOLS_INC} -I${GSL_INC}

GSL_LIB = -L/home/pd3/cpptools/gsl-1.15/.libs/ -L/home/pd3/cpptools/gsl-1.15/cblas/.libs/
CCFITS_LIB = /home/pd3/cpptools/CCfits/.libs/libCCfits.so
TIFF_LIB = /home/pd3/cpptools/tiff-4.0.0/libtiff/.libs/libtiff.so

RUN_TIME_PATHS = -R/home/pd3/cpptools/gsl-1.15/.libs/:/home/pd3/cpptools/gsl-1.15/cblas/.libs/:/home/pd3/cpptools/CCfits/.libs/:/home/pd3/cpptools/tiff-4.0.0/libtiff/.libs/

CFLAGS = -Wall ${INC} -fopenmp 
LFLAGS = ${GSL_LIB} -lgsl -lgslcblas -lm ${CCFITS_LIB} ${TIFF_LIB} -Xlinker ${RUN_TIME_PATHS}

objs =   ../qini/qini_utils.o ../utils/utils.o ../funcs/funcs.o 

testfits: testfits.o fits.o
	g++  $^ ${objs} ${LFLAGS} -o testfits

clean:
	rm -f *.o

.cpp.o:
	indent $<
	g++ ${CFLAGS} $< -c




