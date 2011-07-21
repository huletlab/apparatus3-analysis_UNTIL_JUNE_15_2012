all: testfourier

INI_INC = /lab/software/apparatus3/cpp/qini
CCFITS_INC = /usr/local/include/CCfits
GSL_INC = /lab/software/apparatus3/cpp/gsl/include

INC = -I${INI_INC} -I${CCFITS_INC} -I${GSL_INC} 

GSL_LIB = /lab/software/apparatus3/cpp/gsl/lib

CFLAGS = -Wall ${INC} 
LFLAGS = -L${GSL_LIB} -lgsl -lgslcblas -lm -lCCfits -ltiff


objs =   /lab/software/apparatus3/cpp/utils/utils.o /lab/software/apparatus3/cpp/qini/qini_utils.o 


testfourier: testfourier.o fourier.o
	g++  $^ ${objs} ${LFLAGS} -o testfourier

clean:
	rm -f *.o

.cpp.o:
	indent $<
	g++ ${CFLAGS} $< -c




