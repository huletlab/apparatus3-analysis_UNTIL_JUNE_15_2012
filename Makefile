all: analyze basler 

INI_INC = /lab/software/apparatus3/cpp/qini
CCFITS_INC = /usr/local/include/CCfits
GSL_INC = /lab/software/apparatus3/cpp/gsl/include

INC = -I${INI_INC} -I${CCFITS_INC} -I${GSL_INC} 

GSL_LIB = /lab/software/apparatus3/cpp/gsl/lib

CFLAGS = -Wall ${INC} 
LFLAGS = -L${GSL_LIB} -lgsl -lgslcblas -lm -lCCfits -ltiff

objs = /lab/software/apparatus3/cpp/utils/utils.o /lab/software/apparatus3/cpp/qini/qini_utils.o /lab/software/apparatus3/cpp/fits/fits.o 
 
analyze: analyze.o ${objs} vt100_macros.h Fermions.h
	g++ analyze.o ${objs} ${LFLAGS} -o analyze
	chmod a+w analyze 
	cp -v analyze /lab/software/apparatus3/bin/analyze

basler: basler.o ${objs} vt100_macros.h Fermions.h
	g++ basler.o ${objs} ${LFLAGS} -o basler
	chmod a+w basler
	cp -v basler /lab/software/apparatus3/bin/basler

clean:
	rm -f *.o

.cpp.o: 
	indent $<
	g++ ${CFLAGS} $< -c



