all: analyze basler probe 

CCFITS_INC = /usr/local/include/CCfits
GSL_INC = gsl/include

INC = -Iqini -Ifits -Ifuncs -Iutils -I${CCFITS_INC} -I${GSL_INC} 

GSL_LIB = gsl/lib

CFLAGS =  -Wall ${INC} -fopenmp
LFLAGS = -L${GSL_LIB} -lgsl -lgslcblas -lm -lCCfits -ltiff -fopenmp

objs = funcs/funcs.o utils/utils.o qini/qini_utils.o fits/fits.o 
 
analyze: analyze.o ${objs} vt100_macros.h Fermions.h
	g++ analyze.o ${objs} ${LFLAGS} -o analyze
	chmod a+w analyze 
	cp -v analyze /lab/software/apparatus3/bin/analyze
	

basler: basler.o ${objs} vt100_macros.h Fermions.h
	g++ basler.o ${objs} ${LFLAGS} -o basler
	chmod a+w basler
	cp -v basler /lab/software/apparatus3/bin/basler

probe: probe.o ${objs} vt100_macros.h Fermions.h
	g++ probe.o ${objs} ${LFLAGS} -o probe
	chmod a+w probe
	cp -v probe /lab/software/apparatus3/bin/probe

clean:
	rm -f *.o

.cpp.o: 
	indent $<
	g++ ${CFLAGS} $< -c




