KDT_LIBS=$(BASILISK)/kdt
PPR_LIBS=$(BASILISK)/ppr
GL_LIBS=$(BASILISK)/gl
CFD_LIBS=/home/martin/Documents/master/basilisk/src

.PHONY: default clean

default: bar ;

bar: bar.c
	CC99='gcc -std=c99' qcc bar.c -o bar -fopenmp -L$(CFD_LIBS) -lCFDwavemaker_swd_openmp -lm -L. -L$(GL_LIBS) -lglutils -lfb_osmesa -lOSMesa -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

clean:
	rm *.o ; \
	rm bar;
