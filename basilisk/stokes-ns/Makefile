KDT_LIBS=$(BASILISK)/kdt
PPR_LIBS=$(BASILISK)/ppr
GL_LIBS=$(BASILISK)/gl
CFD_LIBS=/home/martin/Documents/master/basilisk/src

.PHONY: default clean

default: stokes-ns ;

stokes-ns: stokes-ns.c
	CC99='gcc -std=c99' qcc stokes-ns.c -o stokes-ns -L$(CFD_LIBS) -lCFDwavemaker_swd_openmp -lm -L. -L$(GL_LIBS) -lglutils -lfb_osmesa -lOSMesa -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

clean:
	rm *.o ; \
	rm stokes-ns;
	rm vtk/*.vtk
