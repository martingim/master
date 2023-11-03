KDT_LIBS=$(BASILISK)/kdt
PPR_LIBS=$(BASILISK)/ppr
GL_LIBS=$(BASILISK)/gl
CFD_LIBS=$PWD

.PHONY: default clean

default: irregular ;

irregular: irregular.c
	CC99='gcc -std=c99' qcc irregular.c -o irregular -fopenmp -L$(CFD_LIBS) -lCFDwavemaker_swd_openmp -lm -L. -L$(GL_LIBS) -lglutils -lfb_osmesa -lOSMesa -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran

clean:
	rm *.o ; \
	rm irregular;
