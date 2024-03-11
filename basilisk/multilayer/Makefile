KDT_LIBS=$(BASILISK)/kdt
PPR_LIBS=$(BASILISK)/ppr
GL_LIBS=$(BASILISK)/gl
CFD_LIBS=/home/martin/Documents/master/basilisk/src

.PHONY: default clean

default: multilayer ;

multilayer: multilayer.c
	CC99='gcc -std=c99' qcc multilayer.c -o multilayer -L$(CFD_LIBS) -lCFDwavemaker_swd_openmp -lm -L. -L$(GL_LIBS) -lglutils -lfb_osmesa -lOSMesa -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

bar: bar.c
	CC99='gcc -std=c99' qcc bar.c -o bar -L$(CFD_LIBS) -lCFDwavemaker_swd_openmp -lm -L. -L$(GL_LIBS) -lglutils -lfb_osmesa -lOSMesa -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)


clean:
	rm *.o ; \
	rm multilayer;\