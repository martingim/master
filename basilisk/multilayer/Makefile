KDT_LIBS=$(BASILISK)/kdt
PPR_LIBS=$(BASILISK)/ppr
GL_LIBS=$(BASILISK)/gl
#CFD_LIBS=$(HOME)/CFD_wavemaker

.PHONY: default clean

default: multilayer ;

#multilayer: multilayer.c
#	CC99='gcc -std=c99' qcc multilayer.c -o multilayer -fopenmp -L$(CFD_LIBS) -lCFDwavemaker_swd_openmp -lm -L. -L$(GL_LIBS) -lglutils -lfb_osmesa -lOSMesa -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

mpi_multilayer: multilayer.c
	CC99='mpicc -std=c99' qcc -D_MPI=1 multilayer.c -o multilayer  -lm -L. -L$(GL_LIBS) -lglutils -lfb_osmesa -lOSMesa -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran

openmp_multilayer: multilayer.c
	CC99='gcc -std=c99' qcc -autolink -fopenmp multilayer.c -o multilayer -L$(CFD_LIBS) -lm -L. -L$(GL_LIBS) -lglutils -lfb_osmesa -lOSMesa -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

multilayer: multilayer.c
	CC99='gcc -std=c99' qcc multilayer.c -o multilayer -lm -L. -L$(GL_LIBS) -lglutils -lfb_osmesa -lOSMesa -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran



clean:
	rm perfs ; \
	rm multilayer;\
	rm snapshot.png;\
	rm *.csv;\