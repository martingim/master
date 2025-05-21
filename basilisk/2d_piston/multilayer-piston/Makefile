KDT_LIBS=$(BASILISK)/kdt
PPR_LIBS=$(BASILISK)/ppr
GL_LIBS=$(BASILISK)/gl
#CFD_LIBS=$(HOME)/Documents/master/basilisk/src

.PHONY: default clean

default: multilayer ;

#multilayer: multilayer.c
#	CC99='gcc -std=c99' qcc multilayer.c -o multilayer -fopenmp -L$(CFD_LIBS) -lCFDwavemaker_swd_openmp -lm -L. -L$(GL_LIBS) -lglutils -lfb_osmesa -lOSMesa -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

mpi_multilayer: multilayer.c
	CC99='mpicc -std=c99' qcc -D_MPI=3 -autolink multilayer.c -o multilayer -lm -L. -L$(GL_LIBS) -lglutils -lfb_osmesa -lOSMesa -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran

openmp_multilayer: multilayer.c
	CC99='gcc -std=c99' qcc -autolink -fopenmp -g multilayer.c -o multilayer  -fopenmp -lm -L. -L$(CFD_LIBS) -L$(GL_LIBS) -lglutils  -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS) -I$(GL_LIBS)

multilayer: multilayer.c
	CC99='gcc -std=c99' qcc multilayer.c -o multilayer -lm -L. -L$(CFD_LIBS) -L$(GL_LIBS) -lglutils -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

stokes: stokes.c
	CC99='gcc -std=c99' qcc stokes.c -o stokes -lm -L. -L$(CFD_LIBS) -L$(GL_LIBS) -lglutils -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)


clean:
	rm perfs ; \
	rm multilayer;\
	rm snapshot.png;\
	rm *.csv;\
	rm *.vts;
