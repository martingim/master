KDT_LIBS=$(BASILISK)/kdt
PPR_LIBS=$(BASILISK)/ppr
GL_LIBS=$(BASILISK)/gl


.PHONY: default clean

default: stokes ;

#stokes: stokes.c
#	CC99='gcc -std=c99' qcc stokes.c -o stokes -fopenmp -L$(CFD_LIBS) -lCFDwavemaker_swd_openmp -lm -L. -L$(GL_LIBS) -lglutils -lfb_osmesa -lOSMesa -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

mpi_stokes: stokes.c
	CC99='mpicc -std=c99' qcc -D_MPI=1 stokes.c -o stokes  -lm -L$(CFD_LIBS) -L. -L$(GL_LIBS) -lglutils -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

openmp_stokes: stokes.c
	CC99='gcc -std=c99' qcc -DTRACE=2 -fopenmp stokes.c -o stokes -lm -L$(CFD_LIBS) -L. -L$(GL_LIBS) -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

test_stokes: stokes-ns.c
	CC99='gcc -std=c99' qcc -DTRACE=2 -fopenmp stokes-ns.c -o test_stokes -lm -L$(CFD_LIBS) -L. -L$(GL_LIBS) -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)


stokes: stokes.c
	CC99='gcc -std=c99' qcc stokes.c -o stokes -g3 -L$(CFD_LIBS) -lCFDwavemaker_swd_openmp -lm -L. -L$(GL_LIBS) -lglutils -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS) ;\
	./stokes



clean:
	rm perfs ; \
	rm stokes;\
	rm -r results;\
	rm *.vtu
