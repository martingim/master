KDT_LIBS=$(BASILISK)/kdt
PPR_LIBS=$(BASILISK)/ppr
GL_LIBS=$(BASILISK)/gl


.PHONY: default clean

default: moving_piston ;

#moving_piston: moving_piston.c
#	CC99='gcc -std=c99' qcc moving_piston.c -o moving_piston -fopenmp -L$(CFD_LIBS) -lCFDwavemaker_swd_openmp -lm -L. -L$(GL_LIBS) -lglutils -lfb_osmesa -lOSMesa -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

mpi_moving_piston: moving_piston.c
	CC99='mpicc -std=c99' qcc -D_MPI=1 moving_piston.c -o moving_piston  -lm -L$(CFD_LIBS) -L. -L$(GL_LIBS) -lglutils -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

openmp_moving_piston: moving_piston.c
	CC99='gcc -std=c99' qcc -DTRACE=2 -fopenmp moving_piston.c -o moving_piston -lm -L$(CFD_LIBS) -L. -L$(GL_LIBS) -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

moving_piston: moving_piston.c
	CC99='gcc -std=c99' qcc moving_piston.c -o moving_piston -g3 -L$(CFD_LIBS) -lCFDwavemaker_swd_openmp -lm -L. -L$(GL_LIBS) -lglutils -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS) ;\
	./moving_piston



clean:
	rm perfs ; \
	rm moving_piston;\
	rm -r results