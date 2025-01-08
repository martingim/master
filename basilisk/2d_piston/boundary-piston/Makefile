KDT_LIBS=$(BASILISK)/kdt
PPR_LIBS=$(BASILISK)/ppr
GL_LIBS=$(BASILISK)/gl


.PHONY: default clean

default: piston ;

#piston: piston.c
#	CC99='gcc -std=c99' qcc piston.c -o piston -fopenmp -L$(CFD_LIBS) -lCFDwavemaker_swd_openmp -lm -L. -L$(GL_LIBS) -lglutils -lfb_osmesa -lOSMesa -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

mpi_piston: piston.c
	CC99='mpicc -std=c99' qcc -D_MPI=1 piston.c -o piston  -lm -L$(CFD_LIBS) -L. -L$(GL_LIBS) -lglutils -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

openmp_piston: piston.c
	CC99='gcc -std=c99' qcc -DTRACE=2 -fopenmp piston.c -o piston -lm -L$(CFD_LIBS) -L. -L$(GL_LIBS) -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

piston: piston.c
	CC99='gcc -std=c99' qcc piston.c -o piston -g3 -L$(CFD_LIBS) -lCFDwavemaker_swd_openmp -lm -L. -L$(GL_LIBS) -lglutils -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS) ;\
	./piston



clean:
	rm perfs ; \
	rm piston;\
	rm -r results