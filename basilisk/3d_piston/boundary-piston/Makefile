KDT_LIBS=$(BASILISK)/kdt
PPR_LIBS=$(BASILISK)/ppr
GL_LIBS=$(BASILISK)/gl


.PHONY: default clean

default: piston3d ;

#piston3d: piston3d.c
#	CC99='gcc -std=c99' qcc piston3d.c -o piston3d -fopenmp -L$(CFD_LIBS) -lCFDwavemaker_swd_openmp -lm -L. -L$(GL_LIBS) -lglutils -lfb_osmesa -lOSMesa -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

mpi_piston3d: piston3d.c
	CC99='mpicc -std=c99' qcc -D_MPI=1 piston3d.c -o piston3d  -lm -L$(CFD_LIBS) -L. -L$(GL_LIBS) -lglutils -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

openmp_piston3d: piston3d.c
	CC99='gcc -std=c99' qcc -DTRACE=2 -fopenmp piston3d.c -o piston3d -lm -L$(CFD_LIBS) -L. -L$(GL_LIBS) -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

piston3d: piston3d.c
	CC99='gcc -std=c99' qcc piston3d.c -o piston3d -g3 -L$(CFD_LIBS) -lCFDwavemaker_swd_openmp -lm -L. -L$(GL_LIBS) -lglutils -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS) ;\
	./piston3d



clean:
	rm perfs ; \
	rm piston3d;\
	rm -r results