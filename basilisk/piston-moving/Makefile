KDT_LIBS=$(BASILISK)/kdt
PPR_LIBS=$(BASILISK)/ppr
GL_LIBS=$(BASILISK)/gl


.PHONY: default clean

default: wave_wall ;

#wave_wall: wave_wall.c
#	CC99='gcc -std=c99' qcc wave_wall.c -o wave_wall -fopenmp -L$(CFD_LIBS) -lCFDwavemaker_swd_openmp -lm -L. -L$(GL_LIBS) -lglutils -lfb_osmesa -lOSMesa -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

mpi_wave_wall: wave_wall.c
	CC99='mpicc -std=c99' qcc -D_MPI=1 wave_wall.c -o wave_wall  -lm -L$(CFD_LIBS) -L. -L$(GL_LIBS) -lglutils -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

openmp_wave_wall: wave_wall.c
	CC99='gcc -std=c99' qcc -fopenmp wave_wall.c -o wave_wall -lm -L$(CFD_LIBS) -L. -L$(GL_LIBS) -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS)

wave_wall: wave_wall.c
	CC99='gcc -std=c99' qcc wave_wall.c -o wave_wall -g3 -L$(CFD_LIBS) -lCFDwavemaker_swd_openmp -lm -L. -L$(GL_LIBS) -lglutils -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran -I$(CFD_LIBS) ;\
	./wave_wall



clean:
	rm perfs ; \
	rm wave_wall;\
	rm snapshot.png;\
	rm *.csv;\
	rm vtu/ascii/*;\
	rm vtu/matlab/*;\
	rm vtu/*.vtu;
