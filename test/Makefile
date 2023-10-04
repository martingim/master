KDT_LIBS=$(BASILISK)/kdt
PPR_LIBS=$(BASILISK)/ppr
GL_LIBS=$(BASILISK)/gl

.PHONY: default clean

default: bar ;

irregular: bar.c
	CC99='gcc -std=c99' qcc bar.c -o bar -fopenmp -lm -L. -L$(GL_LIBS) -lglutils -lfb_osmesa -lOSMesa -lGLU -ldl -lpthread -Wall -O2 -lstdc++ -L$(PPR_LIBS) -lppr -lgfortran
	rm *.o ; \
	rm bar;

clean:
	rm *.o; \
	rm bar;