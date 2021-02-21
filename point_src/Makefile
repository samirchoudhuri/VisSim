LINKLIB=  -L/home/samir/astro/lib  -fopenmp -lfftw3 -lcfitsio -lnrcp -lm -lgsl -lgslcblas
INCLUDE=-I/home/samir/astro/include/


turget:	clean visfits

clean:
	rm -rf visfits *~

visfits: beam.c visfits.c funcs.c utils.c
	gcc -g -o visfits $(INCLUDE) visfits.c beam.c funcs.c utils.c ${LINKLIB}

