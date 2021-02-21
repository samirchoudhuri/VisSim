LINKLIB=  -L/home/samir/astro/lib  -fopenmp -lfftw3 -lcfitsio -lm -lnrcp -lgsl -lgslcblas
INCLUDE=-I/home/samir/astro/include/


target: grf grf_unity grf_rea grf_unity_rea visfitsgrid gen_bubble

grf:  grf.c fitsprog.c beam.c
	gcc -g -o grf $(INCLUDE) fitsprog.c beam.c grf.c $(LINKLIB)
	rm -rf *~
grf_unity: grf_unity.c fitsprog.c beam.c
	gcc -g -o grf_unity $(INCLUDE) fitsprog.c beam.c grf_unity.c $(LINKLIB)
	rm -rf *~
grf_rea: grf_rea.c fitsprog.c beam.c
	gcc -g -o grf_rea $(INCLUDE) fitsprog.c beam.c grf_rea.c $(LINKLIB)
	rm -rf *~
grf_unity_rea: grf_unity_rea.c fitsprog.c beam.c
	gcc -g -o grf_unity_rea $(INCLUDE) fitsprog.c beam.c grf_unity_rea.c $(LINKLIB)
	rm -rf *~
visfitsgrid: fitsprog.c visfitsgrid.c read_fits_func.c 
	gcc -g -o visfitsgrid $(INCLUDE) fitsprog.c visfitsgrid.c read_fits_func.c $(LINKLIB)
	rm -rf *~
gen_bubble: gen_bubble.c funcs_comov.c
	gcc -g -o gen_bubble $(INCLUDE) funcs_comov.c simp.c gen_bubble.c $(LINKLIB)


clean:
	rm -rf grf grf_unity grf_rea grf_unity_rea visfitsgrid gen_bubble
