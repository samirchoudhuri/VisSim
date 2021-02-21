NTHREAD=25
PS=pttest.dat
make

time ./visfits INP.VISFITS scan.vf source.vf difftest.fits ${PS} visdiff.fits ${NTHREAD}