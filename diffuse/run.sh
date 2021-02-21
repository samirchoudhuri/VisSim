path="/home/samir/simmulti/diffmulti"
export path

# In input.grf alpha=-2, because it for Unit power spectrum, c_l is unity for all frequency and only freq. correction in sp. intensity I=nu^2.... for diff 0.8

for i in `seq 1 1`;
do
    #sed "s/-850/-2$i"79"/" input.grf >input.grfcp
    
    #$path/grf input.grfcp img_diff$i.fits
    #cp R53D22.RES.UVFITS tmp/visdiff$i.fits
    #$path/visfitsgrid tmp/img_diff$i.fits tmp/visdiff$i.fits 0.
done
