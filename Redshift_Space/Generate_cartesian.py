#LSS Clustering Code::
import numpy as np
from astropy.io import fits as pf
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM as FLCDM
#from astropy.cosmology import comoving_distance, Planck13
from astropy import coordinates as coord
from scipy import integrate, constants
import time

start=time.time()

#Import data 

obsdat=pf.open('SpIES_SDSS_Specz20150619_Photz20160612_fixtphotz.fits')[1].data
#obsdat=pf.open('TEST_cut_dat.fits')[1].data
randat=pf.open('SpIES_random_specplphot_highz.fits')[1].data
#randat=pf.open('SpIES_random_specplphot_highz_test.fits')[1].data

end1=time.time()
print end1-start, 'Files read'

#Create arrays of zeros which we use to fill with data
obsra=obsdat.RA
obsdec=obsdat.DEC
obsz=obsdat.Z
randra=randat.RA
randdec=randat.DEC
randz=randat.Z


#Convert redshift to Mpc
#Create the 2SLAQ cosmology (H0=100h^-1, Om=0.26, Ol=0.74, Tcmb=2.725)
slaqcosmo=FLCDM(100,0.26,2.725)
obsX=slaqcosmo.comoving_distance(obsz)
randX=slaqcosmo.comoving_distance(randz)
end2=time.time()
print end2-start, 'Comoving distances calculated'



#Convert RA/DEC to cartesian coords
cutcoord=coord.SkyCoord(ra=obsra*u.degree,dec=obsdec*u.degree, distance=obsX,frame='icrs')
randcoord=coord.SkyCoord(ra=randra*u.degree,dec=randdec*u.degree, distance=randX,frame='icrs')
cx=cutcoord.cartesian.x
cy=cutcoord.cartesian.y
cz=cutcoord.cartesian.z
rx=randcoord.cartesian.x
ry=randcoord.cartesian.y
rz=randcoord.cartesian.z

'''
newfile=open('cut_cartesian','w')
newfile.write('cutx, cuty, cutz')
for i in range(len(cx)):
	newfile.write(str(cx[i])+','+str(cy[i])+','+str(cz[i]))

rnewfile=open('rand_cartesian','w')
rnewfile.write('randx, randy, randz')
for i in range(len(rx)):
	newfile.write(str(rx[i])+','+str(ry[i])+','+str(rz[i]))
'''

tbhdu=pf.BinTableHDU.from_columns([pf.Column(name='datx',format='D',array=cx),
pf.Column(name='daty',format='D',array=cy),pf.Column(name='datz',format='D',array=cz),pf.Column(name='randx',format='D',array=rx),
pf.Column(name='randy',format='D',array=ry), pf.Column(name='randz',format='D',array=rz)])

prihdr=pf.Header()
prihdr['COMMENT']="Catalog of high redshift quasars in the SpIES field"
prihdu=pf.PrimaryHDU(header=prihdr)

hdulist=pf.HDUList([prihdu,tbhdu])
hdulist.writeto('SpIES_highz_cartesian_coords_specplphot.fits')
#hdulist.writeto('SpIES_highz_cartesian_coords_specplphot_test.fits')


