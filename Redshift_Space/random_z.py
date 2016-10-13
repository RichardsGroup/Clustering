import os
import numpy as np
from astropy.io import fits as pf
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.interpolate as interpolate
from astroML.density_estimation import EmpiricalDistribution



#Read in the data
#data=pf.open('GTR-ADM-QSO-ir-testhighz_kdephotoz_lup_s82_quasar_candidates.fits')[1].data
data=pf.open('SpIES_SDSS_Specz20150619_Photz20160612_fixtphotz.fits')[1].data
#data=pf.open('TEST_cut_dat.fits')[1].data


#Read in the random catalog
#randat=open('highz_rand_spies.txt','r') 
randat=open('high_z_SpIES_footprint.txt','r') 
randheader=randat.readline()
randoms=randat.readlines()	
randra=np.zeros(len(randoms))
randdec=np.zeros(len(randoms))
#randz=np.zeros(len(randoms))

for i in range(len(randoms)):
	randra[i]=float(randoms[i].split()[0])
	randdec[i]=float(randoms[i].split()[1])
	#randz[i]=float(randoms[i].split()[2])



#gdx=data.zphotbest>=2.5
gdx=(data.Z>2.9) 
#factor=len(randra)/float(len(data.zphotbest[gdx]))
factor=len(randra)/float(len(data.Z[gdx]))

print factor


#Compute the optimal bin width (and number) using the Freedman-Diaconis rule 
q75, q25 = np.percentile(data.Z[gdx], [75 ,25])
iqr = q75 - q25
FDrule = 2* iqr / (float(len(data.Z[gdx]))**(1/3.0))
bnum = round((max(data.Z[gdx])-min(data.Z[gdx]))/FDrule)

print FDrule

bins=np.linspace(min(data.Z[gdx]),max(data.Z[gdx]),bnum)
binsmid = bins[:-1] + np.diff(bins)/2. #find the center of the bins

bin=np.linspace(min(data.Z[gdx]),max(data.Z[gdx]),len(bins))#*factor)
binmid = bin[:-1] + np.diff(bin)/2. #find the center of the bins

randsz = EmpiricalDistribution(data.Z[gdx]).rvs(len(data.Z[gdx])*factor)
print len(randsz)

#dat,xd = np.histogram(data.zphotbest[gdx], bins=bins)
dat,xd = np.histogram(data.Z[gdx], bins=bins)
rand,xr = np.histogram(randsz, bins=bin)


plt.figure(1)
plt.title('Redshift Distributions')
plt.plot(binsmid,dat,linestyle='steps-mid',label='data', linewidth=3)
plt.plot(binmid,rand/factor,linestyle='steps-mid',label='scaled randoms',linewidth=2)
plt.xlabel('z')
plt.ylabel('Number')
plt.legend()


RA=np.array(randra)
DEC=np.array(randdec)
Z=np.array(randsz)

print len(RA), len(DEC), len(Z)


#Write to Fits file
tbhdu=pf.BinTableHDU.from_columns([pf.Column(name='RA',format='D',array=RA),
pf.Column(name='DEC',format='D',array=DEC),
pf.Column(name='Z',format='D',array=Z)])
	
	

prihdr=pf.Header()
prihdr['COMMENT']="250,000 random points in the SpIES dual-band footprint"
prihdu=pf.PrimaryHDU(header=prihdr)

hdulist=pf.HDUList([prihdu,tbhdu])
hdulist.writeto('SpIES_random_specplphot_highz.fits')
#hdulist.writeto('SpIES_random_specplphot_highz_test.fits')

plt.show()