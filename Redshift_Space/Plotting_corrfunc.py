from astropy.io import fits as pf
from astroML.correlation import two_point, bootstrap_two_point
import numpy as np
from matplotlib import *
from pylab import *
import pylab as plt
from time import time


##########################
###Open data set
##########################


data = open('SpIES_fullcorr_test_v1.txt','rw')
#data = open('SpIES_fullcorr_test_v3.txt','rw')

header=data.readline()

sep=[]
xifull=[]
sigma=[]
RR=[]
DD=[]
dat=data.readlines()
for i in range(len(dat)):
	sep.append(float(dat[i].split()[0]))
	xifull.append(float(dat[i].split()[4]))
	sigma.append(float(dat[i].split()[5]))
	RR.append(float(dat[i].split()[3]))
	DD.append(float(dat[i].split()[1]))
#Calculate Poisson error

Poisson_err = (1+np.array(xifull))/np.sqrt(1.0*np.array(DD))

LS_err = (1+np.array(xifull))/np.sqrt(4996.0)

print 'Poisson Errors calculated'



##############################

#OPEN THE ROSS 2009 DATA TO COMPARE
data=open('k_output_UNI22.dat','r')

r=[]
Xi=[]
errxi = []

for i in data.readlines():
	val=i.split()
	r.append(float(val[0]))
	Xi.append(float(val[2]))
	errxi.append(float(val[3]))
	
R=np.array(r)


#The Shen 2007 data
Shenxi = [-999,-999,-999,16.5,-999,3.54,1.26,0.663,0.191,0.131,0.236,-0.28,0.361,0.101,0.0384,0.0368,0.0101,0.0194,-0.00396,0.0101,-0.00296,0.00214]
Shenr = [2.244,2.825,3.557,4.477,5.637,7.096,8.934,11.25,14.16,17.83,22.44,28.25,35.57,44.77,56.37,70.96,89.34,112.5,141.6,178.3,224.4,282.5]
Shenerr = [0,0,0,12.8,0,3.61,1.88,0.733,0.786,0.472,0.175,0.223,0.170,0.121,0.0862,0.0644,0.0382,0.0250,0.0219,0.0134,0.00672,0.00953]

print 'Plotting'

params = {'legend.fontsize': 16, 'xtick.labelsize': 20, 'ytick.labelsize': 20, 'xtick.major.width':2, 'xtick.minor.width':2, 'ytick.major.width':2, 'ytick.minor.width':2, 'xtick.major.size':8, 'xtick.minor.size':6, 'ytick.major.size':8, 'ytick.minor.size':6}
plt.rcParams.update(params)
plt.rc("axes", linewidth=3.0)

figure(1,figsize=(10,10))
scatter(Shenr,Shenxi,color='r',s=80,label='Shen 2007')
errorbar(Shenr,Shenxi,yerr=Shenerr,linestyle="None",linewidth=2,color='r')
scatter(10**R,Xi, color='b',s=80,label='Ross2009')
errorbar(10**R,Xi,yerr=errxi,linestyle="None",linewidth=2,color='b')
scatter(sep,xifull,s=80,color='g',label='SpIES z>3')
errorbar(sep,xifull,yerr=sigma,linestyle="None",linewidth=2,color='g')
#errorbar(sep,xifull,yerr=Poisson_err,linestyle="None",linewidth=2,color='r')
xscale('log')
yscale('log')
ylim(10**-3,100)
xlim(1,200)
xlabel(r's (h$^{-1}$Mpc)',fontsize=16)
ylabel(r'$\xi$(s)',fontsize=16)
legend(scatterpoints=1)
savefig('SpIES_Jackknife_highz_specplphot_smallscale.png')


figure(2,figsize=(10,10))
scatter(sep,sigma/Poisson_err,color='r',s=80,label='JK/Poisson')
scatter(sep,sigma/LS_err, color='b',s=80,label='JK/LS')
xscale('log')
ylim(0,5)
legend(scatterpoints=1)

show()
