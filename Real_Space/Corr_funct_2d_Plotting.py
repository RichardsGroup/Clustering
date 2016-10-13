import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
#from matplotlib.colors import LogNorm
from palettable.colorbrewer.sequential import Oranges_9 as col
from scipy.special import gamma
from scipy.optimize import curve_fit

Corrdat = open('SpIES_RealSpace_Clustering_v2.txt','rw')

def Power(r,rp,gam):
	realxi = ((np.sqrt(np.pi)*gamma((gam-1.0)/2.0)/gamma(gam/2.0)))*(r/rp)**-gam
	return realxi
	
def Powerfixed(r,rp):
	realxi =((np.sqrt(np.pi)*gamma(1/2.0)/gamma(1.0)))*(r/rp)**-2
	return realxi

header = Corrdat.readline()

P = []
rp = []
xi = []
DD = []
DR = []
RR = []
for i in Corrdat.readlines():
	p = float(i.split()[0])
	r = float(i.split()[1])
	x = float(i.split()[5])
	dd = float(i.split()[2])
	dr = float(i.split()[3])
	rr = float(i.split()[4])
	P.append(p)
	rp.append(r)
	xi.append(x)
	DD.append(dd)
	DR.append(dr)
	RR.append(rr)
	
size = np.sqrt(float(len(rp)))
#X = np.zeros((size,size))
#Y = np.zeros((size,size))
#Z = np.zeros((size,size))



X = np.split(np.asarray(rp),size)
Y = np.split(np.asarray(P),size)
Z = np.split(np.asarray(xi),size)
QQ = np.split(np.asarray(DD),size)
RRpairs=np.split(np.asarray(RR),size)
fact = np.ones(np.shape(X))
totRR=np.asarray(np.sum(RRpairs,axis=0))

z = np.log10(Z)


wp = 2 * np.asarray(np.sum(Z,axis=0))
rp = np.asarray(Y)[:,0]
angcor = wp/rp
WpRp = angcor#[angcor>0]
RP = rp#[angcor>0]

####################################
#POISSON ERROR BAR
####################################

totQQ = np.matrix(QQ).sum(axis=1)
Errors = (1 + wp)/np.sqrt(totQQ/11)
poisson = np.asarray(Errors.T[angcor.T>0])[0]


####################################
#ERROR BARS FROM JACKKNIFE
###################################

jkdat = open('SpIES_RealSpace_Clustering_Jackknifedat_v2.txt','rw')


#Read the header of the jackknife file and assign the number of jackknife resamples done for the sample
jkhead = jkdat.readline()
jknum= np.float(jkhead.split()[8])-1

#Pull out the proper info from the file (RR,Xi,etc.) for each jackknife and put into array
RRjk = np.zeros([jknum,size**2])
Xijk = np.zeros([jknum,size**2])
Pijk = np.zeros([jknum,size**2])
Rpjk = np.zeros([jknum,size**2])

num=jkdat.readlines()
row = 0
k=0


for j in range(len(num)):
	#At the end of a Jackknife run, I create a new line called 'STEP_DONE' which tells me that my next JK routine is running and separates the two.
	if num[j].split()[0] == 'STEP_DONE':
		row +=1
		k=0
	#For each of the JK runs, append the info into arrays
	else:
		Rpjk[row,k] = float(num[j].split()[1])
		Pijk[row,k] = float(num[j].split()[0])
		RRjk[row,k] = float(num[j].split()[4])
		Xijk[row,k] = float(num[j].split()[5])
		k+=1

#Using these arrays, compute Wp(rp) for the JK samples

WpRpjk=[]
RPjk = []
jkRR=[]

for i in range(len(Rpjk)):
	rpv = Rpjk[i]
	piv = Pijk[i]
	rrv = RRjk[i]
	xiv = Xijk[i]
	size = np.sqrt(float(len(rpv)))
	across=np.split(np.asarray(rpv),size)
	along=np.split(np.asarray(piv),size)
	rand=np.split(np.asarray(rrv),size)
	tpcf=np.split(np.asarray(xiv),size)
	jkRR.append(np.asarray(np.sum(rand,axis=0)))
	wpjk = 2 * np.asarray(np.sum(tpcf,axis=0))
	rpjk = np.asarray(across[0])
	angcorjk = wpjk/rpjk
	WpRpjk.append(np.asarray(angcorjk))#[angcorjk>0]))
	RPjk.append(np.asarray(rpjk))#[angcorjk>0]))


#Compute the main diagonal of the covariance matrix
c=0
C=[]
for i in range(len(WpRpjk)):
	#print WpRpjk[i]
	#print WpRp
	sig = (np.sqrt(jkRR[i]/totRR)*(WpRpjk[i]-WpRp))**2
	c+=sig
	if i+1 == len(WpRpjk):
		C.append(c)
print C
#print WpRp




#################################
#FIT USING CURVE_FIT
##################################

x = RP[(RP>3.5) & (RP<30)]
y = WpRp[(RP>3.5) & (RP<30)]


# Fit1
popt, pcov = curve_fit(Power, x, y,p0=[5,2])

perr = np.sqrt(np.diag(pcov))


print 'rp=',popt[0]
print 'gamma=',popt[1]
print pcov
'''
# Fit2
popt, pcov = curve_fit(Powerfixed, x, y,p0=[5])

perr = np.sqrt(np.diag(pcov))
print 'rp=',popt[0]
print pcov
'''




###################### Sarah's 2015 paper data #######################
Er = [2.09,2.53,3.06,3.70,4.47,5.41,6.55,7.92,9.58,11.6,14.03,16.97,20.54,24.85,30.06,36.37,44.0,53.23,64.4]
Ew = [69.353,60.807,57.125,52.908,45.454,39.291,34.392,25.846,23.907,16.056,11.771,13.753,9.404,6.824,6.554,7.341,4.192,4.483,2.576]

EwEr= np.asarray(Ew)/np.asarray(Er)
eerr = [19.294,14.491,10.707,8.654,5.981,10.607,9.142,6.051,5.210,3.002,1.854,5.348,3.130,0.958,2.099,1.422,1.061,1.076,0.159]
##################### Shen 2007 paper data #################
sr = [1.679,2.371,3.350,4.732,6.683,9.441,13.34,18.84,26.61,37.58,53.09,74.99,105.9,149.6,211.3]
swsr=[154,236,78.1,91.3,15.7,10.6,3.06,-0.681,0.516,0.437,0.0675,0.0484,0.0674,0.0228,-0.0183]
serr = [162,195,51.5,41.6,7.81,4.45,2.85,0.913,0.810,0.395,0.259,0.145,0.0592,0.0292,0.00992]

plt.figure(1,figsize=[10,8])
#l = np.linspace(-1,2,8)
l=[-1,-0.5,0,0.5,1,1.5,2]
plt.contourf(X,Y,z,levels = l,cmap = col.mpl_colormap)
plt.contourf(fact*-1*X,Y,z,levels = l,cmap = col.mpl_colormap)
plt.contourf(X,fact*-1*Y,z,levels = l,cmap = col.mpl_colormap)
plt.contourf(fact*-1*X,fact*-1*Y,z,levels = l,cmap = col.mpl_colormap)
plt.xlim(-20,20)
plt.ylim(-20,20)
plt.xlabel(r'r$_{p}$ (h$^{-1}$ Mpc)',fontsize = 14)
plt.ylabel(r'$\pi$ (h$^{-1}$ Mpc)',fontsize = 14)
cb=plt.colorbar(cmap = col.mpl_colormap)
cb.set_label(label=r'log$_{10}(\xi$) (r$_p$,$\pi$)',size=14)
#plt.savefig('Test_2Dcorr_SpIES.pdf')


plt.figure(2,figsize=[10,8])

x1 = np.linspace(4,25,100)
x = np.linspace(4,150,100)
Efit = ((np.sqrt(np.pi)*gamma(1/2.0)/gamma(1.0)))*(x1/8.12)**-2
Sfit = ((gamma(0.5)*gamma(1.33/2.0))/gamma(2.33/2.0)) *(x/16.1)**-2.33

plt.scatter(RP,WpRp,label=('SpIES match to SDSS'))
plt.errorbar(RP,WpRp,yerr=C[0]**0.5/RP,fmt='o')
plt.scatter(Er,EwEr,c='r',marker='o',label=('Eftekharzadeh 2015'))
plt.errorbar(Er,EwEr,yerr=np.asarray(eerr)/np.asarray(Er),c='r',fmt='o')
plt.plot(x1,Efit,c='r',linestyle='--',label=('Eftekharzadeh 2015 fit'))
plt.plot(x,Sfit,c='g',linestyle='--',label=('Shen 2007 fit'))
plt.scatter(sr,swsr,c='g',marker='o',label=('Shen 2007'))
plt.errorbar(sr,swsr,yerr=serr,c='g',fmt='o')
#plot the fit to SpIES data
#plt.plot(RP[(RP>3.5) & (RP<30)],Power(RP[(RP>3.5) & (RP<30)],popt[0],popt[1]))
plt.plot(RP[(RP>3.5) & (RP<30)],Powerfixed(RP[(RP>3.5) & (RP<30)],popt[0]))

#plt.axvline(4)
#plt.axvline(25)
plt.xlim(2,80)
plt.ylim(10**-3,10**3)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'r$_{p}$ (h$^{-1}$ Mpc)',fontsize = 14)
plt.ylabel(r'w$_{p}$(r$_{p}$)/r$_{p}$ ',fontsize = 14)
plt.legend(scatterpoints=1)
#plt.savefig('Test_1Dcorr_JKerr.pdf')

plt.show()
