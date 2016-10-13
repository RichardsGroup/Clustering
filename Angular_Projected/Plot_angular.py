import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
#from matplotlib.colors import LogNorm
from palettable.colorbrewer.sequential import Oranges_9 as col
from scipy.special import gamma
from scipy.optimize import curve_fit

Corrdat = open('SpSh_angular_Clustering_nooutlier_correctedz_v2.txt','rw')

def Power(r,rp,gam):
	realxi = ((np.sqrt(np.pi)*gamma((gam-1.0)/2.0)/gamma(gam/2.0)))*(r/rp)**-gam
	return realxi
	
def Powerfixed(r,rp):
	realxi =((np.sqrt(np.pi)*gamma(1/2.0)/gamma(1.0)))*(r/rp)**-2
	return realxi

header = Corrdat.readline()


th = []
xi = []
DD = []
DR = []
RR = []
for i in Corrdat.readlines():
	t = float(i.split()[0])
	x = float(i.split()[4])
	dd = float(i.split()[1])
	dr = float(i.split()[2])
	rr = float(i.split()[3])
	th.append(t)
	xi.append(x)
	DD.append(dd)
	DR.append(dr)
	RR.append(rr)
	
size = np.sqrt(float(len(th)))
RRpairs=np.asarray(RR)

totRR=np.asarray(np.sum(RRpairs,axis=0))

####################################
#POISSON ERROR BAR
####################################

totQQ = np.sum(DD)
om = np.asarray(xi)
Errors = (1 + om)/np.sqrt(totQQ/11)
poisson = np.asarray(Errors)


####################################
#ERROR BARS FROM JACKKNIFE
###################################

jkdat = open('SpSh_angular_Clustering_Jackknifedat_nooutlier_correctedz_v2.txt','rw')


#Read the header of the jackknife file and assign the number of jackknife resamples done for the sample
jkhead = jkdat.readline()
jknum= np.float(jkhead.split()[7])

#Pull out the proper info from the file (RR,Xi,etc.) for each jackknife and put into array
RRjk = np.zeros([jknum,size**2])
Xijk = np.zeros([jknum,size**2])
Thjk = np.zeros([jknum,size**2])

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
		Thjk[row,k] = float(num[j].split()[0])
		RRjk[row,k] = float(num[j].split()[3])
		Xijk[row,k] = float(num[j].split()[4])
		k+=1


#Compute the main diagonal of the covariance matrix
c=0
C=[]
for i in range(len(Thjk)):

	sig = (np.sqrt(np.sum(RRjk[i])/totRR)*(Xijk[i]-om))**2
	c+=sig
	if i+1 == len(Xijk):
		C.append(c)
print C


'''

#################################
#FIT USING CURVE_FIT
##################################

x = Th[(RP>3.5) & (RP<30)]
y = WpRp[(RP>3.5) & (RP<30)]


# Fit1
popt, pcov = curve_fit(Power, x, y,p0=[5,2])

perr = np.sqrt(np.diag(pcov))


print 'rp=',popt[0]
print 'gamma=',popt[1]
print pcov

# Fit2
popt, pcov = curve_fit(Powerfixed, x, y,p0=[5])

perr = np.sqrt(np.diag(pcov))
print 'rp=',popt[0]
print pcov
'''


#Myers2006 Line
t = np.linspace(-1.40,2.4,15)
thet=10**t
omegaall = (0.056) * thet**(-0.91)
omegasmall=(0.11) * (10**np.linspace(-1.4,0,10))**(-0.4)
omegamid = (0.035) * (10**np.linspace(0,1.69,10))**(-0.55)
omegalarge=(0.066) * (10**np.linspace(0,2.3,10))**(-0.98)
M2007 = (0.0493) * thet**(-0.928)

params = {'legend.fontsize': 16, 'xtick.labelsize': 20, 'ytick.labelsize': 20, 'xtick.major.width':2, 'xtick.minor.width':2, 'ytick.major.width':2, 'ytick.minor.width':2, 'xtick.major.size':8, 'xtick.minor.size':6, 'ytick.major.size':8, 'ytick.minor.size':6}
plt.rcParams.update(params)
plt.rc("axes", linewidth=3.0)


plt.figure(1,figsize=[10,10])

x1 = np.linspace(4,25,100)
x = np.linspace(4,150,100)
Efit = ((np.sqrt(np.pi)*gamma(1/2.0)/gamma(1.0)))*(x1/8.12)**-2
Sfit = ((gamma(0.5)*gamma(1.33/2.0))/gamma(2.33/2.0)) *(x/16.1)**-2.33

plt.scatter(th,xi,s = 100, edgecolor='None',label=('SpIES & SHELA 2.9<z<5'))
plt.errorbar(th,xi,yerr=C[0]**0.5,elinewidth=3,fmt=',')
plt.plot(thet,omegaall,linewidth=2,label = 'Myers2006 all scale')
plt.plot(10**np.linspace(-1.4,0,10),omegasmall,linewidth=2,linestyle='-.',dashes = [8,4,2,4],label = 'Myers 2006 small scale')
plt.plot(10**np.linspace(0,1.69,10),omegamid,linewidth=2,linestyle='--',dashes = [8,4,8,4],label = 'Myers 2006 med scale')
plt.plot(10**np.linspace(0,2.3,10),omegalarge,linewidth=2,linestyle=':',dashes = [2,4,2,4],label = 'Myers 2006 large scale')
plt.plot(thet,M2007,linewidth=2, color = 'm',label = 'Myers2007 all scale')
'''
plt.scatter(Er,EwEr,c='r',marker='o',label=('Eftekharzadeh 2015'))
plt.errorbar(Er,EwEr,yerr=np.asarray(eerr)/np.asarray(Er),c='r',fmt='o')
plt.plot(x1,Efit,c='r',linestyle='--',label=('Eftekharzadeh 2015 fit'))
plt.plot(x,Sfit,c='g',linestyle='--',label=('Shen 2007 fit'))
plt.scatter(sr,swsr,c='g',marker='o',label=('Shen 2007'))
plt.errorbar(sr,swsr,yerr=serr,c='g',fmt='o')
#plot the fit to SpIES data
#plt.plot(RP[(RP>3.5) & (RP<30)],Power(RP[(RP>3.5) & (RP<30)],popt[0],popt[1]))
plt.plot(RP[(RP>3.5) & (RP<30)],Powerfixed(RP[(RP>3.5) & (RP<30)],popt[0]))
'''
#plt.axvline(4)
#plt.axvline(25)
plt.xlim(10**-1,500)
plt.ylim(4*10**-5,10)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$\theta$ (Arcmin)',fontsize = 16)
plt.ylabel(r'$\omega (\theta)$',fontsize = 16)
plt.legend(scatterpoints=1)
plt.savefig('SpSh_angular_corr_JKerr_nooutlier_v2.png')

plt.show()
