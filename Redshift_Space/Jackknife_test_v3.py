from astropy.io import fits as pf
from astroML.correlation import two_point, bootstrap_two_point
import numpy as np
from matplotlib import *
from pylab import *
from time import time
from operator import itemgetter

start=time()

#Import data here
#obs=pf.open('SpIES_highz_cartesian_coords.fits')[1].data
obs=pf.open('SpIES_highz_cartesian_coords_specplphot.fits')[1].data

bdx=(obs.datx!=0)

############################################################################################
#COMPUTE THE CORRELATION FUNCTION
#Combine the x,y,z data from the imported file in the form [[x1,y1,z1],[x2,y2,z2]...]
obsX=np.concatenate((np.array([obs.datx[bdx]]).T,np.array([obs.daty[bdx]]).T,np.array([obs.datz[bdx]]).T),axis=1)
randX=np.concatenate((np.array([obs.randx]).T,np.array([obs.randy]).T,np.array([obs.randz]).T),axis=1)

print 'Data Read in'

#2slaq radii
#rad=np.linspace(-0.8,2.3,32)
#radii=10**rad

#BOSS radii
#radii=np.arange(2,203,8)


#Ross 2009 radii
#rad=np.arange(0.1,4.1,0.1)
#radii=10**rad

#Shen 2007 approx radii
rad = np.arange(0.3,2.6,0.1)
radii=10**rad


TPCF=two_point(obsX,radii,method='landy-szalay',data_R=randX,random_state=None)
xifull=np.asarray(TPCF[0])
DDfull=np.asarray(TPCF[1])
DRfull=np.asarray(TPCF[2])
RRfull=np.asarray(TPCF[3])
print 'Full 2PCF calculated'
end = time()

print end-start,'s'

#############################################################
#calculate Poisson Errors using Equation 13 in Ross 2009
#CALCULATE POISSON ERRORS USING EQUATION 13 IN ROSS 2009
Poisson_err = (1+xifull)/np.sqrt(1.0*DDfull)

print 'Poisson Errors calculated'

########################################################
#CALCULATE THE CENTERS OF THE BINS
sep=[]
	
for i in range(len(radii)-1):
	a=radii[i]
	b=radii[i+1]
	sep.append((a+b)/2.)
print 'Separations Calculated'


############################################################
###JACKKNIFE CODE
##FIRST SEPARATE THE X,Y,X DATA AND RANDOMS
xdat=np.array([obs.datx[bdx]])
ydat=np.array([obs.daty[bdx]])
zdat=np.array([obs.datz[bdx]])
xrand=np.array([obs.randx])
yrand=np.array([obs.randy])
zrand=np.array([obs.randz])

xd = xdat[0]
yd = ydat[0]
zd = zdat[0]
xr = xrand[0]
yr = yrand[0]
zr = zrand[0]

start2=time()

#oX = sorted(obsX, key = itemgetter(2,1,0))
#rX = sorted(randX, key = itemgetter(2,1,0))

#split= np.linspace(0,len(randX),10)

zmin= min(np.hstack((zr,zd)))

zmax = max(np.hstack((zr,zd)))

xmax = max(np.hstack((xr,xd)))

slopes = np.linspace(zmin,zmax,10)/float(xmax)

XIjk=[]
RRjk=[]
ctd = 0
ctr = 0
for i in range(len(slopes)):

	if i == 0:

		continue
	
	elif i == 1:
		print "first",i
		cutd = (obsX[:,2] <= slopes[i]*obsX[:,0])
		cutr = (randX[:,2] <= slopes[i]*randX[:,0])
		 
		ctd += len(obsX[cutd])/float(len(obsX))
		ctr += len(randX[cutr])/float(len(randX))
		
	elif i == len(slopes)-1:
		print "last",i
		
		cutd = (obsX[:,2] >= slopes[i-1]*obsX[:,0])
		cutr = (randX[:,2] >= slopes[i-1]*randX[:,0])
		
		ctd += len(obsX[cutd])/float(len(obsX))
		ctr += len(randX[cutr])/float(len(randX))
		
	else:
		print i
		
		cutd = (obsX[:,2] >= slopes[i-1]*obsX[:,0]) & (obsX[:,2] <= slopes[i]*obsX[:,0])
		cutr = (randX[:,2] >= slopes[i-1]*randX[:,0]) & (randX[:,2] <= slopes[i]*randX[:,0])	
		
		ctd += len(obsX[cutd])/float(len(obsX))
		ctr += len(randX[cutr])/float(len(randX))
		
	print ctd
	print ctr
	
	dcut=np.invert(cutd)
	rcut=np.invert(cutr)
	
	oX = obsX[dcut]
	rX = randX[rcut]
	
	print "data ratio=", 1.0-len(obsX[dcut])/float(len(obsX))
	print "random ratio=", 1.0-len(randX[rcut])/float(len(randX))
	
	#scatter(obsX[:,0],obsX[:,2],s=1, c='g',edgecolor='g')
	#scatter(rX[:,0],rX[:,2],s=1,c='r', edgecolor='r')
	#scatter(oX[:,0],oX[:,2],s=1,c='b')
	#show()
	
	TPCF=two_point(oX,radii,method='landy-szalay',data_R=rX, random_state=None)
	xi=np.asarray(TPCF[0])
	dd=np.asarray(TPCF[1])
	dr=np.asarray(TPCF[2])
	rr=np.asarray(TPCF[3])
	RRjk.append(rr)
	XIjk.append(xi)
	


'''
#SPLIT INTO THE 'SPLITVAL' AMOUNT OF SUBSETS
xsplitval=5
ysplitval=3
xspvals= np.linspace(min(np.hstack((zr,zd))),max(np.hstack((zr,zd))),xsplitval)
yspvals= np.linspace(min(np.hstack((yr,yd))),max(np.hstack((yr,yd))),ysplitval)

XIjk=[]
RRjk=[]



for i in range(xsplitval-1):
	for j in range(ysplitval-1):
		print i, j 
		cut1 = (zd>=xspvals[i])&(zd<=xspvals[i+1])&(yd>=yspvals[j])&(yd<=yspvals[j+1])
		cut2 = (zr>=xspvals[i])&(zr<=xspvals[i+1])&(yr>=yspvals[j])&(yr<=yspvals[j+1])
		dcut=np.invert(cut1)
		rcut=np.invert(cut2)
		xdatsplit=xd[dcut]
		ydatsplit=yd[dcut]
		zdatsplit=zd[dcut]
		xrandsplit=xr[rcut]
		yrandsplit=yr[rcut]
		zrandsplit=zr[rcut]
		oX=np.concatenate((np.array([xdatsplit]).T,np.array([ydatsplit]).T,np.array([zdatsplit]).T),axis=1)
		rX=np.concatenate((np.array([xrandsplit]).T,np.array([yrandsplit]).T,np.array([zrandsplit]).T),axis=1)
		scatter(ydat,zdat,s=1, c='g',edgecolor='g')
		scatter(yrandsplit,zrandsplit,s=1,c='r', edgecolor='r')
		scatter(ydatsplit,zdatsplit,s=1,c='b')
		show()
		
		#TPCF=two_point(oX,radii,method='landy-szalay',data_R=rX, random_state=None)
		#xi=np.asarray(TPCF[0])
		#dd=np.asarray(TPCF[1])
		#dr=np.asarray(TPCF[2])
		#rr=np.asarray(TPCF[3])
		#RRjk.append(rr)
		#XIjk.append(xi)
		


'''

##Compute the variances (update to include the full covariance matrix)
c=0
C=[]
for i in range(len(XIjk)):
	print c
	sig = (np.sqrt(RRjk[i]/(1.0*RRfull))*(XIjk[i]-xifull))**2
	print sig
	c += sig
	if i+1 == len(XIjk):
		C.append(c)
	
stdev=C[0]**0.5



end2=time()
print 'Done Jackknife'
print 'C=', C
print 'stdev=', stdev
print end2-start2,'s'





file=open('SpIES_fullcorr_test_v3.txt','w')
file.write('s DD DR RR Xifull sigmajk\n')
for i in range(len(sep)):
	s= sep[i]
	dd= DDfull[i]
	dr= DRfull[i]
	rr=RRfull[i]
	xif=xifull[i]
	sigma=stdev[i]
	file.write(str(s)+' '+str(dd)+' '+str(dr)+' '+str(rr)+' '+str(xif)+' '+str(sigma)+'\n')
file.close()

f=open('SpIES_jackknife_results_v3.txt','w')
f.write('s RRjk Xijk \n')
for i in range(len(slopes)-1):
	for j in range(len(RRjk[i])):
		rrjk=RRjk[i][j]
		xijk=XIjk[i][j]
		se=sep[j]
		f.write(str(se)+' '+str(rrjk)+' '+str(xijk)+'\n')
f.close()







