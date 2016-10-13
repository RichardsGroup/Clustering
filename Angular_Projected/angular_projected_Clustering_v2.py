import os
from time import time
import astropy.io.fits as pf
import numpy as np
from numpy import sin, cos
from numpy.linalg import norm
from astroML.correlation import two_point, bootstrap_two_point
from multiprocessing import Pool
from functools import partial
import matplotlib.pyplot as plt
import random as rd

# Find the separation angle between quasars
# Calculate pi = |s1-s2|cos(theta/2). For small angles (<2 degrees) cos(theta/2)=1
# Calculate r_p = (s1+s2) * sin(theta/2). For small angles (<2 degrees) sin(theta/2) = theta/2

#Define a function to calculate rp and pi for the data
def data_loop(N,datvect,thetarad):
	DD = np.zeros((len(thetarad)))
	for i in range(N):
		if i % 1000==0:
			print i/float(len(datvect))  #Step counter
		'''
		vect = SkyCoord([datvect[i]],unit='deg')
		rest = SkyCoord(datvect,unit='deg')
		ang = vect.separation(rest)
		angle = np.asarray(ang.degree)
		'''
		
		vect = datvect[i] * np.pi/180.0  #The x,y,z cartesian vector of a quasar in the data wrt the origin [RA1,DEC1]
		rest = datvect * np.pi/180.0 #The remaining array of vectors in the data set (i.e.[[RA2,DEC2],[RA3,DEC3],...
		angle = np.arccos((sin(vect[1])*sin(rest[:,1]))+(cos(vect[1])*cos(rest[:,1])*cos(vect[0]-rest[:,0])))*180*60/np.pi #angle between vector 1 and the rest in arcmin
		#Build a DD array with increasing theta
		for j in range(len(thetarad)):
			gdx = (angle < thetarad[j]) & (np.abs(angle) != 0.0)
			DD[j] += float(len(angle[gdx]))
	dd = np.diff(DD)	
		
	return dd
	
#Define function to calculate rp and pi for randoms (See function above for comments)
def rand_loop(thetarad,randvect,count):
	RR = np.zeros(len(thetarad))
	for i in count:
		s = time()
		'''
		vect = SkyCoord([randvect[i]],unit='deg')
		rest = SkyCoord(randvect,unit='deg')
		ang = vect.separation(rest)
		angle = np.asarray(ang.degree)
		'''

		vect = randvect[i] * np.pi/180.0
		rest = randvect * np.pi/180.0
		angle = np.arccos((sin(vect[1])*sin(rest[:,1]))+(cos(vect[1])*cos(rest[:,1])*cos(vect[0]-rest[:,0])))*180*60/np.pi
		
		#Build an RR matrix with increasing rp across the matrix and increasing pi down the matrix
		for j in range(len(thetarad)):
			gdx = (angle < thetarad[j]) & (np.abs(angle) != 0.0)
			RR[j] += float(len(angle[gdx]))
	
		e = time()
		print "RR loop",i,"done in",e-s,"s"
	rr = np.diff(RR)
	return rr
	

	
#Define function to calculate rp and pi for the data - random pairs (See function above for comments)
def dat_rand_loop(thetarad,datvect,randvect,count):
	DR = np.zeros(len(thetarad))
	for i in count:
		s = time()
		'''
		vect = SkyCoord([datvect[i]],unit='deg')
		rest = SkyCoord(randvect,unit='deg')
		ang = vect.separation(rest)
		angle = np.asarray(ang.degree)
		'''

		vect = datvect[i] * np.pi/180.0
		rest = randvect * np.pi/180.0
		angle = np.arccos((sin(vect[1])*sin(rest[:,1]))+(cos(vect[1])*cos(rest[:,1])*cos(vect[0]-rest[:,0])))*180*60/np.pi
		
		#Build an DR matrix with increasing rp across the matrix and increasing pi down the matrix
		for j in range(len(thetarad)):
			gdx = (angle < thetarad[j])
			DR[j] += float(len(angle[gdx]))
	
		e = time()
		print "DR loop",i,"done in",e-s,"s"
	dr = np.diff(DR)
	return dr

############################################          #         #
#################################################################
############################################          #         #

#Open the positions here

start1=time()

#Home Path
#file1='/home/john/Clustering/Combine_SpIES_Shela/Data_sets/Match_SpSh_Cand_Spec_wz.fits'
#file2='/home/john/Clustering/Combine_SpIES_Shela/Data_sets/SpSh_randoms.fits'

#Triton path
file1='/Users/johntimlin/Clustering/Combine_SpIES_Shela/Data_sets/Match_SpSh_Cand_Spec_wzcorrected_nooutlier.fits'
file2='/Users/johntimlin/Clustering/Combine_SpIES_Shela/Data_sets/SpSh_randoms.fits'

#Dirac Path
#file1='/home/jtimlin/CLUSTERING_IN_SPIES/high_z_20160615/SpIES_SDSS_Specz20150619_Photz20160612.fits'
#file2='/home/jtimlin/CLUSTERING_IN_SPIES/high_z_20160615/SpIES_random_specplphot_highz.fits'

data=pf.open(file1)[1].data
data2=pf.open(file2)[1].data
# Separate the data and random points and stack so that I have an array of [x,y,z] coordinates for each quasar


dx = ((data.ZSPEC>=2.9) | (data.zphotbest>=2.9)) & ((data.ZSPEC<=5) | (data.zphotbest<=5))

datra = data.RA[dx]
datdec = data.DEC[dx]

datvect = np.concatenate((np.array([datra]).T,np.array([datdec]).T),axis = 1)

randra = data2.RA[:]
randdec = data2.DEC[:]

randvect = np.concatenate((np.array([randra]).T,np.array([randdec]).T),axis = 1)

print len(randvect),len(datvect)

#########
#Make array for theta radii (in ARCMIN)

th=np.linspace(np.log10(0.04),np.log10(250),22)
thetarad=10**th



#########
#Run the data conversion loop

datapairs = data_loop(len(datvect),datvect,thetarad)

DDpairs = datapairs  #Returned DD matrix

end1= time()

print "Data conversion done in",end1-start1,"s"

start2=time()


########################################################################
#Begin multiprocess of randoms

nproc = 4  #Number of processors on which to run
vals = np.arange(len(randvect))  #an array of integers the length of the data to split for multiprocessing
ev=len(vals)/nproc #Define an equal split of data...NOTE!!!! THIS MUST BE INT DIVISION FOR THIS TO WORK!
sploc=[]
for i in range(nproc-1):
	sploc.append((i+1)*ev)
ct = np.split(vals,sploc)


#Call the Pool command to multiprocess
if __name__ == '__main__':
	p = Pool(nproc)
	func = partial(rand_loop,thetarad,randvect)
	randoms = p.map(func, ct) #Map the random_loop routine to each processor with a single array in ct to run over
	RRdummy = np.zeros((len(thetarad)-1))
	print randoms
	print np.shape(RRdummy)
	for i in range(len(randoms)):
		RRdummy = np.add(RRdummy,randoms[i])
	RRpairs = RRdummy

end2=time()
print "Random conversion done in",end2-start2,"s"


##################################################
#Begin multiprocess of data-randoms
start3=time()
nproc2 = 4  #Number of processors on which to run
vals2 = np.arange(len(datvect))  #an array of integers the length of the data to split for multiprocessing
c=len(vals2)/nproc2 #Define an equal split of data...NOTE!!!! THIS MUST BE INT DIVISION FOR THIS TO WORK!
array=[]
for i in range(nproc2-1):
	array.append((i+1)*c)
ct2 = np.split(vals2,array)


#Call the Pool command to multiprocess
if __name__ == '__main__':
	p = Pool(nproc2)
	funcdr = partial(dat_rand_loop,thetarad,datvect,randvect)
	datrands = p.map(funcdr, ct2) #Map the data-random_loop routine to each processor with a single array in ct to run over
	DRdummy = np.zeros((len(thetarad)-1))
	for i in range(len(datrands)):
		DRdummy = np.add(DRdummy,datrands[i])
	DRpairs = DRdummy
	
end3=time()

print "Data Random conversion done in",end3-start3,"s"
start4=time()



#Calculate XI matrix 

factor = len(randvect) * 1.0/len(datvect)
print len(randvect)
print len(datvect)
xi = (pow(factor,2)*DDpairs - (2 * factor * DRpairs) + RRpairs) / RRpairs

for i in range(len(thetarad)-1):
	DD = DDpairs[i]
	DR = DRpairs[i]
	RR = RRpairs[i]
	XI = xi[i]
	#Write/Append to file depending on whether file exists or not
	if os.path.exists("SpSh_angular_Clustering_v2.txt"):
		infile=open("SpSh_angular_Clustering_v2.txt" ,"a")
		THavg1 = (thetarad[i] + thetarad[i+1])/2.0
		infile.write(str(THavg1)+ ' '+ str(DD)+ ' ' + str(DR) + ' ' + str(RR) + ' '+str(XI)+'\n')
		infile.close()
	else:
		infile=open("SpSh_angular_Clustering_v2.txt" ,"w")
		infile.write("#Theta DD DR RR XI factor="+str(factor)+" \n")
		THavg1 = (thetarad[i] + thetarad[i+1])/2.0
		infile.write(str(THavg1)+ ' '+ str(DD)+ ' ' + str(DR) + ' ' + str(RR) + ' '+str(XI)+'\n')
		infile.close()
    		

end4=time()

print 'Full 2D 2PCF calculated in',end4-start4,"s"


print 'Begin Jackknife'

jknum = 10
D = datvect[datvect[:,0].argsort()]
jksplit = np.arange(len(datvect))  #an array of integers the length of the data to split for multiprocessing
evensplit=len(jksplit)/jknum #Define an equal split of data...NOTE!!!! THIS MUST BE INT DIVISION FOR THIS TO WORK!
splitloc=[]
for i in range(jknum-1):
	splitloc.append((i+1)*evensplit)
jklist = np.split(D,splitloc)

R =randvect[randvect[:,0].argsort()]

print jklist[0][:,0]
print R[:,0][(R[:,0] > 0.0) & (R[:,0] < max(jklist[0][:,0]))]

jkrand = []
jkdat = []
for i in range(len(jklist)):
	if i == 0:
		cutd = (D[:,0] > 0.0) & (D[:,0] <= min(jklist[i+1][:,0]))
		cutr = (R[:,0] > 0.0) & (R[:,0] <= min(jklist[i+1][:,0]))
		dcut = np.invert(cutd)
		rcut = np.invert(cutr)
		darr = D[dcut]
		rarr = R[rcut]
		jkdat.append(darr)
		jkrand.append(rarr)
		print 'first', i
		print len(jkrand[i])/float(len(randvect))
		print len(jkdat[i])/float(len(datvect))
	elif i == jknum-1:
		cutd = (D[:,0] >= min(jklist[i][:,0])) & (D[:,0] < 360.0)
		cutr = (R[:,0] >= min(jklist[i][:,0])) & (R[:,0] < 360.0)
		dcut = np.invert(cutd)
		rcut = np.invert(cutr)
		darr = D[dcut]
		rarr = R[rcut]
		jkdat.append(darr)
		jkrand.append(rarr)
		print 'last', i
		print len(jkrand[i])/float(len(randvect))
		print len(jkdat[i])/float(len(datvect))
	else:
		cutd = (D[:,0] > max(jklist[i-1][:,0])) & (D[:,0] <= min(jklist[i+1][:,0]))
		cutr = (R[:,0] > max(jklist[i-1][:,0])) & (R[:,0] <= min(jklist[i+1][:,0]))
		dcut = np.invert(cutd)
		rcut = np.invert(cutr)
		darr = D[dcut]
		rarr = R[rcut]
		jkdat.append(darr)
		jkrand.append(rarr)
		print 'mid', i
		print len(jkrand[i])/float(len(randvect))
		print len(jklist[i])/float(len(datvect))

### Plot conditional as a test
'''
for i in range(jknum):
	#plt.subplot(111,projection = 'aitoff')
	plt.scatter(jkrand[i][:,0],jkrand[i][:,1],s = 1,color='g')
	plt.scatter(jkdat[i][:,0],jkdat[i][:,1],s = 1,color='r')
	plt.show()
'''

for i in range(jknum):

	dvect = jkdat[i]
	rvect = jkrand[i]

	#########Run the data conversion loop

	datapairs = data_loop(len(dvect),dvect,thetarad)

	DDpairs = datapairs  #Returned DD matrix

	end1= time()

	print "Data conversion done in",end1-start1,"s"

	start2=time()


	########################################################################
	#Begin multiprocess of randoms

	nproc = 4  #Number of processors on which to run
	vals = np.arange(len(rvect))  #an array of integers the length of the data to split for multiprocessing
	ev=len(vals)/nproc #Define an equal split of data...NOTE!!!! THIS MUST BE INT DIVISION FOR THIS TO WORK!
	sploc=[]
	for i in range(nproc-1):
		sploc.append((i+1)*ev)
	ct = np.split(vals,sploc)
	
	
	#Call the Pool command to multiprocess
	if __name__ == '__main__':
		p = Pool(nproc)
		func = partial(rand_loop,thetarad,rvect)
		randoms = p.map(func, ct) #Map the random_loop routine to each processor with a single array in ct to run over
		RRdummy = np.zeros((len(thetarad)-1))
		print randoms
		print np.shape(RRdummy)
		for i in range(len(randoms)):
			RRdummy = np.add(RRdummy,randoms[i])
		RRpairs = RRdummy
		
	end2=time()
	print "Random conversion done in",end2-start2,"s"


	##################################################
	#Begin multiprocess of data-randoms
	start3=time()
	nproc2 = 4  #Number of processors on which to run
	vals2 = np.arange(len(dvect))  #an array of integers the length of the data to split for multiprocessing
	c=len(vals2)/nproc2 #Define an equal split of data...NOTE!!!! THIS MUST BE INT DIVISION FOR THIS TO WORK!
	array=[]
	for i in range(nproc2-1):
		array.append((i+1)*c)
	ct2 = np.split(vals2,array)


	#Call the Pool command to multiprocess
	if __name__ == '__main__':
		p = Pool(nproc2)
		funcdr = partial(dat_rand_loop,thetarad,dvect,rvect)
		datrands = p.map(funcdr, ct2) #Map the data-random_loop routine to each processor with a single array in ct to run over
		DRdummy = np.zeros((len(thetarad)-1))
		for i in range(len(datrands)):
			DRdummy = np.add(DRdummy,datrands[i])
		DRpairs = DRdummy
	
	end3=time()

	print "Data Random conversion done in",end3-start3,"s"
	start4=time()



	#Calculate XI matrix 
	
	factor = len(randvect) * 1.0/len(datvect)
	print len(randvect)
	print len(datvect)
	xi = (pow(factor,2)*DDpairs - (2 * factor * DRpairs) + RRpairs) / RRpairs
	
	for i in range(len(thetarad)-1):
		DD = DDpairs[i]
		DR = DRpairs[i]
		RR = RRpairs[i]
		XI = xi[i]
		#Write/Append to file depending on whether file exists or not
		if os.path.exists("SpSh_angular_Clustering_Jackknifedat_v2.txt"):
			infile=open("SpSh_angular_Clustering_Jackknifedat_v2.txt" ,"a")
			THavg1 = (thetarad[i] + thetarad[i+1])/2.0
			infile.write(str(THavg1)+ ' '+ str(DD)+ ' ' + str(DR) + ' ' + str(RR) + ' '+str(XI)+'\n')
			infile.close()
		else:
			infile=open("SpSh_angular_Clustering_Jackknifedat_v2.txt" ,"w")
			infile.write("#Theta DD DR RR XI factor="+str(factor)+" "+"jknum="+" "+str(jknum)+ "\n")
			THavg1 = (thetarad[i] + thetarad[i+1])/2.0
			infile.write(str(THavg1)+ ' '+ str(DD)+ ' ' + str(DR) + ' ' + str(RR) + ' '+str(XI)+'\n')
			infile.close()
    
	infile=open("SpSh_angular_Clustering_Jackknifedat_v2.txt" ,"a")
	infile.write("STEP_DONE \n")
	infile.close()
	end5=time()
	
end5=time()
print "Total Time=",end5-start1,"s"

	
	
	
	
	
	
	
	
