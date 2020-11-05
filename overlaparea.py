import matplotlib.pyplot as plt
import numpy as np
import math
from datetime import datetime
from mpl_toolkits.mplot3d import Axes3D
from Functions import SimGrid
from sys import exit

starttime = datetime.now()
"""Define simulation base parameters, N=number of particles to simulate, Psize=
particle size, Pmass=particle mass in amu, EA=the electrochemically tested
area, Loadings=loadings in ng/cm^2 to simulate,
SC= stripcharge per area [muC/m^2]."""
N=2000
Psize = 3.8
Pmass = 370e3
EA = 3.14*(0.25)**2
SC= 340e4
Loading =  np.array([565, 997, 5000])
#np.linspace(10, 5000)
"""Calculate necessary parameters."""
Number_of_particles = Loading*6.022e14/Pmass*EA
Proj_Area_of_particle = math.pi*(Psize*1e-9/2)**2
Area_of_particle = 4*math.pi*(Psize*1e-9/2)**2
Charge_per_particle = SC*Area_of_particle*1e-6
Total_Theory_Charge = Number_of_particles*Charge_per_particle
SApar= 4*math.pi*(3.8/2)**2
SAtotal = 2000*SApar

"""The raster areas defined to be used to calculate the coverages."""
SmallRaster= 6.592e-6
LargeRaster= 3.127e-5

"""Selection of the raster area used."""
Dep_area = LargeRaster

"""Coverage calculation."""
Coverages =Number_of_particles*Proj_Area_of_particle/Dep_area*100
#Coverages = [0.7, 1.0, 2.7, 5.5, 1.0, 2.8, 5.5, 10.9, 10.8, 27.3, 27.4, 55.1, 275.1]


MR = []
BMR = []
BSTD = []
STD =[]
"""Simulate each of the coverages calculated above with 5 (for x in range(5))
repeated simulations of 2000 particles."""
for cov in Coverages:
    #print(cov)
    #Noverlaps = []
    Ratios = []
    Bratios = []
    #Overlapareas = []

    for repeat in range(5):
        Gsize = SimGrid(N, Psize, cov)
        """Generates the x,y positions for all N particles and sets their
        initial z coordinate to 1."""
        Positions = np.random.rand(2,N)*Gsize
        Z = np.full((N),1)

        for loop in range(8):
            for i, a in enumerate(Positions[0]):
                x, y = 0, 1
                Distance1 = np.sqrt((Positions[x, i] - Positions[x, :])**2 \
                 + (Positions[y, i] - Positions[y, :])**2 + (Z[i]-Z[:])**2)
                index = np.where(Distance1 < 35)
                #print(index[0])
                for g in index[0]:
                    if g !=i:
                        if Z[g]!=Z[i]:
                            Z[g]=Z[i]
                            Distance1 = np.sqrt((Positions[x, i] - \
                             Positions[x, :])**2 + (Positions[y, i] - \
                              Positions[y, :])**2 + (Z[i]-Z[:])**2)
                        Z[g] = Z[g]+ np.sqrt((Psize*10)**2 - Distance1[g]**2)

        counter=[]
        posx =np.array([])
        posy =np.array([])
        olapsim = []
        btomsim = []
        sumlosses = []
        Bases = []
        base = 0
        for i, a in enumerate(Positions[0]):
            olap =[]
            btom =[]
            x, y = 0, 1
            Distance1 = np.sqrt((Positions[x, i] - Positions[x, :])**2 + \
             (Positions[y, i] - Positions[y, :])**2 + (Z[i]-Z[:])**2)
            index = np.where(Distance1< Psize*10+5)
            index2 = np.where(Distance1<30)
            if Z[i] < 5:
                base += 1
                btom.append(12)
            for t in index2[0]:
                if t!=i:
                    print(Distance1[t])
            for b in index[0]:
                if b !=i:
                    #print(Distance1[b])
                    olap.append(8.4)
                    posx=np.append(posx,Positions[0,b])
                    posy=np.append(posy,Positions[1,b])

            olapsum = np.sum(olap)
            btomsum = np.sum(btom)
            sumloss = olapsum + btomsum
            if sumloss >= SApar:
                sumloss = SApar
            olapsim.append(olapsum)
            btomsim.append(btomsum)
            sumlosses.append(sumloss)
        #print(base)
        Bases.append(base)
        #print(np.sum(sumlosses))
        #print(len(counter))
        Area = np.sum(sumlosses)
#        Bareas = np.sum(btomsim)
        Ratio = (SAtotal-Area)/SAtotal
#        Bratio = (SAtotal-Bareas)/SAtotal
        Ratios.append(Ratio)
#        Bratios.append(Bratio)
    #print(np.mean(Bases))
    MR.append(np.mean(Ratios))
    STD.append(np.std(Ratios))
#    BMR.append(np.mean(Bratios))
#    BSTD.append(np.std(Bratios))
#print(datetime.now()-starttime)
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(Positions[0], Positions[1], Z, c='r', marker='o')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()
TheoryCharge = [6.5e-6, 9.4e-6, 2.4e-5, 4.9e-5, 8.6e-6, 2.5e-5, 4.9e-5, 9.8e-5, 9.7e-5, 2.44e-4, 2.5e-4, 4.9e-4, 0.00246]
RealCharge = np.array([3.5e-6, 6e-6, 1.7e-5, 2.9e-5, 5.4e-6, 1.8e-5, 3.1e-5, 6.8e-5, 5e-5, 1.6e-4, 1.5e-4, 2.9e-4, 1.1e-3])/10
Loadingmes = [13, 19, 48, 100, 17, 51, 100, 199, 197, 495, 498, 1002, 5000]
#BT = []
#BST = []
T = []
ST = []
for i, a in enumerate(Total_Theory_Charge):
    T.append(a * MR[i])
    ST.append(a*STD[i])
#    BT.append(a*BMR[i])
#    BST.append(a*BSTD[i])
T = np.array(T)
ST =np.array(ST)
np.savetxt('SimCharge.txt', T)
np.savetxt('SimLoad.txt', Loading)
#BT = np.array(BT)
#BST =np.array(BST)
#T-(T[0]-RealCharge[0])
plt.figure(2)
plt.errorbar(Loading, T, ST)
#plt.errorbar(Loading, BT, BST)
#plt.scatter(Loading, T, color = 'red')
plt.scatter(Loadingmes, RealCharge, color = 'k')
plt.xlabel('Loading [ng]', Fontsize = 18)
plt.ylabel('CO Strip Charge [C]', Fontsize = 18)
#plt.xlim([12, 50])
#plt.ylim([0,0.000026])
#Axes.ticklabel_format(self, axis='both', style='sci')
#plt.scatter(Loading, TheoryCharge, color = 'blue')
plt.show()
#plt.errorbar(Coverages, CorrAs/TotalApar, STD/TotalApar)


#plt.errorbar([500, 1000, 5000], MeanR, Stddev,fmt='o', capsize = 3)
#plt.xlabel('Loading [ng]', Fontsize = 18)
#plt.ylabel('Overlapping Particles [%]', Fontsize = 18)
#plt.show()
