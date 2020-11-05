import matplotlib.pyplot as plt
import numpy as np
import math

def BoundaryFilter(Positions, Gridsize, Boundary):
    """ This function returns two lists of x-y particle coordinates for
    particles that were more than the Boundary value (in nanometers, Bvalue
    defined below converts to Ångstrøm) from the edge of the grid. The
    "Positions" parameter is a 2D array of x-y particle positions and the
    gridsize is the sidelength of the simulation grid."""

    """Creating empty lists for appending later"""
    dim1 =[]
    dim2 =[]
    """Converting the desired boundary value from nanometers to Ångstrøm"""
    Bvalue = Boundary*10

    """Loop through each position and determine if it is filtered or not."""
    for i, a in enumerate(Positions[0]):
        if Positions[0][i]>Bvalue and Positions[0][i]<Gridsize-Bvalue and \
         Positions[1][i]>Bvalue and Positions[1][i]<Gridsize-Bvalue:
            dim1.append(Positions[0][i])
            dim2.append(Positions[1][i])
    """Returns the two lists of filtered coordinates."""
    return (dim1, dim2)

def FindDistances(FiltDim1, FiltDim2, Unfiltered, Npar):
    """Get the length of X, corresponding to the number of filtered particles
    and create two 2D arrays filled with zeros for distances in the x and
    y direction."""
    D=len(FiltDim1)
    Distance1 = np.zeros((D,Npar))
    Distance2 = np.zeros((D,Npar))

    """Find distances in the x-direction from each filtered particle to all
    other particles (also the unfiltered particles, otherwise the boundary
    effect is not avoided)."""
    for i, pos1 in enumerate(FiltDim1):
        for d ,pos2 in enumerate(Unfiltered[0]):
            Distance1[i][d]=(pos1-pos2)

    for i, pos1 in enumerate(FiltDim2):
        for d ,pos2 in enumerate(Unfiltered[1]):
            Distance2[i][d]=(pos1-pos2)

    """Creates a single 2D array where pythagoras is used to calculate the
    actual interparticle distances from the distances in the x and
    y directions."""
    Distance = np.sqrt(Distance1**2+Distance2**2)
    return Distance

def SimGrid(NumPar, ParSize, SimCoverage):
    """Calculation of the simulation area from the desired coverage and area
     per particle."""
    TotalApar = math.pi*(ParSize/2)**2*NumPar
    Area = TotalApar/(SimCoverage/100)
    """The grid size is set in units of an Å to increase the accuracy of the
     simulation."""
    Asize = np.sqrt(Area)
    Gridsize = round(Asize)*10
    return Gridsize

def SmallestDistNonNeg(Distances, ParSize):
    """Runs through each row in Distances and finds the smallest nonzero distance.
    Each row contains the distances from one particle to all others. Thus it finds
    the smallest interparticle distance for all particles."""
    SmallDist = []
    for a in Distances:
        SmallDist.append(np.min(a[np.nonzero(a)]))

    """Convert distances back to nanometers."""
    SmallDist = np.array(SmallDist)
    SmallDist = SmallDist/10

    """Find the edge to edge interparticle distance instead of center to center."""
    SmallDist = SmallDist-ParSize

    """Copy the absolute distances before removing negatives. Not used for anything
    currently, but can be interesting for other purposes."""
    SmallDistNeg=SmallDist

    """Convert to numpy array and set all distances smaller than 0 to 0 since for
    the current purpose any particle overlapping is the same."""
    SmallDist[SmallDist < 0] = 0

    return SmallDist



"""Projected coverage in % monolayer."""
coverage = 0.6
""" Particle diameter in nm. """
Psize = 3.8
"""Number of particles to include in the simulation."""
N=2000
Gsize = SimGrid(N, Psize, coverage)

""" Creates a loop around the simulation that allows for multiple simulations
of 2000 particles"""
AvgDists=[]
StdDevs=[]
TotalSmallDists=[]
for i in range(5):

    """A 2D array of random x-y coordinates for N particles is generated."""
    Positions = np.random.rand(2,N)*Gsize

    """The positions are filtered to avoid boundary effects. Any positions
    that lie more than 30 nm (300 Å) from the edge of the simulated area are
    saved into the X-Y lists. These positions are the ones that will be used
    to find the mean interparticle distance."""
    (X,Y) = BoundaryFilter(Positions, Gsize, 30)

    """Finds the distances between the particles and returns them in a
     2D array."""
    Distances = FindDistances(X, Y, Positions, N)

    """Finds the smallest distances for each particle and sets all distances
     below zero to zero."""
    SmallDist = SmallestDistNonNeg(Distances, Psize)

    """Find the mean interparticle distance and the standard deviation.
    Print them out in the terminal."""
    AverageInterparticleDistance = np.mean(SmallDist)
    Std = np.std(SmallDist)
    print(AverageInterparticleDistance, Std)

    AvgDists.append(AverageInterparticleDistance)
    StdDevs.append(Std)
    TotalSmallDists.append(SmallDist)



TotalSmallDistslist= [item for sublist in TotalSmallDists for item in sublist]

MeanIPD= np.mean(AvgDists)
IPDStdDev= np.std(AvgDists)
MeanDistStdDev= np.mean(StdDevs)
TotalNpar=len(TotalSmallDistslist)

print('Number of Particles included:', TotalNpar)
print( 'The mean interparticle distance is', MeanIPD, \
'nm. With a standard deviation of', IPDStdDev,'nm.')
print(' The distribution has a standard deviation of', MeanDistStdDev, 'nm.')

"""Creates a histogram of the interparticle distance distribution.
Note that the bar at 0 may be very large at higher coverages since all negative
interparticle distances are set to 0."""
plt.hist(TotalSmallDistslist,30, edgecolor = 'black')
#plt.xlim(0,50)
#plt.ylim(0,160)
plt.xlabel('Interparticle Distance [nm]', Fontsize = 18)
plt.ylabel('Counts', Fontsize = 18)
plt.tick_params(labelsize=14)
plt.tight_layout()


"""Saves the figure with the given name so it can be used elsewhere."""
savename = 'SimofSize{}nmCoverage{}.png'.format(Psize,coverage)
plt.savefig(savename, dpi= 1200, format = 'png')
plt.show()


"""Gives the possibility to plot the particle positions so that one can see
them."""
#plt.figure()
#plt.scatter(Positions[0], Positions[1])
#plt.show()
