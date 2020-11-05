import math
import numpy as np

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
