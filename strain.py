from pymatgen.io.vasp import Poscar
from pymatgen.core import Structure
import numpy as np 

#SETTINGS:
scalelength=0.4 #scale length of vectors to better visualise
minmom=0.1 #minimum magnetic moment to show
vecradius=0.5 #vector radius
veccolour=[255,0,0] #vector RGB colour

#Read in ext defect supercell and work out matrix to transform vectors
gb = Structure.from_file(filename="./POSCAR-GB.vasp")

#Read in grain ref structure
grain = Structure.from_file(filename="./POSCAR-grain.vasp")

#x-coordinates defining range of grain to be analysed
Gmin=0.0*gb.lattice.matrix[2][2]
Gmax=0.5*gb.lattice.matrix[2][2]

#Atom ids to be aligned in GB and grain structures (note these are ids that start at 1)
GrainRefSite=234 
GBRefSite=54
MaxStrain=1.0

#extract sites from structures and shift so coincident in centre of grain
gbsites=list(gb.sites)
p1=gbsites[GBRefSite-1].coords
grainsites=list(grain.sites)
p2=grainsites[GrainRefSite-1].coords
shift=p1-p2

for site in grainsites:
    site.coords=site.coords+shift

#Sort GB sites according to z coord: 
sortedgbsites = sorted(gbsites, key=lambda site: site.coords[2])

#Calculate atom displacements wrt bulk ref within specific grain
for i, site in enumerate(sortedgbsites):
    if site.coords[2] > Gmin and site.coords[2] < Gmax:
        for siteref in grainsites:
          delta=site.coords-siteref.coords
          dist=np.linalg.norm(delta)
          if dist<MaxStrain:
              print(site.coords[2],dist)
              #print(i,delta[0],delta[1],delta[2])



