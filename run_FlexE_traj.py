#! /usr/bin/env python

from prody import *
import hamiltonian
import argparse
import MDAnalysis
import numpy


def load_trajectory(top,dcd):
    """
    """
    universe = MDAnalysis.Universe(top, dcd)
    calpha = universe.selectAtoms('name CA')
    ca_coords = universe.trajectory.timeseries(calpha,format='fac')
    ref_coords = calpha.coordinates()

    ensemble = prody.dynamics.Ensemble('MD snapshots')
    ensemble.addCoordset(ca_coords)
    return ensemble

def main():
    #no log messages:
    #prody.ProDySetVerbosity('none')
    changeVerbosity('none')
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Calculate MDENM energies from a pdb \
                                    will calculate energy using modes from pdb\
                                    and then from reference--> crystal should\
                                    be the reference')
    parser.add_argument('--pdb', help='Molecule we want to examine.')
    parser.add_argument('--top', help='It will be used as topology')
    parser.add_argument('--dcd', help='File containing trajectory coordinates.')
    parser.add_argument('--ref', help='compare trajectory to pdb, if none, compare trajectory with itself',default='')
    args = parser.parse_args() 

    #Load the structures
    pdb = parsePDB(args.pdb)
    calphas = pdb.select('calpha')

    if args.ref:
        ref = parsePDB(args.ref)
        ref_alpha = ref.select('calpha')

    ensemble = load_trajectory(args.top,args.dcd)

    ncoor = ensemble.getNumOfCoordsets()

    if args.ref:
        native = ref_alpha.copy()
        h = hamiltonian.EDENMHamiltonian( native.getCoordinates() )
        for pred in ensemble.iterCoordsets():
            Forw_E_ED = h.evaluate_energy( pred.getCoordinates())
            h = hamiltonian.EDENMHamiltonian( pred.getCoordinates() )
            Back_E_ED = h.evaluate_energy( native.getCoordinates())
            print Forw_E_ED,Back_E_ED
    else:
        energies = numpy.zeros((ncoor,ncoor))
        maximum = numpy.zeros((ncoor,ncoor))
        for i in range(0,ncoor):
            native = ensemble.getCoordsets(indices=i)
            h = hamiltonian.EDENMHamiltonian( native )
            for j in range(i+1,ncoor):
                pred = ensemble.getCoordsets(indices=j)
                Forw_E_ED = h.evaluate_energy( pred)
                h = hamiltonian.EDENMHamiltonian( pred )
                Back_E_ED = h.evaluate_energy( native)
                energies[i,j] = (Forw_E_ED+Back_E_ED)/2
                energies[j,i] = energies[i,j]
		max_number = Forw_E_ED
		if (Back_E_ED > max_number):
			max_number = Back_E_ED
		maximum[i,j] = max_number
		maximum[j,i] = max_number
        numpy.savetxt('Flex_matrix.txt',energies,fmt='%8.2f')
        numpy.savetxt('Flex_max.txt',maximum,fmt='%8.2f')



if __name__ == '__main__':
    main()

