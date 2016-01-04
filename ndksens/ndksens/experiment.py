"""docstring for experiment"""

import source
import mode
import units
import conflimits

import math



class Experiment(object):
    """An experimental search for a specific nucleon decay mode"""

    def __init__(self, name, source, mass, mode, eff, bgr):
        self.name = name ### Experiment's name
        self.source = source ### Source used by the experiment
        self.mass = mass ### Detector mass
        self.mode = mode ### Nucleon decay (NDK) mode searched for
        self.eff = eff ### Detection efficiency for a given NDK mode
        self.bgr = bgr ### Background rate for a given NDK mode

    def Nbkg(self, exposure):
        """Return the number of background events expected 
        for a given exposure."""
        return (self.bgr * exposure)

    def sensitivity(self, exposure, clc):
        """Return the experiment's partial lifetime sensitivity for this NDK decay mode."""
        aul = clc.AverageUpperLimit(self.Nbkg(exposure))
        return self.source.protonsperunitmass*self.eff*exposure/aul


############################################################
### DUNE NUCLEON DECAY SEARCHES --- DATABASE

### The pre-defined efficiency/background assumptions for an argon-based detector are taken from: A. Bueno et al., "Nucleon decay searches with large liquid argon TPC detectors at shallow depths: atmospheric and cosmogenic backgrounds," JHEP 0704 (2007) 041, http://arxiv.org/abs/hep-ph/0701101

source = source.Ar
mass = 40*units.kton

### dune_p1
name = 'dune_p1'
eff = 0.45
bgr = 1./(units.Mton*units.year)
dune_p1 = Experiment(name, source, mass, mode.p1, eff, bgr)

### dune_n1
name = 'dune_n1'
eff = 0.44
bgr = 8./(units.Mton*units.year)
dune_n1 = Experiment(name, source, mass, mode.n1, eff, bgr)

### dune_p16
name = 'dune_p16'
eff = 0.47
bgr = 2./(units.Mton*units.year)
dune_p16 = Experiment(name, source, mass, mode.p16, eff, bgr)

### dune_p19
name = 'dune_p19'
eff = 0.97
bgr = 1./(units.Mton*units.year)
dune_p19 = Experiment(name, source, mass, mode.p19, eff, bgr)

### dune_n34
name = 'dune_n34'
eff = 0.96
bgr = 2./(units.Mton*units.year)
dune_n34 = Experiment(name, source, mass, mode.n34, eff, bgr)

### dune_p41
name = 'dune_p41'
eff = 0.98
bgr = 1./(units.Mton*units.year)
dune_p41 = Experiment(name, source, mass, mode.p41, eff, bgr)


##############################
### DICTIONARY ###############

experiments = {dune_p1.name: dune_p1, dune_n1.name: dune_n1, 
               dune_p16.name: dune_p16, dune_p19.name: dune_p19, 
               dune_n34.name: dune_n34, dune_p41.name: dune_p41}

############################################################




if __name__ == '__main__':

    print "experiment.py executed directly"

    name = "dune_test"
    mass = 40.*units.kton 
    eff  = 1.
    bgr  = 1. /(units.Mton*units.year)

    dune_test = Experiment(name, source.Ar, mass, mode.n1, eff, bgr)

    FC = conflimits.FCMemoizer(90)
    FC.ReadTableAverageUpperLimits()
    sens = dune_test.sensitivity(400.*units.kton*units.year, FC)

    print "Sensitivity (at 90% CL) of the dune_test experiment: {0:.1e} years.".format(sens/units.year)

