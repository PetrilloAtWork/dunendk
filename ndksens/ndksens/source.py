"""docstring for module source"""

import units
import math



class Source(object):
    """docstring for Source"""

    def __init__(self, symbol, protonsperunitmass, neutronsperunitmassA):
        super(Source, self).__init__()
        self.symbol = symbol
        self.protonsperunitmass = protonsperunitmass # protons per unit mass
        self.neutronsperunitmass = neutronsperunitmass # neutrons per unit mass


############################################################
### SOURCES --- DATABASE

### Ar
symbol = 'Ar'
protonsperunitmass = 2.71E32*(1/units.kton)
neutronsperunitmass = 3.31E32*(1/units.kton)
Ar = Source(symbol, protonsperunitmass, neutronsperunitmass)

### H20 
symbol = 'H20'
protonsperunitmass = 3.34E32*(1/units.kton)
neutronsperunitmass = 2.68E32*(1/units.kton)
H20 = Source(symbol, protonsperunitmass, neutronsperunitmass)



##############################
### DICTIONARY ###############

sources = {Ar.symbol: Ar, H20.symbol: H20}

############################################################


if __name__ == '__main__':

    print "source.py executed directly"
    print "Protons per kton of argon: {0:.1e}".format(Ar.protonsperunitmass/(1/units.kton))


