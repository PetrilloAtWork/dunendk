"""docstring for module mode"""

import units
import math



class Mode(object):
    """docstring for Mode"""

    def __init__(self, symbol, name, pdgnumber, nucleonpdg, pdglimit): 
        super(Mode, self).__init__()
        self.symbol = symbol
        self.name = name
        self.pdgnumber = pdgnumber
        self.nucleonpdg = nucleonpdg
        self.pdglimit = pdglimit


############################################################
### NUCLEON DECAY MODES --- DATABASE

### The decay mode numbering scheme adopted matches the one of the Particle Data Group, see: http://pdg.lbl.gov/2015/listings/rpp2015-list-p.pdf. The existing limits are taken from the same source.

### ANTILEPTON + SINGLE MESON MODES

### p1
symbol = 'p1'
name = 'p2epluspizero'
pdglimit = 8200.E30*units.year
p1 = Mode(symbol, name, 1, 2212, pdglimit)

### n1
symbol = 'n1'
name = 'n2epluspiminus'
pdglimit = 2000.E30*units.year
n1 = Mode(symbol, name, 1, 2112, pdglimit)

### p16
symbol = 'p16'
name = 'p2mupluskzero'
pdglimit = 1600.E30*units.year
p16 = Mode(symbol, name, 16, 2212, pdglimit)

### p19
symbol = 'p19'
name = 'p2nubarkplus'
pdglimit = 5900.E30*units.year
p19 = Mode(symbol, name, 19, 2212, pdglimit)


### ANTILEPTON + MULTIPLE MESONS MODES

### LEPTON + SINGLE MESON MODES

### LEPTON + MULTIPLE MESONS MODES

### n34
symbol = 'n34'
name = 'n2eminuskplus'
pdglimit = 32.E30*units.year
n34 = Mode(symbol, name, 34, 2112, pdglimit)

### p41
symbol = 'p41'
name = 'p2muminuspipluskplus'
pdglimit = 245.E30*units.year 
p41 = Mode(symbol, name, 41, 2212, pdglimit)

##############################
### DICTIONARY ###############

modes = {p1.symbol: p1, n1.symbol: n1, 
         p16.symbol: p16, p19.symbol: p19, 
         n34.symbol: n34, p41.symbol: p41}

############################################################




if __name__ == '__main__':

    print "mode.py executed directly"
