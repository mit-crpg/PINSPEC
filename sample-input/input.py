import matplotlib.pyplot as plt
from pinspec import *

def main():
    
    print 'Starting PIN SPEC...'

    # Define isotopes
    h1 = Isotope('H-1')
    o16 = Isotope('O-16')
    u235 = Isotope('U-235')
    u238 = Isotope('U-238')
    
    print 'Made isotopes'
    
    # Define materials
    moderator = Material()
    moderator.setDensity(1.0, 'g/cc')
    moderator.addIsotope(h1, 2.0)
    moderator.addIsotope(o16, 1.0)
    
    print 'Added isotopes to materials'
    
    fuel = Material()
    fuel.setDensity(10.0, 'g/cc')
    fuel.addIsotope(u235, 0.02)
    fuel.addIsotope(u238, 0.98)
    fuel.addIsotope(o16, 2.0)
    
    print 'Made materials'

    # Define regions
    region_mod = Region()
    region_mod.setMaterial(moderator)
    region_mod.setRegionType(MODERATOR)
    region_mod.setRegionName('moderator')
    
    print 'Set region_mod type'

    region_fuel = Region()
    region_fuel.setMaterial(fuel)
    region_fuel.setRegionType(FUEL)
    region_fuel.setRegionName('fuel')
    
    print 'Made regions'

    # Define geometry
    geometry = Geometry()
    geometry.addRegion(region_mod)
    geometry.addRegion(region_fuel)

    print 'Made geometry'

    del geometry

    # Set energy bins
    

    # Specify data to output
    
    
    # Run code



if __name__ == '__main__':
    
    main()







    def __iter__(self):                 # ???
        for r in self.regions:
            yield r
