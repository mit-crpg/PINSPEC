import matplotlib.pyplot as plt
from isotope import *

def main():
    
    print 'Starting PIN SPEC...'

    # Define isotopes
    h1 = Isoptope('H-1')
    h1.setAO(2.0)
    o16 = Isotope('O-16')
    o16.setAO(1.0)

    u235 = Isoptope('U-235')
    U235.setAO(3.0)
    u238 = Isotope('U-238')
    U238.setAO(97.0)
    
    # Define materials
    #moderator = Material()
    #moderator.setDensity(1.0, 'g/cc')
    #moderator.addIsotope(h1)
    #moderator.addIsotope(o16)

    #fuel = Material()
    #fuel.setDensity(10.0, 'g/cc')
    #fuel.addIsotope(u235)
    #fuel.addIsotope(u238)
    
    #moderator.loadXS()
    #fuel.loadXS()
    
    # Define regions
    #region_mod = Region()
    #region_mod.addMaterial(moderator)
    #region_mod.setType('moderator')

    #region_fuel = Region()
    #region_fuel.addMaterial(fuel)
    #region_fuel.setType('fuel')

    # Define geometry
    #geometry = Geometry()
    #geometry.addRegion(region_mod)
    #geometry.addRegion(region_fuel)
    
    # Set energy bins
    

    # Specify data to output
    
    
    # Run code



if __name__ == '__main__':
    
    main()
