surface = ZCylinder('testing zcylinders')
surface.setX0(0.0)
surface.setY0(0.0)
surface.setRadius(0.5)

neutron = createNewNeutron()

x = []
y = []

# Send neutron traveling towards the right
neutron._u = random.random()
neutron._v = random.random()
neutron._w = random.random()

norm = math.sqrt(neutron._u**2 + neutron._v**2 + neutron._w**2)
neutron._mu = math.cos(neutron._w / norm)
neutron._phi = math.atan2(neutron._v, neutron._u)
neutron._u /= norm
neutron._v /= norm
neutron._w /= norm

theta = math.acos(neutron._mu)

for i in range(10):
    
    neutron._x = 0.5*random.random() - 0.25
    neutron._y = 0.5*random.random() - 0.25
    neutron._z = random.random()

    dist = surface.computeParametrizedDistance(neutron)
    print 'dist = ' + str(dist)

    if dist != float('inf'):
        neutron._x += dist * neutron._u
        neutron._y += dist * neutron._v
        neutron._z += dist * neutron._w

    print 'new x = ' + str(neutron._x) + ' y = ' + str(neutron._y) + ' z =' + \
        str(neutron._z)

    dist = surface.computeParametrizedDistance(neutron)
    print 'dist = ' + str(dist)

#    geometry.initializeSourceNeutron(neutron)

#    region = neutron._region
#    region.collideNeutron(neutron)

#    if region_fuel.onBoundary(neutron):
#        x.append(neutron._x)
#        y.append(neutron._y)


#plt.scatter(x,y)
#plt.plot(x,y,'wo')
#plt.show()
