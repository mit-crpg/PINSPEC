#!/bin/bash -ex

CFLAGS="-O3 -march=native -fPIC -I/usr/include/python2.7"
LFLAGS="-bundle `python-config --ldflags`"

if [ "$1" = '--release' ]
then
    CFLAGS="$CFLAGS"
elif [ "$1" = "--debug" ]
then
    CFLAGS="-lf2c -lm -DDEBUG -g"
elif [ "$1" = '--profile' ]
then
    CFLAGS="$CFLAGS -pg -g"
elif [ "$1" = '--benchmark' ]
then
    CFLAGS="$CFLAGS -DBENCHMARK"
fi


if [ ! -d "objs" ]; then
	mkdir objs
fi


# Comment and uncomment the std::map container in Isotope as a
# workaround for SWIG
cp src/Isotope.h src/Isotope.bak
sed 's/	std::map<collisionType/\/\/	std::map<collisionType/' src/Isotope.bak > src/Isotope.h
swig -python -c++ python/pinspec/Geometry.i
cp src/Isotope.h src/Isotope.bak
sed 's/\/\/	std::map<collisionType/	std::map<collisionType/' src/Isotope.bak > src/Isotope.h
rm src/Isotope.bak

# Compile and link code
g++ src/arraycreator.h -c $CFLAGS
g++ src/interpolate.h -c $CFLAGS
g++ src/integrate.h -c $CFLAGS
g++ src/log.cpp -c $CFLAGS
g++ src/xsreader.cpp -c $CFLAGS	
g++ src/Isotope.cpp -c $CFLAGS
g++ src/Material.cpp -c $CFLAGS
g++ src/Neutron.cpp -c $CFLAGS
g++ src/Tally.cpp -c $CFLAGS
g++ src/Fissioner.cpp -c $CFLAGS
g++ src/Region.cpp -c $CFLAGS
g++ src/Geometry.cpp -c $CFLAGS
g++ python/pinspec/Geometry_wrap.cxx -c $CFLAGS
g++ $LFLAGS Isotope.o xsreader.o log.o Material.o Tally.o Neutron.o Region.o Fissioner.o Geometry.o Geometry_wrap.o -o python/pinspec/_pinspec.so

# Cleanup
mv *.o objs/
