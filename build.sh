#!/bin/bash -ex

os_type=`uname`;

if [[ $os_type == *Linux* ]]; then
    echo "Using linux flags"
    CFLAGS="-O3 -march=native -fPIC -fopenmp `pkg-config python-2.7 --cflags` -I/usr/lib64/python2.7/site-packages/numpy/core/include/"
    LFLAGS="-shared `pkg-config python-2.7 --libs` -fPI -fopenmp -lgomp -I/usr/lib64/python2.7/site-packages/numpy/core/include/"
elif [[ $os_type == *Darwin* ]]; then
    echo "Using Mac OSX flags"
    CFLAGS="-O3 -march=native -fPIC -I/usr/include/python2.7"
    LFLAGS="-bundle `python-config --ldflags` -fPI -lgomp -L/opt/local/lib/gcc47"    
fi

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


swig -python -c++ -threads pinspec/Geometry.i

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
g++ pinspec/Geometry_wrap.cxx -c $CFLAGS
g++ $LFLAGS Isotope.o xsreader.o log.o Material.o Tally.o Neutron.o Region.o Fissioner.o Geometry.o Geometry_wrap.o -o pinspec/_pinspec.so

# Cleanup
mv *.o objs/
