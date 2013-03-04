#!/bin/bash -ex

CFLAGS="-O3 -march=native -fPIC `pkg-config python-2.7 --cflags`"
LFLAGS="-shared `pkg-config python-2.7 --libs` -fPI"

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

g++ src/arraycreator.h -c $CFLAGS
g++ src/interpolate.h -c $CFLAGS
g++ src/integrate.h -c $CFLAGS
g++ src/log.cpp -c $CFLAGS
g++ src/xsreader.cpp -c $CFLAGS
g++ src/Isotope.cpp -c $CFLAGS
g++ src/Material.cpp -c $CFLAGS
g++ src/Neutron.cpp -c $CFLAGS
g++ src/Tally.cpp -c $CFLAGS
g++ src/Region.cpp -c $CFLAGS
g++ SWIG/Region_wrap.cxx -c $CFLAGS
g++ $LFLAGS Isotope.o xsreader.o log.o Material.o Tally.o Neutron.o Region.o Region_wrap.o -o SWIG/_region.so
