/**
 * @file xsreader.h
 * @brief Utility functions to parse in ENDF cross-sections from text files
 *        into data arrays.
 * @author William Boyd (wboyd@mit.edu)
 * @date March 20, 2012
 *
 */

#ifndef XSREADER_H_
#define XSREADER_H_

#ifdef __cplusplus
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include "log.h"
#endif

void setXSLibDirectory(const char* directory);
const char* getXSLibDirectory();
int restoreXSLibrary();
int parseCrossSections(const char* file, float* energies, float* xs_values);
int getNumCrossSectionDataPoints(const char* filename);


#endif /* XSREADER_H_ */
