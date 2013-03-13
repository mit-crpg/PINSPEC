/*
 * xsreader.h
 *
 *  Created on: Mar 20, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
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
#endif

#include "log.h"

static char* _xs_directory = (char*)"../xs-lib/";

void setXSLibDirectory(char* xs_directory);
char* getXSLibDirectory();
int parseCrossSections(const char* file, float* energies, float* xs_values);
int getNumCrossSectionDataPoints(const char* file);


#endif /* XSREADER_H_ */
