/*
 * log.h
 *
 *  Created on: Jan 22, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 *  Level-based logging module
 */

#ifndef LOG_H_
#define LOG_H_

#ifdef __cplusplus
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stdexcept>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

extern void set_err(const char *msg);

typedef enum logLevels {
	DEBUG,
	INFO,
	NORMAL,
    TITLE,
	WARNING,
	CRITICAL,
	RESULT,
	ERROR
} logLevel;

void setOutputDirectory(char* directory);
const char* getOutputDirectory();
void setLogfileName(char* filename);
void log_setlevel(logLevel newlevel);
void log_setlevel(const char* newlevel);
int get_loglevel();

void log_printf(logLevel level, const char *format, ...);

#ifndef LOG_C
	extern logLevel log_level;
    extern std::string logfile_name;
    extern std::string output_directory;
    extern bool logging;
#endif

#endif


#endif /* LOG_H_ */
