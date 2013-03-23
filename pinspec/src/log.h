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
#include <iostream>
#include <string.h>
#include <stdexcept>

extern void set_err(const char *msg);

typedef enum logLevels {
	DEBUG,
	INFO,
	NORMAL,
	WARNING,
	CRITICAL,
	RESULT,
	ERROR
} logLevel;

void log_setlevel(logLevel newlevel);
void log_setlevel(const char* newlevel);
int get_loglevel();

void log_printf(logLevel level, const char *format, ...);

#ifndef LOG_C
	extern logLevel log_level;
#endif

#endif


#endif /* LOG_H_ */