/*
 * log.cpp
 *
 *  Created on: Jan 22, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#define LOG_C

#include "log.h"


/* Default logging level is the lowest (most verbose) level */
logLevel log_level = NORMAL;


/**
 * Set the minimum level which will be printed to the console
 * @param newlevel the logging level
 */
void log_setlevel(logLevel newlevel) {
    log_level = newlevel;

    switch (newlevel) {
    case DEBUG:
	log_printf(INFO, "Logging level set to DEBUG");
	break;
    case INFO:
	log_printf(INFO, "Logging level set to INFO");
	break;
    case NORMAL:
	log_printf(INFO, "Logging level set to NORMAL");
	break;
    case WARNING:
	log_printf(INFO, "Logging level set to WARNING");
	break;
    case CRITICAL:
	log_printf(INFO, "Logging level set to CRITICAL");
	break;
    case RESULT:
	log_printf(INFO, "Logging level set to RESULT");
	break;
    case ERROR:
	log_printf(INFO, "Logging level set to ERROR");
	break;
    }
}


/**
 * Set the logging level from a character string which must correspond
 * to one of the logLevel enum types (NORMAL, INFO, CRITICAL, WARNING,
 * ERROR, DEBUG, RESULT)
 * @param newlevel a character string loglevel
 */
void log_setlevel(const char* newlevel) {

    if (strcmp("DEBUG", newlevel) == 0) {
	log_level = DEBUG;
	log_printf(INFO, "Logging level set to DEBUG");
    }
    else if (strcmp("INFO", newlevel) == 0) {
	log_level = INFO;
	log_printf(INFO, "Logging level set to INFO");
    }
    else if (strcmp("NORMAL", newlevel) == 0) {
	log_level = NORMAL;
	log_printf(INFO, "Logging level set to NORMAL");
    }
    else if (strcmp("WARNING", newlevel) == 0) {
	log_level = WARNING;
	log_printf(INFO, "Logging level set to WARNING");
    }
    else if (strcmp("CRITICAL", newlevel) == 0) {
	log_level = CRITICAL;
	log_printf(INFO, "Logging level set to CRITICAL");
    }
    else if (strcmp("RESULT", newlevel) == 0) {
	log_level = RESULT;
	log_printf(INFO, "Logging level set to RESULT");
    }
    else if (strcmp("ERROR", newlevel) == 0) {
	log_level = ERROR;
	log_printf(INFO, "Logging level set to ERROR");
    }

    return;
}


/**
 * Print a formatted message to the console. If logging level is
 * ERROR, this function will end the program
 * @param level the logging level for this message
 * @param *format variable list of C++ formatted i/o
 */
void log_printf(logLevel level, const char *format, ...) {

    char msg[512];
    std::string msg_string;

    if (level >= log_level) {
    	va_list args;

        va_start(args, format);
        vsprintf(msg, format, args);
        va_end(args);

    	/* Append the log level to the message */
    	switch (level) {
	        case (DEBUG):
                msg_string = std::string("[  DEBUG  ]  ") + msg + "\n";
	            break;
	        case (INFO):
                msg_string = std::string("[  INFO   ]  ") + msg + "\n";
	            break;
	        case (NORMAL):
                msg_string = std::string("[  NORMAL ]  ") + msg + "\n";
	            break;
	        case (WARNING):
                msg_string = std::string("[ WARNING ]  ") + msg + "\n";
	            break;
	        case (CRITICAL):
                msg_string = std::string("[ CRITICAL]  ") + msg + "\n";
	            break;
	        case (RESULT):
                msg_string = std::string("[  RESULT ]  ") + msg + "\n";
	            break;
	        case (ERROR):
	            va_start(args, format);
	            vsprintf(msg, format, args);
 	            va_end(args);
                set_err(msg);
                throw std::runtime_error(msg);
	            break;
          }
    }

   std::cout << msg_string;

}


/**
 * Return the log_level
 * @return log_level
 */
int get_loglevel(){
	return log_level;
}

