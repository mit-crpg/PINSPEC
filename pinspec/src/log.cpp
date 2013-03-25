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
std::string logfile_name = "";
std::string output_directory = ".";
bool logging = false;


void setOutputDirectory(char* directory) {

    output_directory = std::string(directory);

    /* Check to see if directory exists - if not, create it */
    struct stat st;
    if (!stat(directory, &st) == 0) {
		mkdir(directory, S_IRWXU);
        mkdir((output_directory+"/log").c_str(), S_IRWXU);
    }

    logfile_name = output_directory + "/" + logfile_name;

    return;
}


const char* getOutputDirectory() {
    return output_directory.c_str();
}


void setLogfileName(char* filename) {
    logfile_name = output_directory + "/" + filename;
}


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
    case TITLE:
	    log_printf(INFO, "Logging level set to TITLE");
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
    else if (strcmp("TITLE", newlevel) == 0) {
	    log_level = NORMAL;
	    log_printf(INFO, "Logging level set to TITLE");
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
	        case (TITLE):
            {
                int size = strlen(msg);
                int halfpad = (67 - size) / 2;
                std::string pad = std::string(halfpad, ' ');
                std::string prefix = std::string("[  TITLE  ]  ");
                std::stringstream ss;
                ss << prefix << std::string(67, '*') << "\n";
                ss << prefix << pad << msg << pad << "\n";
                ss << prefix << std::string(67, '*') << "\n";
                msg_string = ss.str();
	            break;
            }
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


        /* If this is our first time logging, add a header with date, time */
        if (!logging) {

            /* If output directory was not defined by user, then log file is
             * written to a "log" subdirectory. Create it if it doesn't exist */
            if (output_directory.compare(".") == 0) {
                struct stat st;
                if (!stat("log", &st) == 0)
                    mkdir("log", S_IRWXU);
            }

            /* Write the message to the output file */
            std::ofstream logfile;
            logfile.open (logfile_name.c_str(), std::ios::app); 

            /* Append date, time to the top of log output file */
            time_t rawtime;
            struct tm * timeinfo;
            time (&rawtime);
            timeinfo = localtime (&rawtime);
            logfile << "Current local time and date: " << asctime(timeinfo);
            logging = true;

            logfile.close();
        }

        /* Write the log message to the logfile */
        std::ofstream logfile;
        logfile.open (logfile_name.c_str(), std::ios::app); 
        logfile << msg_string;
        logfile.close();

        /* Write the log message to the shell */
        std::cout << msg_string;
    }
}


/**
 * Return the log_level
 * @return log_level
 */
int get_loglevel(){
	return log_level;
}

