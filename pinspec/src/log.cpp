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
char separator_char = '-';
char header_char = '*';
char title_char = '*';
int line_length = 67;           /* Must account for log level prefix */


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


void setSeparatorCharacter(char c) {
    separator_char = c;
}


void setHeaderCharacter(char c) {
    header_char = c;
}


void setTitleCharacter(char c) {
    title_char = c;
}


void setLineLength(int length) {
    line_length = length;
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
    case SEPARATOR:
	    log_printf(INFO, "Logging level set to SEPARATOR");
        break;
    case HEADER:
	    log_printf(INFO, "Logging level set to HEADER");
        break;
    case TITLE:
	    log_printf(INFO, "Logging level set to TITLE");
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
    else if (strcmp("SEPARATOR", newlevel) == 0) {
	    log_level = HEADER;
	    log_printf(INFO, "Logging level set to SEPARATOR");
    }
    else if (strcmp("HEADER", newlevel) == 0) {
	    log_level = HEADER;
	    log_printf(INFO, "Logging level set to HEADER");
    }
    else if (strcmp("TITLE", newlevel) == 0) {
	    log_level = TITLE;
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

    char message[512];
    std::string msg_string;

    if (level >= log_level) {
    	va_list args;

        va_start(args, format);
        vsprintf(message, format, args);
        va_end(args);

    	/* Append the log level to the message */
    	switch (level) {
	        case (DEBUG):
            {
                std::string msg = std::string(message);
                std::string level_prefix = "[  DEBUG  ]  ";

                /* If message is too long for a line, split into many lines */
                if (int(msg.length()) > line_length)
                    msg_string = createMultilineMsg(level_prefix, msg);

                /* Puts message on single line */
                else
                    msg_string = level_prefix + msg + "\n";

	            break;
            }
	        case (INFO):
            {
                std::string msg = std::string(message);
                std::string level_prefix = "[  INFO   ]  ";

                /* If message is too long for a line, split into many lines */
                if (int(msg.length()) > line_length)
                    msg_string = createMultilineMsg(level_prefix, msg);

                /* Puts message on single line */
                else
                    msg_string = level_prefix + msg + "\n";

	            break;
            }
	        case (NORMAL):
            {
                std::string msg = std::string(message);
                std::string level_prefix = "[  NORMAL ]  ";

                /* If message is too long for a line, split into many lines */
                if (int(msg.length()) > line_length)
                    msg_string = createMultilineMsg(level_prefix, msg);

                /* Puts message on single line */
                else
                    msg_string = level_prefix + msg + "\n";

	            break;
            }
	        case (SEPARATOR):
            {
                std::string pad = std::string(line_length, separator_char);
                std::string prefix = std::string("[SEPARATOR]  ");
                std::stringstream ss;
                ss << prefix << pad << "\n";
                msg_string = ss.str();
	            break;
            }
	        case (HEADER):
            {
                int size = strlen(message);
                int halfpad = (line_length - 4 - size) / 2;
                std::string pad1 = std::string(halfpad, header_char);
                std::string pad2 = std::string(halfpad + 
                                    (line_length - 4 - size) % 2, header_char);
                std::string prefix = std::string("[  HEADER ]  ");
                std::stringstream ss;
                ss << prefix << pad1 << "  " << message << "  " << pad2 << "\n";
                msg_string = ss.str();
	            break;
            }
	        case (TITLE):
            {
                int size = strlen(message);
                int halfpad = (line_length - size) / 2;
                std::string pad = std::string(halfpad, ' ');
                std::string prefix = std::string("[  TITLE  ]  ");
                std::stringstream ss;
                ss << prefix << std::string(line_length, title_char) << "\n";
                ss << prefix << pad << message << pad << "\n";
                ss << prefix << std::string(line_length, title_char) << "\n";
                msg_string = ss.str();
	            break;
            }
	        case (WARNING):
            {
                std::string msg = std::string(message);
                std::string level_prefix = "[ WARNING ]  ";

                /* If message is too long for a line, split into many lines */
                if (int(msg.length()) > line_length)
                    msg_string = createMultilineMsg(level_prefix, msg);

                /* Puts message on single line */
                else
                    msg_string = level_prefix + msg + "\n";

	            break;
            }
	        case (CRITICAL):
            {
                std::string msg = std::string(message);
                std::string level_prefix = "[ CRITICAL]  ";

                /* If message is too long for a line, split into many lines */
                if (int(msg.length()) > line_length)
                    msg_string = createMultilineMsg(level_prefix, msg);

                /* Puts message on single line */
                else
                    msg_string = level_prefix + msg + "\n";

	            break;
            }
	        case (RESULT):
                msg_string = std::string("[  RESULT ]  ") + message + "\n";
	            break;
	        case (ERROR):
	            va_start(args, format);
	            vsprintf(message, format, args);
 	            va_end(args);
                set_err(message);
                throw std::runtime_error(message);
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


std::string createMultilineMsg(std::string level, std::string message) {

    int size = message.length();

    std::string substring;
    int start = 0;
    int end = line_length;

    std::string msg_string;

    /* Loop over msg creating substrings for each line */
    while (end < size + line_length) {

        /* Append log level to the beginning of each line */
        msg_string += level;

        /* Begin multiline messages with ellipsis */
        if (start != 0)
            msg_string += "... ";

        /* Find the current full length substring for line*/
        substring = message.substr(start, line_length);

        /* Truncate substring to last complete word */
        if (end < size-1) {
            int endspace = substring.find_last_of(" ");
            if (message.at(endspace+1) != ' ' && 
                           endspace != int(std::string::npos)) {
                end -= line_length - endspace;
                substring = message.substr(start, end-start);
            }
        }

        /* concatenate substring to output message */
        msg_string += substring + "\n";

        /* Reduce line length to account for ellipsis prefix */
        if (start == 0)
            line_length -= 4;

        /* Update substring indices */
        start = end + 1;
        end += line_length + 1;
    }

    /* Reset line length */
    line_length += 4;

    return msg_string;
}


/**
 * Return the log_level
 * @return log_level
 */
int get_loglevel(){
	return log_level;
}

