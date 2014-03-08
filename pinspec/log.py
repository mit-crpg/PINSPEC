## 
# @file log.py
# @package pinspec.log 
# @brief Utility functions for writing log messages to the screen.
#
# @details This module includes a set of wrapper functions for
#          the logging routines provided by PINSPEC's C++ source 
#          code. These Python methods provide an interface for 
#          creating formatted log messages using level-based loggin
#          and to print them to the screen as well as a logfile.
# @author Samuel Shaner
# @date March 15, 2013


import numpy
import pinspec


##
# @brief Function to print a log message to the screen
# @details This method is a wrapper to the log_printf C++ routine. It
#          allows for formatted messages to be printed to the screen 
#          in a similar fashion to the C/C++ printf method, but with
#          additional formatting provided by the PINSPEC logging utilities.
#          An example of how this might be used in a PINSPEC Python script
#          is as follows:
#
# @code
#          value1 = 25
#          value2 = 26.0
#          log.py_printf('NORMAL', 'My name is Will and I am %d going on'\
#                               ' %f years of age', value1, value2)
# @endcode
#
# @param level the logging level for this message
# @param my_str the string to print to the screen
# @param *args a variable length list of values for the message string
def py_printf(level, my_str, *args):
    if level == 'DEBUG':
        pinspec.log_printf(pinspec.DEBUG, my_str % args)
    elif level == 'INFO':
        pinspec.log_printf(pinspec.INFO, my_str % args)
    elif level == 'NORMAL':
        pinspec.log_printf(pinspec.NORMAL, my_str % args)
    elif level == 'SEPARATOR':
        pinspec.log_printf(pinspec.SEPARATOR, my_str % args)
    elif level == 'HEADER':
        pinspec.log_printf(pinspec.HEADER, my_str % args)
    elif level == 'TITLE':
        pinspec.log_printf(pinspec.TITLE, my_str % args)
    elif level == 'WARNING':
        pinspec.log_printf(pinspec.WARNING, my_str % args)
    elif level == 'CRITICAL':
        pinspec.log_printf(pinspec.CRITICAL, my_str % args)
    elif level == 'RESULT':
        pinspec.log_printf(pinspec.RESULT, my_str % args)
    elif level == 'UNITTEST':
        pinspec.log_printf(pinspec.UNITTEST, my_str % args)
    elif level == 'ERROR':
        pinspec.log_printf(pinspec.ERROR, my_str % args)


##
# @brief Assigns the lowest level logging message.
# @details Sets the lowest level logging message to print to the screen.
#          This controls the lowest level for both logging messages in the
#          C++ source code as well as the user's PINSPEC Python input file.
#          This function would be called at the beginning of the input file
#          as follows:
#
# @code
#          log.py_setlevel('INFO')
# @endcode
#
# @param level the minimum logging level ('DEBUG', 'INFO', etc)
def py_set_log_level(level):

    if level == 'DEBUG':
        pinspec.set_log_level('DEBUG')
    elif level == 'INFO':
        pinspec.set_log_level('INFO')
    elif level == 'NORMAL':
        pinspec.set_log_level('NORMAL')
    elif level == 'SEPARATOR':
        pinspec.set_log_level('SEPARATOR')
    elif level == 'HEADER':
        pinspec.set_log_level('HEADER')
    elif level == 'TITLE':
        pinspec.set_log_level('TITLE')
    elif level == 'WARNING':
        pinspec.set_log_level('WARNING')
    elif level == 'CRITICAL':
        pinspec.set_log_level('CRITICAL')
    elif level == 'RESULT':
        pinspec.set_log_level('RESULT')
    elif level == 'UNITTEST':
        pinspec.set_log_level('UNITTEST')
    elif level == 'ERROR':
        pinspec.set_log_level('ERROR')
    else:
        py_printf('Cannot set log level to unsupported log level %s', str(level))


